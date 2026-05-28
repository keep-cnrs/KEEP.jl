using Pkg
Pkg.activate("CT_keep")

import OrdinaryDiffEq as ODE
using BifurcationKit
using LinearAlgebra
using ComponentArrays: ComponentArray as CA
using StaticArrays
using Plots

#############

using KEEP: TAU0
using KEEP.PointMass4: dynamics
using KEEP.PointMassPara: build_vbpara, build_para, lmt
using KEEP.Optimization: optimize
# function optimize(p, optim_para_syms, optim_para_lower, optim_para_upper; sense=+, tol=10DEFAULT_TOLERANCE, max_wall_time=120., linear_solver="mumps", initial_guess=(;))

const vbp0 = build_vbpara()
const syms = (:r, :I_eq, :torque_slope)

L, M, T = lmt(vbp0)
p0 = build_para(vbp0)
units = p0[syms] ./ vbp0[syms]
factor = 5
lb = p0[syms] ./ factor
ub = p0[syms] .* factor
solution, stats, model = optimize(p0, syms, lb, ub)

using KEEP.LimitCycle: shoot as lc_shoot

vbp = build_vbpara(CA(p0; solution.params...))
shooting = solution[1:4]
solution_sim = lc_shoot(shooting, vbp, save_everystep=true)

x_direct = t -> solution_sim(t, idxs=1:4)
tf_direct = solution_sim.t[end]

##################

function record_from_solution(x, p; k...)
    # u = BifurcationKit.get_time_slices(x[1:end], STATE_SIZE, disc.m, disc.Ntst)
    return (max_u=norm(x, 2), s=sum(x))
end

function plot_solution(x, p; kwargs...)
    u = BifurcationKit.get_time_slices(x[1:end], STATE_SIZE, disc.m, disc.Ntst)
    plot!(u[4, :]; ylabel="u(t)", title="KEEP Solution (p₁=)", legend=false, marker=:d, markersize=1, kwargs...)
end


function F(u, params, t=0)
    α, τ, dα, dτ, tf = u
    u_dyn = SA[α, τ, dα, dτ, 0]
    # vbp = @set nt_vbp0[keys(p)] = values(p)
    _, _, ddα, ddτ, _ = dynamics(u_dyn, params)
    dtf = 0
    return tf .* SA[dα, dτ, ddα, ddτ, dtf]
end

function g(u0, uT, p)
    return SA[
        uT[1] - u0[1]       # α loop
        uT[3] - u0[3]       # dα loop
        uT[4] - u0[4]       # dτ loop
        u0[2] - TAU0        # τ init
        uT[2] - TAU0 - 2π   # τ end
    ]
end

nt_vbp0 = (; zip(keys(vbp0), vbp0)...)  # Named Tuple with same parameters
u0_bif = SA[x_direct(0)..., tf_direct]  # α, τ, dα, dτ, tf

## Collocation
const STATE_SIZE = length(u0_bif)
model = BifurcationKit.BVP.BVPModel(F, g; n=STATE_SIZE)
# Rename to grid_size and degree
grid_size, degree = 30, 5
const disc = BifurcationKit.BVP.Collocation(Ntst=grid_size, m=degree, meshadapt=true)
bvp = BifurcationKit.BVP.discretize(model, disc)

params = nt_vbp0
# we could also do x0 = BK.BVP.generate_solution(bvp, my_guess_function)
x0 = BifurcationKit.BVP.generate_solution(bvp, s -> vcat(x_direct(tf_direct * s), tf_direct))

prob = BifurcationKit.BVP.BVPBifProblem(bvp, x0, params, (@optic _.v_ref);
    jacobian=BifurcationKit.DenseAnalytical(),
    record_from_solution,
    plot_solution
)

# @assert false
optn = NewtonPar(tol=1e-10, verbose=true)

# J = BifurcationKit.jacobian(prob, prob.u0, prob.params)

sol = @time BifurcationKit.solve(prob, Newton(), optn);

plot();
plot_solution(sol.u, prob.params);

optc = ContinuationPar(
    p_min=0.1,
    p_max=50.05,
    dsmax=0.1,
    ds=0.01,
    detect_bifurcation=0,
    # detect_fold = false,
    newton_options=optn,
    max_steps=100,
    nev=20,
    n_inversion=6
)

br = continuation(prob, PALC(), optc;
    plot=true,
    verbosity=1,
    normC=norminf,
    bothside=true,
)
plot(br)

# plot(br)



## Multiple shooting

odeprob = ODE.ODEProblem(F, u0_bif, (0, 1), nt_vbp0)
model = BifurcationKit.BVP.BVPModel(odeprob, g; n=2)
# 4. Discretize using Collocation method
# Using 201 points for better accuracy
const disc = BifurcationKit.BVP.Shooting(10, ODE.Vern9(), true)
bvp = BifurcationKit.BVP.discretize(model, disc; abstol=1e-12, reltol=1e-10)

nothing

#######
## Sparse test
#######

using DifferentiationInterface
using ADTypes
import ForwardDiff
using SparseArrays

# 1. Define the residual
f_bvp_init(u) = BifurcationKit.residual(prob, u, params)

# 2. Compute sparsity pattern
@info "Computing initial dense Jacobian for BVP sparsity pattern..."
J_dense = ForwardDiff.jacobian(f_bvp_init, x0)
sparsity_pattern = sparse(J_dense)
droptol!(sparsity_pattern, 1e-12)

# 3. Setup backend with qualified KnownJacobianSparsityDetector
backend = AutoSparse(
    AutoForwardDiff();
    sparsity_detector=ADTypes.KnownJacobianSparsityDetector(sparsity_pattern)
)

# 4. Prepare the cache
extras = prepare_jacobian(f_bvp_init, backend, x0; strict=Val(false))

# 5. Define the sparse Jacobian
function J_sparse_bvp(u, p)
    f_current(v) = BifurcationKit.residual(prob, v, p)
    return jacobian(f_current, extras, backend, u)
end

# 6. Rebuild the BVPBifProblem using the sparse Jacobian function
prob_sparse = BifurcationKit.BVP.BVPBifProblem(bvp, x0, params, (@optic _.v_ref);
    jacobian=J_sparse_bvp,
    record_from_solution,
    plot_solution
)

# function BifurcationKit.jacobian(prob::BifurcationKit.BVP.BVPBifProblem, u, p)
#     return prob.jacobian(u, p)
# end
function BifurcationKit.jacobian(prob::BifurcationKit.BVP.BVPBifProblem{Tbvp,<:Function}, u, p) where Tbvp
    return prob.jacobian(u, p)
    @info "Jacobian metadata:" typeof(J) size(J) nnz(J) density = (nnz(J) / length(J))
end

# 7. Define Newton parameters
optn_sparse = NewtonPar(
    tol=1e-10,
    verbose=false,
    linsolver=DefaultLS()
)

# Test sparse Newton solve
@info "Solving BVP with sparse Jacobian..."

@benchmark sol_sparse = BifurcationKit.solve(prob_sparse, Newton(), optn_sparse)  # 170 ms
@benchmark sol = BifurcationKit.solve(prob, Newton(), optn)  # 30 ms
#https://github.com/bifurcationkit/MultiParamContinuation.jl
