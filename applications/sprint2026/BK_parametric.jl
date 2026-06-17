using Pkg
Pkg.activate("CT_keep")

import OrdinaryDiffEq as ODE
using BifurcationKit
using LinearAlgebra
using ComponentArrays: ComponentArray as CA
using StaticArrays
using Plots

using KEEP: TAU0
using KEEP.PointMass4: dynamics
using KEEP.PointMassPara: build_vbpara, build_para, lmt
using KEEP.Optimization: optimize
using KEEP.LimitCycle: shoot as lc_shoot

## Solve the optimization problem
const vbp0 = build_vbpara()
const p0 = build_para(vbp0)
const nt_p0 = (; zip(keys(p0), p0)...)
const syms = (:r, :I_eq, :torque_slope)

factor = 5
lb = p0[syms] ./ factor
ub = p0[syms] .* factor
solution, stats, model = optimize(p0, syms, lb, ub)

## Reconstruct the dense ODE solution
vbp = build_vbpara(CA(p0; solution.params...))
shooting = solution[1:4]
solution_sim = lc_shoot(shooting, vbp, save_everystep=true)

x_optimization = t -> solution_sim(t, idxs=1:4)
tf_optimization = solution_sim.t[end]

## BVP model
function F(u, params, t=0)
    α, τ, dα, dτ, tf = u
    u_dyn = SA[α, τ, dα, dτ, 0]
    _, _, ddα, ddτ, _ = dynamics(u_dyn, build_vbpara(params))
    dtf = 0
    return tf .* [dα, dτ, ddα, ddτ, dtf]
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

u0_bif = SA[x_optimization(0)..., tf_optimization]  # α, τ, dα, dτ, tf

const STATE_SIZE = length(u0_bif)
model = BifurcationKit.BVP.BVPModel(F, g; n=STATE_SIZE)

## Collocation
grid_size, degree = 30, 5
const disc = BifurcationKit.BVP.Collocation(Ntst=grid_size, m=degree, meshadapt=true)
bvp = BifurcationKit.BVP.discretize(model, disc)

params = nt_p0
# we could also do x0 = BK.BVP.generate_solution(bvp, my_guess_function)
x0 = BifurcationKit.BVP.generate_solution(bvp, s -> vcat(x_optimization(tf_optimization * s), tf_optimization))

prob = BifurcationKit.BVP.BVPBifProblem(bvp, x0, params, (@optic _.v_ref);
    jacobian=BifurcationKit.DenseAnalytical(),
    # record_from_solution,
    # plot_solution_col
)

# @assert false
optn = NewtonPar(tol=1e-10, verbose=true)

# J = BifurcationKit.jacobian(prob, prob.u0, prob.params)

sol = @time BifurcationKit.solve(prob, Newton(), optn);

# plot();
# plot_solution_col(sol.u, prob.params);

# Continuation
optc = ContinuationPar(
    p_min=0.0,
    p_max=25.0,
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


## Multiple shooting
function plot_solution_ms(x, p; kwargs...)
    sol = BifurcationKit._get_shooting_solution(bvp.cache, reshape(x, 5, disc.M), 1, @set params.v_ref = p)

    plot!(sol.t, sol.u[4, :]; ylabel="u(t)", title="Bratu Solution (p₁=)", kwargs...)
end


odeprob = ODE.ODEProblem(F, u0_bif, (0, 1), nt_p0)
model = BifurcationKit.BVP.BVPModel(odeprob, g; n=5)
# 4. Discretize using Collocation method
# Using 201 points for better accuracy
disc2 = BifurcationKit.BVP.Shooting(10, ODE.Vern9(), true)
bvp = BifurcationKit.BVP.discretize(model, disc2; abstol=1e-12, reltol=1e-10)
#https://github.com/bifurcationkit/MultiParamContinuation.jl

x0 = BifurcationKit.BVP.generate_solution(bvp, s -> vcat(x_optimization(tf_optimization * s), tf_optimization))

prob = BifurcationKit.BVP.BVPBifProblem(bvp, x0, params, (@optic _.v_ref);
)

optn = NewtonPar(tol=1e-10, verbose=true)

@error "un objet est mal initialisé lors du solve, exécuter plusieurs fois donne différents résultats"
sol = @time BifurcationKit.solve(prob, Newton(), optn);

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
