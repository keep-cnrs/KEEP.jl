using Pkg
Pkg.activate("CT_keep")

using StaticArrays
using ForwardDiff: derivative
using OptimalControl, NLPModelsIpopt
using ComponentArrays: ComponentArray as CA
using Plots
using NonlinearSolve

using KEEP.PointMassPara: build_vbpara, lmt
using KEEP.PointMass4: dynamics
using KEEP.TorqueFunction: torque_function
using KEEP.LimitCycle: compute_limit_cycle
using KEEP: TAU0

const vbp0 = build_vbpara()
const syms = (:r, :I_eq, :torque_slope)

begin  # Helper Functions
    f(x, p) = begin
        T = promote_type(eltype(vbp0), eltype(p))
        vbp = CA(T.(vbp0); (syms .=> p)...)
        x_dyn = vcat(x, 0)
        return dynamics(x_dyn, vbp)[SOneTo(4)]
    end

    generated_power(dα, p) = begin
        T = promote_type(eltype(vbp0), eltype(p))
        vbp = CA(T.(vbp0); (syms .=> p)...)
        L, M, T = lmt(vbp)
        dα_adim = dα * T
        return dα_adim * torque_function(dα_adim, vbp) * M * L^2 * T^-3
    end

    function calc_vars(lc)
        tf = last(lc.t)
        p = lc.prob.p[syms]
        return [tf, p...]
    end

    function build_optim_ocp(lc; factor)
        ocp = @def begin
            vars ∈ R⁴, variable
            tf = vars[1]
            p = vars[2:4]
            t ∈ [0, tf], time
            X = (α, τ, dα, dτ) ∈ R⁴, state

            calc_vars(lc) ./ factor <= vars <= calc_vars(lc) .* factor

            τ(0) == TAU0
            α(tf) - α(0) == 0
            τ(tf) - τ(0) == 2π
            dα(tf) - dα(0) == 0
            dτ(tf) - dτ(0) == 0

            Ẋ(t) == f(X(t), p)
            ∫(generated_power(dα(t), p) / tf) → max
        end

        init = @init ocp begin
            vars := calc_vars(lc)
            X(t) := lc(t)[1:4]
        end
        return ocp, init
    end
end

# Define structures for Indirect Multiple Shooting
struct MultipleShooting{F,G,P}
    flow::F
    dynamics::G
    power_func::P
    N::Int
    nx::Int
    nv::Int
    tau0::Float64
end

struct StateInterpolator{F}
    flow::F
    N::Int
    tf::Float64
    v::Vector{Float64}
    nodes_x::Vector{Vector{Float64}}
    nodes_p::Vector{Vector{Float64}}
end

# Callable implementation for the interpolator
function (interp::StateInterpolator)(t)
    t_clamped = clamp(t, 0.0, interp.tf)
    dt = interp.tf / interp.N

    if t_clamped == interp.tf
        i = interp.N - 1
    else
        i = floor(Int, t_clamped / dt)
        i = clamp(i, 0, interp.N - 1)
    end

    t_start = i * dt
    x_node = interp.nodes_x[i+1]
    p_node = interp.nodes_p[i+1]

    if t_clamped == t_start
        return x_node
    end

    xf, _, _ = interp.flow(t_start, x_node, p_node, t_clamped, interp.v; augment=true)
    return xf
end

@views function shoot!(res, x0, p0, xf, pf, tf, λ, pλf)
    res[1] = x0[2] - TAU0
    res[2:5] = xf - x0 - [0, 2π, 0, 0]
    res[6] = pf' * f(xf, λ) - generated_power(xf[3], λ) / tf # Hf
    res[7:9] = (pf-p0)[[1, 3, 4]]
    res[10:12] = pλf
    return res
end

@views function shoot!(res, x0, p0, tf, λ)
    xf, pf, pvf = flow(0.0, x0, p0, tf, [tf; λ]; augment=true)
    pλf = pvf[2:4]
    return shoot!(res, x0, p0, xf, pf, tf, λ, pλf)
end

# Multiple shooting residual evaluation
function (ms::MultipleShooting)(res, Y)
    N = ms.N
    nx = ms.nx
    nv = ms.nv
    nw = 2 * nx + nv

    v = Y[N*nw+1:end]
    tf = v[1]
    dt = tf / N

    w0 = Y[1:nw]
    x_curr = w0[1:nx]
    p_curr = w0[nx+1:2*nx]
    pv_curr = w0[2*nx+1:2*nx+nv]

    x0 = x_curr
    p0 = p_curr
    pv0 = pv_curr

    for i in 0:N-2
        t_curr = i * dt
        t_next = (i + 1) * dt

        xf, pf, dpv = ms.flow(t_curr, x_curr, p_curr, t_next, v; augment=true)
        pv_next_pred = pv_curr + dpv

        idx_next = (i + 1) * nw
        w_next = Y[idx_next+1:idx_next+nw]
        x_next = w_next[1:nx]
        p_next = w_next[nx+1:2*nx]
        pv_next = w_next[2*nx+1:2*nx+nv]

        res[idx_next-nw+1:idx_next-nw+nx] .= x_next .- xf
        res[idx_next-nw+nx+1:idx_next-nw+2*nx] .= p_next .- pf
        res[idx_next-nw+2*nx+1:idx_next] .= pv_next .- pv_next_pred

        x_curr = x_next
        p_curr = p_next
        pv_curr = pv_next
    end

    t_curr = (N - 1) * dt
    t_next = tf
    xf, pf, dpv = ms.flow(t_curr, x_curr, p_curr, t_next, v; augment=true)
    pvf = pv_curr + dpv

    # The continuity matching equations filled up to (N - 1) * nw
    offset = (N - 1) * nw

    # 1. 12 Boundary conditions (cyclic constraints and Hamiltonian transversality)
    shoot!(view(res, offset+1:offset+12), x0, p0, xf, pf, tf, v[2:end], pvf[2:end])

    # 2. 4 Parameter adjoint initial conditions (pv0 == 0)
    res[offset+13:offset+12+nv] .= pv0

    return res
end

# Automatic size-inferring initialization
function init_multiple_shooting(flow, optim_sol, N)
    tf, λ... = variable(optim_sol)
    v = [tf; λ...]
    dt = tf / N

    x0_direct = state(optim_sol)(0)
    p0_direct = costate(optim_sol)(0)
    nx = length(x0_direct)
    nv = length(v)
    nw = 2 * nx + nv

    x_curr = copy(x0_direct)
    p_curr = copy(p0_direct)
    pv_curr = zeros(nv)

    Y_init = zeros(N * nw + nv)

    for i in 0:N-1
        idx = i * nw
        Y_init[idx+1:idx+nx] .= x_curr
        Y_init[idx+nx+1:idx+2*nx] .= p_curr
        Y_init[idx+2*nx+1:idx+nw] .= pv_curr

        if i < N - 1
            t_curr = i * dt
            t_next = (i + 1) * dt
            x_curr, p_curr, dpv = flow(t_curr, x_curr, p_curr, t_next, v; augment=true)
            pv_curr .+= dpv
        end
    end

    Y_init[N*nw+1:end] .= v
    return Y_init
end

# Reconstruct trajectory and generate a callable interpolator
function reconstruct_trajectory(ms::MultipleShooting, Y; num_points=200)
    N = ms.N
    nx = ms.nx
    nv = ms.nv
    nw = 2 * nx + nv

    v = Y[N*nw+1:end]
    tf = v[1]
    dt = tf / N

    nodes_x = [Y[i*nw+1:i*nw+nx] for i in 0:N-1]
    nodes_p = [Y[i*nw+nx+1:i*nw+2*nx] for i in 0:N-1]

    interpolator = StateInterpolator(ms.flow, N, tf, v, nodes_x, nodes_p)

    ts = range(0.0, tf, length=num_points)
    xs = [interpolator(t) for t in ts]

    return ts, reduce(hcat, xs), interpolator
end


# --- Execution and Setup ---

# Solve direct problem first
lc = compute_limit_cycle(vbp0; sense=+, save_everystep=true)
factor = 5
optim_ocp, optim_init = build_optim_ocp(lc; factor=factor)
oc_kwargs = (grid_size=30, backend=:generic)
optim_init_sol = solve(optim_ocp; oc_kwargs..., init=optim_init, max_iter=0)
optim_sol = solve(optim_ocp; oc_kwargs..., init=optim_init)

const flow_instance = Flow(optim_ocp)

# Infer dimensions automatically
nx_dim = length(state(optim_sol)(0))
nv_dim = length(variable(optim_sol))

# Setup indirect multiple shooting
N_segments = 8
ms_problem = MultipleShooting(flow_instance, f, generated_power, N_segments, nx_dim, nv_dim, TAU0)
Y_guess = init_multiple_shooting(flow_instance, optim_sol, N_segments)

# 1. Define the in-place residual matching the (res, u, p) signature
residual_f!(res, Y, p) = ms_problem(res, Y)

# 2. Define the NonlinearProblem
prob_indirect = NonlinearProblem(residual_f!, Y_guess)

# 3. Solve the problem using TrustRegion (the most stable workhorse for BVPs)
sol_indirect = solve(prob_indirect, TrustRegion(); show_trace=Val(true),)

# 4. Extract solution and reconstruct
if sol_indirect.retcode == ScalarSymbolic.Success # or check SciMLBase.successful_retcode(sol_indirect)
    Y_opt = sol_indirect.u
    ts, xs, state_interp = reconstruct_trajectory(ms_problem, Y_opt; num_points=100)

    # Access parameters
    tf_opt = Y_opt[12*N_segments+1]
    λ_opt = Y_opt[12*N_segments+2:end]
    println("Optimized tf: ", tf_opt)
    println("Optimized parameters: ", λ_opt)
else
    println("Solver failed to converge with retcode: ", sol_indirect.retcode)
end