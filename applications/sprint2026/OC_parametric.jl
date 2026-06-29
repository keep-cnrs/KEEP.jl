## Setup
using Pkg
Pkg.activate("CT_keep")

using StaticArrays
using ForwardDiff: ForwardDiff
using ForwardDiff: derivative
using OptimalControl, NLPModelsIpopt
using ComponentArrays: ComponentArray as CA
using Plots
using LinearAlgebra
using OrdinaryDiffEq
using NonlinearSolve
using Test
using Printf

using KEEP.PointMassPara: build_vbpara, lmt, build_para
using KEEP.PointMass4: dynamics
using KEEP.TorqueFunction: torque_function
using KEEP.LimitCycle: compute_limit_cycle, shoot as lc_shoot
using KEEP.Optimization: optimize
using KEEP: TAU0

begin
    const vbp0 = build_vbpara()
    const nt_vbp0 = (; zip(keys(vbp0), vbp0)...)
    const syms = (:r, :I_eq, :torque_slope)

    function f(x, p, t=0)
        x_dyn = vcat(x, 0)
        p_tuple = ntuple(i -> p[i], Val(length(syms)))
        vbp = merge(nt_vbp0, NamedTuple{syms}(p_tuple))
        return dynamics(x_dyn, vbp)[SOneTo(4)]
    end

    function generated_power(dα, p)
        p_tuple = ntuple(i -> p[i], Val(length(syms)))
        vbp = merge(nt_vbp0, NamedTuple{syms}(p_tuple))
        L, M, T = lmt(vbp)
        dα_adim = dα * T
        return dα_adim * torque_function(dα_adim, vbp) * M * L^2 * T^-3
    end

    f(rand(4), [15, 50, 40])
    generated_power(rand(), [15, 50, 40])

    function calc_vars(lc)
        tf = last(lc.t)
        p = lc.prob.p[syms]
        return [tf, p...]
    end

    function build_lc_ocp(; tf_init=2.8)
        ocp_lc = @def begin
            vars ∈ R², variable  # tf, dummy
            tf = vars[1]
            t ∈ [0, tf], time
            X = (α, τ, dα, dτ) ∈ R⁴, state
            τ(0) == TAU0
            α(tf) - α(0) == 0
            τ(tf) - τ(0) == 2π
            dα(tf) - dα(0) == 0
            dτ(tf) - dτ(0) == 0
            Ẋ(t) == f(X(t), collect(vbp0[syms]))
            vars[2]^2 → min
        end
        estimate_τ(t, tf) = 2π * t / tf + TAU0
        estimate_α(t, tf) = sin(estimate_τ(t, tf) - π / 4)
        estimate_dα(t, tf) = derivative(t_ -> estimate_α(t_, tf), t)
        estimate_dτ(t, tf) = derivative(t_ -> estimate_τ(t_, tf), t)
        init_lc = @init ocp_lc begin
            vars := [tf_init, 0]
            α(t) := estimate_α(t, tf_init)
            τ(t) := estimate_τ(t, tf_init)
            dα(t) := estimate_dα(t, tf_init)
            dτ(t) := estimate_dτ(t, tf_init)
        end
        return ocp_lc, init_lc
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

    function solve_ocp(ocp, init; kwargs...)
        sol_init = solve(ocp; kwargs..., init=init, max_iter=0)
        sol = solve(ocp; kwargs..., init=init)
        return sol, sol_init
    end

    function make_var_plot(sol, sol_init, factor)
        x_ticks = 1:(length(syms) + 1)
        ratios = variable(sol) ./ variable(sol_init)
        x_labels = (x_ticks, (:tf, syms...))

        plot(;
            yscale=:log10,
            xticks=x_labels,
            ylabel="Ratio (Optimized / Initial)",
            title="Parameter Optimization Results",
            framestyle=:box,
            legend=:outertopright,
        )

        hline!(
            [1 / factor]; fillrange=[factor], color=:grey, alpha=0.15, label="Allowed Range"
        )
        hline!([1]; c=:black, ls=:dash, lw=1.5, label="Baseline (1.0)")

        for i in x_ticks
            plot!([i, i], [1.0, ratios[i]]; c=:grey, lw=1.0, label="")
        end

        scatter!(
            x_ticks,
            ones(length(x_ticks));
            m=:circle,
            mc=:white,
            msc=:black,
            ms=6,
            label="Initial",
        )
        scatter!(x_ticks, ratios; m=:circle, mc=:crimson, ms=7, label="Optimized")

        ylims!(1 / factor / 1.2, factor * 1.2)
        return plot!()
    end

    function make_plots(sol_init, sol, factor; legend=true)
        plot(sol_init; label="Initial")
        state_plot = plot!(sol; label="Optimized")
        plot!(; legend)
        var_plot = make_var_plot(sol, sol_init, factor)
        return state_plot, var_plot
    end
end

oc_kwargs = (grid_size=30, backend=:generic)

## Computing a limit cycle (feasability)
# Code from KEEP
lc_keep = compute_limit_cycle(vbp0; sense=(+), save_everystep=true);

# With ControlToolbox
lc_ocp, lc_init = build_lc_ocp();
lc_sol, lc_init_sol = solve_ocp(lc_ocp, lc_init; oc_kwargs...);

vec_keep = [lc_keep.u[1][1:4]; lc_keep.t[end]]
vec_CT = [state(lc_sol)(0); time_grid(lc_sol)[end]]
vec_sum = vec_keep + vec_CT
vec_diff = vec_keep - vec_CT
diff = sum(abs, (2vec_diff ./ vec_sum))

println("############")
@printf "Finding a limit cycle -- norm of relative difference = %.2e\n" diff
println("############")

## Solving the parametric optimization problem
# divide/multiply initial guess parameters by factors to make lower/upper bound
factor = 5

# With KEEP

# dimensionless -> dimension parameters
p0 = build_para(vbp0)

# Parameter box for optimization
lb = p0[syms] ./ factor
ub = p0[syms] .* factor

# Optimize
solution, stats, model = optimize(p0, syms, lb, ub)

# I use lc_keep because it has all the necessary elements for initializing
optim_ocp, optim_init = build_optim_ocp(lc_keep; factor=factor);

optim_sol, optim_init_sol = solve_ocp(optim_ocp, optim_init; oc_kwargs...);

# blue: initial, green: optimized
display.(make_plots(optim_init_sol, optim_sol, factor; legend=false))

param_units = p0[syms] ./ vbp0[syms]

vec_keep = [solution.tf; solution.params ./ param_units]
vec_CT = variable(optim_sol)
vec_sum = vec_keep + vec_CT
vec_diff = vec_keep - vec_CT
diff = sum(abs, (2vec_diff ./ vec_sum))

println("############")
@printf "Parametric optimization -- norm of relative difference = %.2e\n" diff
println("############")

throw("broken from this point on")

## Indirect single shooting
const flow = OptimalControl.Flow(optim_ocp; abstol=1e-12);

x_direct = state(optim_sol)
p_direct = costate(optim_sol)
tf_direct, λ_direct... = variable(optim_sol)

# https://control-toolbox.org/MagneticResonanceImaging.jl/stable/saturation.html
@views function shoot!(res, x0, p0, xf, pf, tf, λ, pλf)
    res[1] = x0[2] - TAU0
    res[2:5] = xf - x0 - [0, 2π, 0, 0]
    res[6] = pf' * f(xf, λ) - generated_power(xf[3], λ) / tf # Hf
    res[7:9] = (pf - p0)[[1, 3, 4]]
    res[10:12] = pλf
    return res
end

@views function shoot!(res, x0, p0, tf, λ)
    # λ is previous p
    xf, pf, pvf = flow(0.0, x0, p0, tf, [tf; λ]; augment=true)
    pλf = pvf[2:4]
    return shoot!(res, x0, p0, xf, pf, tf, λ, pλf)
end

# shoot!(res, y, _) = shoot!(res, y[SA[1, 2, 3, 4]], y[SA[5, 6, 7, 8]], y[9], y[SA[10, 11, 12]])  # for NL-solve
@views shoot!(res, y, _) = shoot!(res, y[1:4], y[5:8], y[9], y[10:12])  # for NL-solve
const ss_res = similar(y_guess)
shoot(y) = shoot!(ss_res, y, nothing)  # For debugging

x0_guess = x_direct(0)
p0_guess = p_direct(0)
tf_guess, λ_guess = tf_direct, λ_direct

y_guess = [x0_guess; p0_guess; tf_guess; λ_guess]
@time @profview shoot(y_guess)

## Single shooting demonstration
# The shooting problem is unstable
prob_single = NonlinearProblem(shoot!, y_guess)
sol_single = solve(prob_single, NewtonRaphson(); maxiters=5)
println("Single shooting retcode: ", sol_single.retcode)
println("Residual norm: ", norm(sol_single.resid))

using KEEP.PointMass4: integrate

## Monodromy matrix
"""
    flow_jacobian(u0, tf, vbp; idxs=1:length(u0))

Integrates the system with ForwardDiff.Dual numbers seeded at `u0`.
Returns two functions and the time grid:
- `f(t)`: state vector at time t for indices in `idxs`.
- `J(t)`: Jacobian matrix (∂u_{idxs}(t)/∂u0_{idxs}) at time t.
- `t_grid`: time points from the ODE solution.
"""
function flow_jacobian(u0, tf, vbp; idxs=1:length(u0))
    N_total = length(u0)
    N = length(idxs)
    T = eltype(u0)
    tag = ForwardDiff.Tag(integrate, T)

    partial_of = zeros(Int, N_total)
    for (j, i) in enumerate(idxs)
        partial_of[i] = j
    end

    u0_dual = SVector{N_total}(
        ForwardDiff.Dual{typeof(tag)}(
            u0[i],
            ForwardDiff.Partials(
                ntuple(j -> j == partial_of[i] ? one(T) : zero(T), Val(N))
            ),
        ) for i in 1:N_total
    )

    ode_sol_dual = integrate(u0_dual, tf, vbp; save_everystep=true)

    f(t) = SVector{N}(ForwardDiff.value(ode_sol_dual(t)[i]) for i in idxs)

    J(t) =
        let ut = ode_sol_dual(t)
            SMatrix{N,N}(ForwardDiff.partials(ut[i], j) for i in idxs, j in 1:N)
        end

    return f, J, ode_sol_dual.t
end

tf, p... = variable(optim_sol)
u0 = SA[vcat(x_direct(0), 0)...]
optim_vbp = CA(vbp0; (syms .=> p)...)

# 4×4 monodromy (ignore accumulated work W)
_, J4, ts4 = flow_jacobian(u0, tf, optim_vbp; idxs=1:4)
svals4 = stack([svdvals(J4(t)) for t in ts4])'
plot(
    ts4,
    svals4;
    yscale=:log10,
    xlabel="Time (t)",
    ylabel="Singular values",
    label=["σ₁" "σ₂" "σ₃" "σ₄"],
    title="Monodromy matrix SVD over time (4×4)",
)

# 5×5 monodromy (includes W: ∂W(t)/∂W(0) = 1)
_, J5, ts5 = flow_jacobian(u0, tf, optim_vbp; idxs=1:5)
svals5 = stack([svdvals(J5(t)) for t in ts5])'
plot(
    ts5,
    svals5;
    yscale=:log10,
    xlabel="Time (t)",
    ylabel="Singular values",
    label=["σ₁" "σ₂" "σ₃" "σ₄" "σ₅"],
    title="Monodromy matrix SVD over time (5×5, σ₅≡1)",
)
