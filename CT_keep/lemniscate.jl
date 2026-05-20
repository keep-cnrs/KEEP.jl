using Pkg
Pkg.activate("CT_keep")

using StaticArrays
using ForwardDiff: derivative
using OptimalControl, NLPModelsIpopt
using ComponentArrays: ComponentArray as CA
using Plots

using KEEP.PointMassPara: build_vbpara, lmt
using KEEP.PointMass4: dynamics
using KEEP.TorqueFunction: torque_function
using KEEP.LimitCycle: compute_limit_cycle
using KEEP: TAU0

const vbp0 = build_vbpara()
const syms = (:r, :I_eq, :torque_slope)

begin  # Functions
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

    f(rand(4), [15, 50, 40])
    generated_power(rand(), [15, 50, 40])

    function build_lc_ocp(; tf_init=2.8)
        ocp_lc = @def begin
            vars ∈ R², variable  # tf, dummy
            tf = vars[1]
            t ∈ [0, tf], time
            X = (α, τ, dα, dτ) ∈ R⁴, state

            # Phase fixing
            τ(0) == TAU0

            # Cyclic condition
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
            # Initialiser variable, state, control avec des vecteurs vides et ajouter au fur et a mesure ?
            # tf ∈ R, variable
            # p ∈ R³, variable
            vars ∈ R⁴, variable
            tf = vars[1]
            p = vars[2:4]
            t ∈ [0, tf], time
            # X = (α, τ, dα, dτ, W) ∈ R⁵, state
            X = (α, τ, dα, dτ) ∈ R⁴, state

            calc_vars(lc) ./ factor <= vars <= calc_vars(lc) .* factor

            τ(0) == TAU0
            α(tf) - α(0) == 0
            τ(tf) - τ(0) == 2π
            dα(tf) - dα(0) == 0
            dτ(tf) - dτ(0) == 0
            # W(0) == 0

            # Meilleur message d'erreur pour = au lieu de ==: did you mean == ?
            # X \dot et pas \dot X
            Ẋ(t) == f(X(t), p)

            # alias `->` to `→`
            ∫(generated_power(dα(t), p) / tf) → max
            # ∫(generated_power(dα(t), p) / tf) → max
        end

        init = @init ocp begin
            vars := calc_vars(lc)
            X(t) := lc(t)[1:4]
        end
        return ocp, init
    end

    function build_scaled_optim_ocp(lc; factor)
        ocp = @def begin
            vars ∈ R⁴, variable
            tf = vars[1]
            p = vars[2:4]
            t ∈ [0, 1], time                      # Fixed time interval [0, 1]
            X = (α, τ, dα, dτ) ∈ R⁴, state

            calc_vars(lc) ./ factor <= vars <= calc_vars(lc) .* factor

            τ(0) == TAU0
            α(1) - α(0) == 0
            τ(1) - τ(0) == 2π
            dα(1) - dα(0) == 0
            dτ(1) - dτ(0) == 0

            Ẋ(t) == tf * f(X(t), p)               # Explicitly scaled dynamics

            ∫(generated_power(dα(t), p)) → max     # Directly yields P_avg
        end

        # Interpolate the initial guess from [0, tf] to [0, 1]
        init = @init ocp begin
            vars := calc_vars(lc)
            X(t) := lc(t * calc_vars(lc)[1])[1:4]
        end
        return ocp, init
    end

    function solve_ocp(ocp, init; kwargs...)
        sol_init = solve(ocp; kwargs..., init=init, max_iter=0)
        sol = solve(ocp; kwargs..., init=init)
        return sol_init, sol
    end

    function make_var_plot(sol, sol_init)
        x_ticks = 1:length(syms)+1
        ratios = variable(sol) ./ variable(sol_init)
        x_labels = (x_ticks, (:tf, syms...))

        plot(
            yscale=:log10,
            xticks=x_labels,
            ylabel="Ratio (Optimized / Initial)",
            title="Parameter Optimization Results",
            framestyle=:box,
            legend=:outertopright
        )

        hline!([1 / factor], fillrange=[factor], color=:grey, alpha=0.15, label="Allowed Range")

        hline!([1], c=:black, ls=:dash, lw=1.5, label="Baseline (1.0)")

        for i in x_ticks
            plot!([i, i], [1.0, ratios[i]], c=:grey, lw=1.0, label="")
        end

        scatter!(x_ticks, ones(length(x_ticks)),
            m=:circle, mc=:white, msc=:black, ms=6, label="Initial")

        scatter!(x_ticks, ratios,
            m=:circle, mc=:crimson, ms=7, label="Optimized")

        ylims!(1 / factor / 1.2, factor * 1.2)

        return plot!()
    end

    function make_plots(sol_init, sol)
        plot(sol_init)
        state_plot = plot!(sol)
        var_plot = make_var_plot(sol, sol_init)
        return state_plot, var_plot
    end
end

oc_kwargs = (grid_size=30, backend=:generic)

lc = compute_limit_cycle(vbp0; sense=+, save_everystep=true);

factor = 5
lc_ocp, lc_init = build_lc_ocp();
optim_ocp, optim_init = build_optim_ocp(lc; factor=factor);
scaled_optim_ocp, scaled_optim_init = build_scaled_optim_ocp(lc; factor=factor);

lc_init_sol, lc_sol = solve_ocp(lc_ocp, lc_init; oc_kwargs...);
optim_init_sol, optim_sol = solve_ocp(optim_ocp, optim_init; oc_kwargs...);
scaled_optim_init_sol, scaled_optim_sol = solve_ocp(scaled_optim_ocp, scaled_optim_init; oc_kwargs...);
# optim_init_sol, optim_sol = solve_ocp(optim_ocp, optim_init; oc_kwargs..., method=:shooting, shooting_method=:single);


display.(make_plots(optim_init_sol, optim_sol))
display.(make_plots(scaled_optim_init_sol, scaled_optim_sol))

using KEEP.PointMassPara: build_para
using KEEP.Optimization: optimize
# function optimize(p, optim_para_syms, optim_para_lower, optim_para_upper; sense=+, tol=10DEFAULT_TOLERANCE, max_wall_time=120., linear_solver="mumps", initial_guess=(;))
L, M, T = lmt(vbp0)
p0 = build_para(vbp0)
units = p0[syms] ./ vbp0[syms]
lb = p0[syms] ./ factor
ub = p0[syms] .* factor
solution, stats, model = optimize(p0, syms, lb, ub)
solution.params ./ variable(optim_sol)[2:end] ./ units  # Optimum is OK

using KEEP.LimitCycle: shoot as lc_shoot

vbp = build_vbpara(CA(p0; solution.params...))
shooting = solution[1:4]
solution_sim = lc_shoot(shooting, vbp, save_everystep=true)


###########
# Same values
###########


P_keep = solution.P
P_ct = optim_sol.objective


##########
# Indirect single shooting
##########
using OrdinaryDiffEq

flow = Flow(optim_ocp);

function shoot!(res, x0, p0, tf, λ)
    # λ is previous p
    xf, pf, pλf = flow(0, x0, p0, tf, [tf, λ...]; augment=true)
    Hf = pf' * f(xf, λ) - generated_power(xf[3], λ) / tf

    res[1] = x0[2] - TAU0
    res[2:5] = xf - x0 - [0, 2π, 0, 0]
    res[6] = Hf
    res[7:9] = (pf-p0)[[1, 3, 4]]
    res[10:12] = pλf[2:4]
    return res
end

shoot!(res, y, _) = shoot!(res, y[1:4], y[5:8], y[9], y[10:12])
shoot(y) = shoot!(similar(y), y, nothing)

# Extract solution from direct method for initialization
p_direct = costate(optim_sol)
x_direct = state(optim_sol)
tf_direct, λ_direct... = variable(optim_sol)

# Initial guess
x0_guess = x_direct(0)
p0_guess = p_direct(0)
tf_guess, λ_guess = tf_direct, λ_direct

y_guess = [x0_guess; p0_guess; tf_guess; λ_guess]
shoot(y_guess)
# NLE problem with initial guess (2 unknowns: p0, λ)
prob_indirect = NonlinearProblem(shoot!,)

# Solve shooting equations
shooting_sol = solve(prob_indirect; show_trace=Val(true))
p0_sol, λ_sol = shooting_sol.u

println("Indirect solution:")
println("Initial costate: p0 = ", p0_sol)
println("Parameter: λ = ", λ_sol)