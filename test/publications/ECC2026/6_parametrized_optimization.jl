using Optimization:
    solve,
    AutoForwardDiff,
    OptimizationFunction,
    OptimizationProblem,
    remake
import OptimizationNLopt: NLopt
using ComponentArrays: ComponentArray as CA
using Logging
using StaticArrays
using Plots
using SplitApplyCombine
using LinearAlgebra: norm
using ProgressMeter
using LaTeXStrings

import KEEP.PointMass4 as PM4
import KEEP.PointMassPara as PMP
import KEEP.Optimization as OPT
import KEEP.LimitCycle as LC

includet("plots_default.jl")

function extrema_tension(sol)
    tmin, tmax = sol.prob.tspan
    tmax = tmin + (tmax - tmin) / 2  # Line tension is symmetrical, we optimize on one lobe only
    f = t -> PM4.compute_line_tension(sol(t), sol.prob.p)
    optf_min = OptimizationFunction((t, _) -> f(first(t)), AutoForwardDiff())
    optf_max = OptimizationFunction((t, _) -> -f(first(t)), AutoForwardDiff())
    prob_min = OptimizationProblem(optf_min, [(tmin + tmax) / 2], lb=[tmin], ub=[tmax])
    prob_max = OptimizationProblem(optf_max, [(tmin + tmax) / 2], lb=[tmin], ub=[tmax])
    min_ = solve(prob_min, NLopt.GN_DIRECT(), maximize=false, reltol=1e-3, maxiters=100).objective
    max_ = -solve(prob_max, NLopt.GN_DIRECT(), maximize=true, reltol=1e-3, maxiters=100).objective
    return (min_, max_)
end

function trajectory_stats(sol)
    min_tension, max_tension = extrema_tension(sol)
    Wf = sol.u[end][5]
    tf = sol.t[end]
    power = Wf / tf
    ode_plot = plot(sol, idxs=[1, 2, 3, 4], ylim=(-4, 8), xticks=0:tf, formatter=:plain, label=["α" "τ" "dα" "dτ"])
    return power, min_tension, max_tension, tf, ode_plot
end

# function compute_results(wind; syms, p0, sense, tol)
#     p0_wind = CA(p0; v_ref=wind)
#     lower, upper = OPT.make_bounds(p0_wind, syms)
#     stats, model = with_logger(NullLogger()) do
#         OPT.optimize(p0_wind, syms, lower, upper; sense=sense, tol=tol)
#     end

#     para_dims, state_dims, iterate_dims, power_dim = compute_dims(p0, syms)
#     solution_ = stats.solution
#     vbp0 = build_vbpara(p0)
#     vbp = @set vbp0[syms] = solution[end-length(syms)+1:end]

#     solution = solution_ .* iterate_dims
#     theta = solution[6:end]

#     vbp0 = PMP.build_vbpara(p0_wind)
#     vbp_optimized = @set vbp0[syms] = solution[end-length(syms)+1:end]
#     s = LC.compute_limit_cycle(PMP.build_vbpara(p); sense=sense, save_everystep=true, tol=tol)
#     s_opt = LC.compute_limit_cycle(vbp_optimized; sense=sense, save_everystep=true, tol=tol)

#     data_no_opt = trajectory_stats(s)
#     data_opt = trajectory_stats(s_opt)

#     (α0, dα0, dτ0), tf, P = s_opt.u[1][1, 3, 4], s_opt.t[end], s_opt.u[end][5] / s_opt.t[end]
#     println("Wind = {wind} m/s
#     opt error = $(norm(model.c!(model.meta.ucon, solution_)))
#     LC error =  $(norm(solution[1:5] - [α0, dα0, dτ0, tf, P]))")

#     return power, min_tension, max_tension, theta .* para_dims, ode_plot
# end

function compute_results(param, value; syms, p0, sense, tol)
    p0 = CA(p0; param => value)

    # Compute limit cycle and corresponding stats
    vbp0 = PMP.build_vbpara(p0)
    lc = LC.compute_limit_cycle(vbp0; sense=sense, save_everystep=true, tol=tol)
    stats_no_opt = trajectory_stats(lc)

    # optimize
    lower, upper = OPT.make_bounds(p0, syms)
    stats, model = with_logger(NullLogger()) do
        OPT.optimize(p0, syms, lower, upper; sense=sense, tol=tol)
    end
    para_dims, state_dims, iterate_dims, power_dim = OPT.compute_dims(p0, syms)
    solution = stats.solution
    theta = solution[6:end] .* para_dims

    # Compute limit cycle and corresponding stats for optimized parameters
    vpb = CA(vbp0; (syms .=> solution[end-length(syms)+1:end])...)
    lc_opt = LC.compute_limit_cycle(vpb; sense=sense, save_everystep=true, tol=tol)
    stats_opt = trajectory_stats(lc_opt)

    # Compute error between limit cycle found by ipopt and by own method
    (α0, dα0, dτ0), tf, P = lc_opt.u[1][[1, 3, 4]], lc_opt.t[end], lc_opt.u[end][5] / lc_opt.t[end]
    println("$param = $value
    opt error = $(norm(model.c!(model.meta.ucon, solution)))
    LC error =  $(norm((solution .* iterate_dims)[1:5] - [α0, dα0, dτ0, tf, P]))")

    return stats_no_opt, stats_opt, theta
end

function range_optim(param, values; syms, p0, sense, tol)
    stats_no_opt_vec, stats_opt_vec, theta_vec = invert(@showprogress y = [
        compute_results(param, v; syms=syms, p0=p_ref, sense=sense, tol=tol)
        for v in values
    ])
end


minmaxnormalize(x) = (x .- minimum(x)) ./ (maximum(x) - minimum(x))
maxnormalize(x) = x ./ maximum(x)

function main()
    sense = +
    syms = [:r, :I_eq, :torque_slope]  # \ell, I, b
    wind_range = range(5, 15, length=11)
    i_ref = findfirst(x -> x >= 10, wind_range)
    wind_ref = wind_range[i_ref]
    p0 = PMP.build_para(v_ref=wind_ref)  # Reference parameters at 10 m/s
    tol = 1e-7
    _, _, theta0 = compute_results(:v_ref, wind_ref; syms=syms, p0=p0, sense=sense, tol=tol)
    p_ref = CA(p0; (syms .=> theta0)...)

    stats_no_opt_vec, stats_opt_vec, theta_vec = range_optim(:v_ref, wind_range; syms=syms, p0=p_ref, sense=sense, tol=tol)

    ## Plot
    my_defaults()
    power_vec, min_tension_vec, max_tension_vec, tf_vec, ode_plot_vec = invert(stats_opt_vec)
    power_no_opt_vec, min_tension_no_opt_vec, max_tension_no_opt_vec, tf_no_opt_vec, _ = invert(stats_no_opt_vec)


    default(xformatter=x -> "", xticks=4:2:16, xlim=(4, 16))
    f1 = plot(yticks=0:20:40, ylim=(-5, 55), ylabel="Power (kW)")
    f2 = plot(yticks=0:20:40, ylim=(-5, 40), ylabel="Tension range (kN)")
    f3 = plot(yticks=((0, 1), ("0", "max")), ylim=(-0.1, 1.1), ylabel="Parameters", xlabel="Wind speed (m/s)", formatter=:plain)

    theta_vec_normalized = [maxnormalize(theta) for theta in invert(theta_vec)]

    c = :black
    plot!(f1, wind_range, power_no_opt_vec ./ 1000, c=c)
    plot!(f2, wind_range, [max_tension_no_opt_vec ./ 1000, min_tension_no_opt_vec ./ 1000], c=:black)
    # plot!(f3, wind_range, [[p for _ in wind_range] for p in invert(theta_vec_normalized)[i_ref]], c = c)

    # Plot optimization results
    plot!(f1, wind_range, power_vec ./ 1000, m=:circle, c=[1], label="Power (kW)")
    plot!(f2, wind_range, max_tension_vec ./ 1000, fillrange=min_tension_vec ./ 1000, color=1, fillcolor=1, label="Tension range (kN)", fillalpha=0.1)
    plot!(f2, wind_range, min_tension_vec ./ 1000, color=1)
    plot!(f3, wind_range, theta_vec_normalized, marker=[:rect :diamond :star6], ls=[:dash :dot :dashdot], legend=:topright, c=[1 2 3], label=[L"\ell" L"I" L"b"])

    # Legend and final plot
    legend = [:ylabel, :legend, :annotate][3]

    if legend != :ylabel
        plot!(f1, ylabel="")
        plot!(f2, ylabel="")
        plot!(f3, ylabel="")
    end
    if legend != :legend
        plot!(f1, legend=false)
        plot!(f2, legend=false)
        plot!(f3, legend=false)
    end
    if legend == :annotate
        annotate!(f1, wind_range[end] - 2, power_vec[end] ./ 1000, ("Power (kW)", TICKFONTSIZE, :right))
        annotate!(f2, wind_range[end] - 2, max_tension_vec[end] ./ 1000, ("Tension range (kN)", TICKFONTSIZE, :right))
        annotate!(f3, wind_range[end] - 2, 1.05, ("Line length", TICKFONTSIZE, :right, palette(:auto)[1]))
        annotate!(f3, wind_range[end] - 0.2, 0.7, ("Arm inertia", TICKFONTSIZE, :right, palette(:auto)[2]))
        annotate!(f3, wind_range[2], 0.3, ("Arm braking coefficient", TICKFONTSIZE, :left, palette(:auto)[3]))
    end

    fig = plot(f1, f2, f3; layout=(3, 1), size=plot_size(3 / 4))

    plot!(f1, xformatter=:plain, xlabel="Wind speed (m/s)", size=plot_size(4 / 3)) |> display
    plot!(f2, xformatter=:plain, xlabel="Wind speed (m/s)", size=plot_size(4 / 3)) |> display
    plot!(f3, size=plot_size(4 / 3)) |> display

    display(fig)

    savefig(fig, "test/publications/ECC2026/figs/parametrized_optimization.pdf")
    savefig(f1, "test/publications/ECC2026/figs/parametrized_optimization_power.pdf")
    savefig(f2, "test/publications/ECC2026/figs/parametrized_optimization_tension.pdf")
    savefig(f3, "test/publications/ECC2026/figs/parametrized_optimization_params.pdf")

    ##
    plot(xlabel="Wind speed (m/s)")
    plot!(wind_range, 1 ./ tf_no_opt_vec)
    plot!(wind_range, 1 ./ tf_vec, c=:black)
    plot!() |> display
    nothing
    #= Observation:
    For each wind speed, we maximize average power (top) using parameters (bottom), while line tension (min-max envelope, lower is better for wearing the lines less) is displayed in the middle.
    For comparison, we also compute power and tension for fixed parameters (optimized for $10$ m/s) at each wind speed (in black).
    We observed that the parameters obtained for $10$ m/s are strongly resilient to change in wind speed, as further optimizing leads to close to no improvement in power and tension.
    We observe that optimizing the parameters for each wind speed does not yield significant improvements, but we notice that the optimized parameters put less strain on the lines while producing more power.
    The system behaves very linearly along this parameterized optimization even though the length of the arm is kept constant (and other parameters), thus making it impossible to perfectly scale the model to wind speed.
    We have the following approximate correlations: power is in cube of wind speed, tensions are in square of wind speed, line length $l$ and arm inertia $I$ are inverse proportional to wind speed, arm braking coefficient $b$ is proportional to wind speed, the period (not shown) is inverse proportional to wind speed.

    \begin{itemize}
        \item power is in cube of wind speed
        \item tensions is in square of wind speed
        \item line length $l$ and arm inertia $I$ are inverse proportional to wind speed
        \item arm braking coefficient $b$ is proportional to wind speed
        \item (not shown) the period is inverse proportional to wind speed
    \end{itemize}
    =#


    ## animation of state
    # @gif for ode_plot in ode_plots
    #     plot(ode_plot)
    # end fps = 5

    my_defaults()
    plot(wind_range, [sqrt.(max_tension_vec ./ 1000) ./ 6, sqrt.(min_tension_vec ./ 1000) ./ 4], c=:blue, label="tension")
    plot!(wind_range, cbrt.(power_vec ./ 1000) ./ 4, label="power")
    plot!(wind_range, [theta[3] / 7000 for theta in theta_vec], label="arm braking coefficient")

    plot!(wind_range, [1 / (theta[1] - 45) + 0.2 for theta in theta_vec], label="line length")
    plot!(wind_range, [(125 / (theta[2] - 3270) + 0.2) / 1.2 for theta in theta_vec], label="arm inertia")
    plot!([0], [0])

    ## Optimize I for range of r

    r_range = 10:5:100
    _, _, theta_vec_ = range_optim(:r, r_range; syms=[:I_eq], p0=p_ref, sense=sense, tol=tol)

    Is = [theta[1] for theta in theta_vec_]
    inds = 30 .< r_range .< 80
    my_defaults()
    ##
    plot(r_range[inds], @. (Is[inds] / 8000 + 0.3))
    plot!([0], [0])
end

main()