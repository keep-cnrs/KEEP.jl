## Define plot defaults
using Plots

default()

const OUTPATH = mkpath("out/2026_03_SMAI_MODE")

const A4_WIDTH_MM = 297
const PT_PER_INCH = 72.27  # DPI of Latex
const DPI = 72  # DPI of plotting library

# Calculate physical A4 width in LaTeX points (1 inch = 25.4 mm)
const A4_WIDTH_IN = A4_WIDTH_MM / 25.4
const A4_WIDTH_PT = A4_WIDTH_IN * PT_PER_INCH

# Define column widths
const WIDTH_FULL_COL_PT = A4_WIDTH_PT / 2
const WIDTH_HALF_COL_PT = A4_WIDTH_PT / 4

# Scaled Font Sizes
const FONTSIZE = 12
const LABELFONTSIZE = FONTSIZE
const TICKFONTSIZE = 8  # Scaled proportionally (was 6 when base was 8)
const LEGENDFONTSIZE = TICKFONTSIZE
const ANNOTATIONFONTSIZE = FONTSIZE

const DEFAULT_ASPECT_RATIO = 3

const FONT = "Computer Modern"

# Pass the target width to calculate dimensions (defaults to full column)
function plot_size(aspect_ratio=DEFAULT_ASPECT_RATIO; width_pt=WIDTH_FULL_COL_PT, dpi=DPI)
    width_in = width_pt / PT_PER_INCH
    width_px = round(Int, width_in * dpi)
    height_px = round(Int, width_px / aspect_ratio)
    return (width_px, height_px)
end

function my_defaults(; width_pt=WIDTH_FULL_COL_PT)
    Plots.default() # Reset defaults before applying
    return Plots.default(;
        dpi=DPI,
        fontfamily=FONT,
        size=plot_size(; width_pt=width_pt),
        formatter=:plain,
        linewidth=1,           # INCREASED: Default lines will look too thin on A4 otherwise
        markerstrokewidth=0,
        markersize=5,
        label="",
        # framestyle=:grid,      # semi, zerolines, none
        labelfontsize=LABELFONTSIZE,
        tickfontsize=TICKFONTSIZE,
        legendfontsize=LEGENDFONTSIZE,
        annotationfontsize=ANNOTATIONFONTSIZE,
        # margin=-1Plots.mm,
    )
end

# Apply defaults initially for full-column plots (A4/2)
my_defaults()

const PALETTE = Plots.palette(:auto);

## Eight plot

τ = range(0, 360; length=501)
θ = τ -> sind(2 * τ)
φ = τ -> sind(τ)

plot(; xmirror=true, size=plot_size())
# plot!(xlabel="\$ϕ(τ)\$", ylabel="\$θ(τ)\$")
plot!(;
    xticks=([-1, 0, 1], ("\$-Δφ\$", "\$φ = 0\$", "\$Δφ\$")),
    yticks=([-1, 0, 1], ("\$θ_0 + Δθ\$", "\$θ = θ_0\$", "\$θ_0 - Δθ\$")),
)

plot!(φ.(τ), θ.(τ); c=:grey)

delta_col = :black
τ_φ = 90.0
plot!(φ.([0, τ_φ]), θ.([τ_φ, τ_φ]); arrow=false, color=delta_col, label="")
annotate!(
    sum(φ.([0, τ_φ])) / 2,
    sum(θ.([τ_φ, τ_φ])) / 2 - 0.1,
    ("\$Δφ\$", ANNOTATIONFONTSIZE, delta_col, :top),
)

τ_θ = 225.0
plot!(φ.([τ_θ, τ_θ]), θ.([τ_θ, 0]); arrow=false, color=delta_col, label="")
annotate!(
    sum(φ.([τ_θ, τ_θ])) / 2 + 0.05,
    sum(θ.([τ_θ, 0])) / 2,
    ("\$Δθ\$", ANNOTATIONFONTSIZE, delta_col, :left),
)

dx, dy = 1 / 4, 1 / 2
up_col, down_col = 1, 2
quiver!([0], [0]; quiver=([dx], [-dy]), c=up_col)
quiver!([0], [0]; quiver=([-dx], [dy]), c=down_col)
annotate!(
    dx - 0.05,
    -dy - 0.25,
    ("\$\\dot{\\tau} > 0\$", ANNOTATIONFONTSIZE, :center, palette(:auto)[up_col]),
)
annotate!(
    -dx + 0.05,
    dy + 0.25,
    ("\$\\dot{\\tau} < 0\$", ANNOTATIONFONTSIZE, :center, palette(:auto)[down_col]),
)

display(plot!())

savefig(joinpath(OUTPATH, "eight.pdf"))

## Bare figure eight

plot(; xmirror=true, size=plot_size(), xticks=nothing, yticks=nothing, framestyle=:none)

plot!(φ.(τ), θ.(τ); c=:grey)

dx, dy = 1 / 4, 1 / 2
up_col, down_col = 1, 2
quiver!([0], [0]; quiver=([dx], [-dy]), c=up_col)
quiver!([0], [0]; quiver=([-dx], [dy]), c=down_col)
annotate!(
    dx - 0.05,
    -dy - 0.25,
    ("\$\\dot{\\tau} > 0\$", ANNOTATIONFONTSIZE, :center, palette(:auto)[up_col]),
)
annotate!(
    -dx + 0.05,
    dy + 0.25,
    ("\$\\dot{\\tau} < 0\$", ANNOTATIONFONTSIZE, :center, palette(:auto)[down_col]),
)

savefig(joinpath(OUTPATH, "bare-eight.pdf"))

## 3D Plots

using StaticArrays: SA
using SplitApplyCombine: invert

import KEEP.PointMassPara as PMP
import KEEP.PointMass4 as PM4
import KEEP.Visualization as VIS

p = PMP.build_vbpara()
α0, τ0, dα0, dτ0 = 0, 1, 0, 0
tf = 100
sol = PM4.integrate(SA[α0, 1, dα0, dτ0, 0], tf, p; save_everystep=true)

VIS.plot_trajectory_4D(sol; tspan=10)
plot!(; size=plot_size(4 / 3), margin=-5Plots.mm)
savefig(joinpath(OUTPATH, "trajectory.pdf"))

## Power plot

VIS.plot_avg_power_4D(sol; tspan=10)
plot!(; size=plot_size(4 / 3), bottommargin=-5Plots.mm)
savefig(joinpath(OUTPATH, "avg_power.pdf"))

## yz_limit_cycles
import KEEP.LimitCycle as LC

include("ECC2026/state_funcs.jl")

lc_p = LC.compute_limit_cycle(p; sense=(+), save_everystep=true)
lc_m = LC.compute_limit_cycle(p; sense=(-), save_everystep=true)

function yz_split(sol, n_splits)
    t = range(sol.prob.tspan..., 300)
    yz = [r_vec(t, sol)[2:3] for t in t]
    Δt = (t[end] - t[1]) / (n_splits + 1)
    for i in 1:n_splits
        t_split = (i - 0.2) * Δt
        ind = searchsortedfirst(t, t_split)
        insert!(yz, ind + i - 1, [NaN, NaN])
    end
    return yz
end

n_splits = 4
arrows = [permutedims(repeat([true], n_splits - 1)) false]
plot(;
    aspect_ratio=:equal,
    xlabel="\$y\$ (m)",
    ylabel="\$z\$ (m)",
    xtick=-20:5:20,
    ytick=5:2.5:10,
    size=plot_size(3),
)
plot!(invert(yz_split(lc_p, n_splits))...; arrows=true, label="\$\\dot{\\tau} > 0\$")
plot!(invert(yz_split(lc_m, n_splits))...; arrows=true, label="\$\\dot{\\tau} < 0\$")
plot!([0], [13])
lens!([-11, -10], [8, 9]; inset=(1, bbox(0.3, 0, 0.3, 0.4)))
l = plot!().inset_subplots[1]
xticks!(l, -11:0.5:-10)
yticks!(l, 8:0.5:9)
plot!(l; tickfontsize=TICKFONTSIZE - 1)
plot!(; bottommargin=2Plots.mm)
display(plot!())

savefig(joinpath(OUTPATH, "yz_limit_cycles.pdf"))

## Optimization
using PyFormattedStrings
using LinearAlgebra: norm
using Plots
using StaticArrays
using Underscores
using ComponentArrays: ComponentArray as CA
using LaTeXStrings

using KEEP: TAU0
import KEEP.LimitCycle as LC
import KEEP.Optimization as OPT
import KEEP.PointMass4 as PM4
import KEEP.PointMassPara as PMP
import KEEP.TorqueFunction as TF

syms = [:r, :I_eq, :torque_slope]
p0 = PMP.build_para()
lower, upper = OPT.make_bounds(p0, syms)
tol = 1e-6
solution, stats, model = OPT.optimize(p0, syms, lower, upper; sense=(+), tol=tol)

para_dims, state_dims, iterate_dims, power_dim = OPT.compute_dims(p0, syms)
initial_guess_ = model.meta.x0
initial_guess = initial_guess_ .* iterate_dims
solution_ = stats.solution

## Plot power vs time before and after optimizing

vbp0 = PMP.build_vbpara(p0)
p_opt = CA(p0; solution.params...)
vbp_opt = PMP.build_vbpara(p_opt)
lc_0 = LC.compute_limit_cycle(vbp0; sense=(+))
tf_0 = last(lc_0.t)
tf_opt = solution[4]
tf = max(tf_0, tf_opt)
sol_0 = PM4.integrate(first(lc_0.u), tf, vbp0; save_everystep=true)
u0_opt = SA[solution[1], TAU0, solution[2], solution[3], 0]
sol_opt = PM4.integrate(u0_opt, tf, vbp_opt; save_everystep=true)

avg_0 = sol_0(tf_0; idxs=5) / 1000tf_0
avg_opt = sol_opt(tf_opt; idxs=5) / 1000tf_opt

step = 0.01
t_0 = 0:step:tf_0
t_opt = 0:step:tf_opt

c_0 = :black
c_opt = :red

tf_0_str = "\$t_f: $(round(tf_0, sigdigits=2))\$ s"
tf_opt_str = "\$t_f: $(round(tf_opt, sigdigits=2))\$ s"
avg_0_str = "Avg: $(round(avg_0, sigdigits=2)) kW"
avg_opt_str = "Avg: $(round(avg_opt, sigdigits=2)) kW"
label_0 = "Initial guess, " * tf_0_str * ", " * avg_0_str
label_opt = "Optimized, " * tf_opt_str * ", " * avg_opt_str
# label = ff"{:s}, $t_f$: {:.2f}, avg: {:.2f}"
label = ff"{lbl:s}, $t_f$: {tf:.2f}, avg: {avg:.2f}"
label_0 = label((lbl="Initial guess", tf=tf_0, avg=avg_0))
label_opt = label((lbl="Optimized", tf=tf_opt, avg=avg_opt))

plot(;
    size=plot_size(4 / 3),
    xticks=(0:0.25:1, ["0", "0.25", "0.5", "0.75", "1"]),
    xlabel="Normalized time",
    yticks=0:5:20,
    ylim=(0, 20),
    ylabel="Power (kW)",
    legend=:outerbottom,
)
plot!(; bottom_margin=-7.5Plots.mm)

# Shenanigans to plot Optimized on top of Initial guess, while it being first in the legend
plot!([], []; label=label_opt, c=c_opt)
@_ plot!(_ / tf_0, P(_, sol_0) / 1000, t_0, label=label_0, c=c_0)
@_ plot!(_ / tf_opt, P(_, sol_opt) / 1000, t_opt, c=c_opt)

# annotate!(0.72, avg_0 - 2, Plots.text(f"+{100*(avg_opt/avg_0 - 1)}%", ANNOTATIONFONTSIZE, c_0, :right, rotation=0, family=FONT))
# annotate!(0.72, avg_0 - 2, Plots.text("Avg: \$$(round(avg_0, sigdigits=2))\$ kW", ANNOTATIONFONTSIZE, c_0, :right, rotation=0, family=FONT))
# annotate!(0.28, avg_opt + 2, Plots.text("Avg: \$$(round(avg_opt, sigdigits=2))\$ kW", ANNOTATIONFONTSIZE, c_opt, :left, rotation=0, family=FONT))
annotate!(
    0.375,
    avg_opt + 2,
    Plots.text(
        f"+{100*(avg_opt/avg_0 - 1):.0f}%",
        ANNOTATIONFONTSIZE,
        c_opt,
        :center;
        rotation=0,
        family=FONT,
    ),
)

n = 100
# scatter!(range(0, 1, n), fill(avg_0, n), ms=.5, c=c_0)
# scatter!(range(0, 1, n), fill(avg_opt, n), ms=.5, c=c_opt)
plot!([0, 1], fill(avg_0, 2); c=c_0, ls=:dot, lw=1)
plot!([0, 1], fill(avg_opt, 2); c=c_opt, ls=:dot, lw=1)
plot!(; size=plot_size(2))

display(plot!())
savefig(joinpath(OUTPATH, "optimized_power.pdf"))

## Parametric Optimization

using Optimization:
    solve, AutoForwardDiff, OptimizationFunction, OptimizationProblem, remake
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

function extrema_tension(sol)
    tmin, tmax = sol.prob.tspan
    tmax = tmin + (tmax - tmin) / 2  # Line tension is symmetrical, we optimize on one lobe only
    f = t -> PM4.compute_line_tension(sol(t), sol.prob.p)
    optf_min = OptimizationFunction((t, _) -> f(first(t)), AutoForwardDiff())
    optf_max = OptimizationFunction((t, _) -> -f(first(t)), AutoForwardDiff())
    prob_min = OptimizationProblem(optf_min, [(tmin + tmax) / 2]; lb=[tmin], ub=[tmax])
    prob_max = OptimizationProblem(optf_max, [(tmin + tmax) / 2]; lb=[tmin], ub=[tmax])
    min_ =
        solve(prob_min, NLopt.GN_DIRECT(); maximize=false, reltol=1e-3, maxiters=100).objective
    max_ =
        -solve(prob_max, NLopt.GN_DIRECT(); maximize=true, reltol=1e-3, maxiters=100).objective
    return (min_, max_)
end

function trajectory_stats(sol)
    min_tension, max_tension = extrema_tension(sol)
    Wf = sol.u[end][5]
    tf = sol.t[end]
    power = Wf / tf
    ode_plot = plot(
        sol;
        idxs=[1, 2, 3, 4],
        ylim=(-4, 8),
        xticks=0:tf,
        formatter=:plain,
        label=["α" "τ" "dα" "dτ"],
    )
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
    _, stats, model = with_logger(NullLogger()) do
        return OPT.optimize(p0, syms, lower, upper; sense=sense, tol=tol)
    end
    para_dims, state_dims, iterate_dims, power_dim = OPT.compute_dims(p0, syms)
    solution = stats.solution
    theta = solution[6:end] .* para_dims

    # Compute limit cycle and corresponding stats for optimized parameters
    vpb = CA(vbp0; (syms .=> solution[(end - length(syms) + 1):end])...)
    lc_opt = LC.compute_limit_cycle(vpb; sense=sense, save_everystep=true, tol=tol)
    stats_opt = trajectory_stats(lc_opt)

    # Compute error between limit cycle found by ipopt and by own method
    (α0, dα0, dτ0), tf, P = lc_opt.u[1][[1, 3, 4]],
    lc_opt.t[end],
    lc_opt.u[end][5] / lc_opt.t[end]
    println("$param = $value
    opt error = $(norm(model.c!(model.meta.ucon, solution)))
    LC error =  $(norm((solution .* iterate_dims)[1:5] - [α0, dα0, dτ0, tf, P]))")

    return stats_no_opt, stats_opt, theta
end

function range_optim(param, values; syms, p_ref, sense, tol)
    return stats_no_opt_vec, stats_opt_vec, theta_vec = invert(
        @showprogress [
            compute_results(param, v; syms=syms, p0=p_ref, sense=sense, tol=tol) for
            v in values
        ]
    )
end

minmaxnormalize(x) = (x .- minimum(x)) ./ (maximum(x) - minimum(x))
maxnormalize(x) = x ./ maximum(x)

sense = +
syms = [:r, :I_eq, :torque_slope]  # \ell, I, b
wind_range = range(5, 15; length=11)
i_ref = findfirst(x -> x >= 10, wind_range)
wind_ref = wind_range[i_ref]
p0 = PMP.build_para(; v_ref=wind_ref)  # Reference parameters at 10 m/s
tol = 1e-7
_, _, theta0 = compute_results(:v_ref, wind_ref; syms=syms, p0=p0, sense=sense, tol=tol)
p_ref = CA(p0; (syms .=> theta0)...)

stats_no_opt_vec, stats_opt_vec, theta_vec = range_optim(
    :v_ref, wind_range; syms=syms, p_ref=p_ref, sense=sense, tol=tol
)

## Plot
power_vec, min_tension_vec, max_tension_vec, tf_vec, ode_plot_vec = invert(stats_opt_vec)
power_no_opt_vec, min_tension_no_opt_vec, max_tension_no_opt_vec, tf_no_opt_vec, _ = invert(
    stats_no_opt_vec
)

default(; xformatter=x -> "", xticks=4:2:16, xlim=(4, 16))
f1 = plot(; yticks=0:20:40, ylim=(-5, 55), ylabel="Power (kW)")
f2 = plot(; yticks=0:20:40, ylim=(-5, 40), ylabel="Tension range (kN)")
f3 = plot(;
    yticks=((0, 1), ("0", "max")),
    ylim=(-0.1, 1.1),
    ylabel="Parameters",
    xlabel="Wind speed (m/s)",
    formatter=:plain,
)

theta_vec_normalized = [maxnormalize(theta) for theta in invert(theta_vec)]

ref_style = (; ls=:dot, c=:black)
plot!(f1, wind_range, power_no_opt_vec ./ 1000; ref_style...)
plot!(
    f2,
    wind_range,
    [max_tension_no_opt_vec ./ 1000, min_tension_no_opt_vec ./ 1000];
    ref_style...,
)
plot!(
    f3,
    wind_range,
    [[p for _ in wind_range] for p in invert(theta_vec_normalized)[i_ref]];
    lw=0.5,
    ref_style...,
)

# Plot optimization results
plot!(f1, wind_range, power_vec ./ 1000; m=:circle, c=[1], label="Power (kW)")
plot!(
    f2,
    wind_range,
    max_tension_vec ./ 1000;
    fillrange=min_tension_vec ./ 1000,
    color=1,
    fillcolor=1,
    label="Tension range (kN)",
    fillalpha=0.1,
)
plot!(f2, wind_range, min_tension_vec ./ 1000; color=1)
plot!(
    f3,
    wind_range,
    theta_vec_normalized;
    marker=[:rect :diamond :star6],
    ls=[:dash :dot :dashdot],
    legend=:topright,
    c=[1 2 3],
    label=[L"\ell" L"I" L"b"],
)

# Legend and final plot
legend = [:ylabel, :legend, :annotate][3]

if legend != :ylabel
    plot!(f1; ylabel="")
    plot!(f2; ylabel="")
    plot!(f3; ylabel="")
end
if legend != :legend
    plot!(f1; legend=false)
    plot!(f2; legend=false)
    plot!(f3; legend=false)
end
if legend == :annotate
    annotate!(
        f1,
        wind_range[end] - 2,
        power_vec[end] ./ 1000,
        ("Power (kW)", TICKFONTSIZE, :right),
    )
    annotate!(
        f2,
        wind_range[end] - 2,
        max_tension_vec[end] ./ 1000,
        ("Tension range (kN)", TICKFONTSIZE, :right),
    )
    annotate!(
        f3,
        wind_range[end] - 2,
        1.0,
        ("Line length", TICKFONTSIZE, :right, palette(:auto)[1]),
    )
    annotate!(
        f3,
        wind_range[end] - 0.2,
        0.74,
        ("Arm inertia", TICKFONTSIZE, :right, palette(:auto)[2]),
    )
    annotate!(
        f3,
        wind_range[2],
        0.35,
        ("Arm braking coefficient", TICKFONTSIZE, :left, palette(:auto)[3]),
    )
end

display(plot!(f1; xformatter=:plain, xlabel="Wind speed (m/s)", size=plot_size(4 / 3)))
display(plot!(f2; xformatter=:plain, xlabel="Wind speed (m/s)", size=plot_size(4 / 3)))
display(plot!(f3; size=plot_size(4 / 3)))

fig = plot(f1, f2, f3; layout=(1, 3), size=plot_size(3; width_pt=2WIDTH_FULL_COL_PT))
plot!(fig; bottommargin=5Plots.mm)
display(fig)

savefig(fig, joinpath(OUTPATH, "parametrized_optimization.pdf"))
# savefig(f1, joinpath(OUTPATH, "parametrized_optimization_power.pdf"))
# savefig(f2, joinpath(OUTPATH, "parametrized_optimization_tension.pdf"))
# savefig(f3, joinpath(OUTPATH, "parametrized_optimization_params.pdf"))
