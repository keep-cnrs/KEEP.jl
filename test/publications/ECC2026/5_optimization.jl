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

include("state_funcs.jl")
include("plots_default.jl")

struct OptVar
    sym::String
    unit::String
    low::Float64
    up::Float64
    guess::Float64
    sol::Float64
end

function main()
    outpath = mkpath("out/2026_07_ECC")

    syms = [:r, :I_eq, :torque_slope]
    p0 = PMP.build_para()
    lower, upper = OPT.make_bounds(p0, syms)
    tol = 1e-6
    solution, stats, model = OPT.optimize(p0, syms, lower, upper; sense=+, tol=tol)

    para_dims, state_dims, iterate_dims, power_dim = OPT.compute_dims(p0, syms)
    initial_guess_ = model.meta.x0
    initial_guess = initial_guess_ .* iterate_dims
    solution_ = stats.solution

    g_fmt = (x) -> begin
        e = floor(Int, log10(x))
        m = x / exp10(e)
        return f"${m:.1f}" * raw"\times10^{" * f"{e}" * raw"}$"
    end

    # 1. Populate the data into a Vector of OptVar
    # This acts as your single source of truth
    data = [
        OptVar(raw"$\ell$", raw"\si{\meter}", lower[1], upper[1], initial_guess[end-2], solution[end-2]),
        OptVar(raw"$I$", raw"\si{\kilogram\meter\squared}", lower[2], upper[2], initial_guess[end-1], solution[end-1]),
        OptVar(raw"$b$", raw"\si{\newton\meter\second\per\radian}", lower[3], upper[3], initial_guess[end], solution[end])
    ]

    # 2. Build Table 1 Rows (Bounds)
    bounds_rows = [f"    {v.sym} ({v.unit}) & {v.low:.1f} & {v.up:.1f} \\\\ " for v in data]

    # 3. Build Table 2 Rows (Results)
    results_rows = [f"    {v.sym} ({v.unit}) & {v.guess:.1f} & {v.sol:.1f} \\\\ " for v in data]

    # Add special rows to Results
    push!(results_rows, raw"    \midrule")
    push!(results_rows, f"    Average power (\\si{{\\kilo\\watt}}) & {initial_guess[5]/1000:.2f} & {solution[5]/1000:.2f} \\\\ ")

    res_guess = g_fmt(norm(model.c!(model.meta.ucon, initial_guess_)))
    res_sol = g_fmt(norm(model.c!(model.meta.ucon, stats.solution)))
    push!(results_rows, raw"    Norm of residual & {" * res_guess * raw"} & {" * res_sol * raw"} \\ ")

    # 4. Generate LaTeX Strings
    table_params_str = raw"""
\begin{table}[thpb]
\centering
\caption{Description of the Optimized Parameters}
\begin{tabular}{l S S}
    \toprule
    {Symbol (Unit)} & {Lower} & {Upper} \\
    \midrule
""" * join(bounds_rows, "\n") * raw""" 
    \bottomrule
\end{tabular}
\label{tab:opt_descr}
\end{table}"""

    table_opt_str = raw"""
\begin{table}[thpb]
\centering
\caption{Optimization Results}
\begin{tabular}{l S S}
    \toprule
    {Symbol (Unit)} & {Initial Guess} & {Solution} \\
    \midrule
""" * join(results_rows, "\n") * raw""" 
    \bottomrule
\end{tabular}
\label{tab:opt_results}
\end{table}"""

    println(join(["%% BEGIN TABLE 2", table_params_str, table_opt_str, "%% END TABLE 2"], "\n\n"))

    ## Plot power vs time before and after optimizing

    vbp0 = PMP.build_vbpara(p0)
    p_opt = CA(p0; solution.params...)
    vbp_opt = PMP.build_vbpara(p_opt)
    lc_0 = LC.compute_limit_cycle(vbp0; sense=+)
    tf_0 = last(lc_0.t)
    tf_opt = solution[4]
    tf = max(tf_0, tf_opt)
    sol_0 = PM4.integrate(first(lc_0.u), tf, vbp0, save_everystep=true)
    u0_opt = SA[solution[1], TAU0, solution[2], solution[3], 0]
    sol_opt = PM4.integrate(u0_opt, tf, vbp_opt, save_everystep=true)

    avg_0 = sol_0(tf_0, idxs=5) / 1000tf_0
    avg_opt = sol_opt(tf_opt, idxs=5) / 1000tf_opt

    step = 0.01
    t_0 = 0:step:tf_0
    t_opt = 0:step:tf_opt

    c_0 = :black
    c_opt = :red

    label = ff"{lbl:s}, $t_f$: {tf:.2f} s, avg: {avg:.2f} kW"
    label_0 = label((lbl="Initial guess", tf=tf_0, avg=avg_0))
    label_opt = label((lbl="Optimized", tf=tf_opt, avg=avg_opt))
    println(label_0)
    println(label_opt)

    plot(size=plot_size(4 / 3), xticks=(0:0.25:1, ["0", "0.25", "0.5", "0.75", "1"]), xlabel="Normalized time", yticks=0:5:20, ylim=(0, 20), ylabel="Power (kW)", legend=:outerbottom)
    plot!(bottom_margin=-7.5Plots.mm)

    # Shenanigans to plot Optimized on top of Initial guess, while it being first in the legend
    plot!([], [], label=label_opt, c=c_opt)
    @_ plot!(_ / tf_0, P(_, sol_0) / 1000, t_0, label=label_0, c=c_0)
    @_ plot!(_ / tf_opt, P(_, sol_opt) / 1000, t_opt, c=c_opt)

    # annotate!(0.72, avg_0 - 2, Plots.text(f"+{100*(avg_opt/avg_0 - 1)}%", ANNOTATIONFONTSIZE, c_0, :right, rotation=0, family=FONT))
    # annotate!(0.72, avg_0 - 2, Plots.text("Avg: \$$(round(avg_0, sigdigits=2))\$ kW", ANNOTATIONFONTSIZE, c_0, :right, rotation=0, family=FONT))
    # annotate!(0.28, avg_opt + 2, Plots.text("Avg: \$$(round(avg_opt, sigdigits=2))\$ kW", ANNOTATIONFONTSIZE, c_opt, :left, rotation=0, family=FONT))
    annotate!(0.9, avg_opt + 1, Plots.text(f"+{100*(avg_opt/avg_0 - 1):.0f}%", ANNOTATIONFONTSIZE, c_opt, :center, rotation=0, family=FONT))

    n = 100
    # scatter!(range(0, 1, n), fill(avg_0, n), ms=.5, c=c_0)
    # scatter!(range(0, 1, n), fill(avg_opt, n), ms=.5, c=c_opt)
    plot!([0, 1], fill(avg_0, 2), c=c_0, ls=:dot, lw=1)
    plot!([0, 1], fill(avg_opt, 2), c=c_opt, ls=:dot, lw=1)
    display(plot!())

    @show avg_0, avg_opt
    println(f"Power gain: {avg_opt:.2f}/{avg_0:.2f} = {100*avg_opt/avg_0:.2f} %")
    println("tf_opt = \$$(round(tf_opt, sigdigits=3))\$ s")

    savefig(plot!(), joinpath(outpath, "optimized_power.pdf"))
end

main()