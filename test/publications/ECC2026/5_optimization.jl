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

function main()
    syms = [:r, :I_eq, :torque_slope]
    p0 = PMP.build_para()
    lower, upper = OPT.make_bounds(p0, syms)
    tol = 1e-6
    stats, model = OPT.optimize(p0, syms, lower, upper; sense=+, tol=tol)

    para_dims, state_dims, iterate_dims, power_dim = OPT.compute_dims(p0, syms)
    initial_guess_ = model.meta.x0
    initial_guess = initial_guess_ .* iterate_dims
    solution_ = stats.solution
    solution = solution_ .* iterate_dims


    table_params_str = raw"\begin{table}[thpb]
    \centering
    \begin{tabular}{l l S S}
        \toprule
        {Name} & {Symbol [Unit]} & {Lower} & {Upper}  \\
        \midrule
        Line length & $\ell$ [\si{\meter}] & " * f"{lower[1]:.1f}" * " & " * f"{upper[1]:.1f}" * raw" \\
        Arm braking & $b$ [\si{\newton\metre\per\radian\per\second}] & " * f"{lower[2]:.1f}" * " & " * f"{upper[2]:.1f}" * raw" \\
        Arm inertia & $I$ [\si{\kilogram\metre\squared}] & " * f"{lower[3]:.1f}" * " & " * f"{upper[3]:.1f}" * raw" \\
        \bottomrule
    \end{tabular}
    \caption{Description of the optimized parameters}
    \end{table}"

    g_fmt = (x) -> begin
        e = floor(Int, log10(x))
        m = x / exp10(e)
        return f"{m:.1f}" * raw"\times10^{" * f"{e}" * "}"
    end

    table_opt_str = raw"\begin{table}[thpb]
    \centering
    \begin{tabular}{l S S}
        \toprule
        {Symbol} & {Initial guess} & {Solution} \\
        \midrule"

    table_opt_str *= "\n    " * join([
                             raw"$\ell$",
                             f"{initial_guess[end-2]:.1f}",
                             f"{solution[end-2]:.1f}"
                         ], " & ") * raw" \\ "

    table_opt_str *= "\n    " * join([
                             raw"$b$",
                             f"{initial_guess[end-1]:.1f}",
                             f"{solution[end-1]:.1f}"
                         ], " & ") * raw" \\ "

    table_opt_str *= "\n    " * join([
                             raw"$I$",
                             f"{initial_guess[end]:.1f}",
                             f"{solution[end]:.1f}"
                         ], " & ") * raw" \\ "

    table_opt_str *= raw"\midrule"

    table_opt_str *= "\n    " * join([
                             raw"Average power (W)",
                             f"{initial_guess[5]:.1f}",
                             f"{solution[5]:.1f}"
                         ], " & ") * raw" \\ "

    table_opt_str *= "\n    " * join([
                             raw"Norm of residual",
                             "{" * g_fmt(norm(model.c!(model.meta.ucon, initial_guess_))) * "}",
                             "{" * g_fmt(norm(model.c!(model.meta.ucon, solution_))) * "}"
                         ], " & ") * raw" \\ "

    table_opt_str *= raw"
        \bottomrule
    \end{tabular}
    \caption{Optimization results}
    \parbox{0.8\linewidth}{\footnotesize
      We optimize for a kite that ascends on the
      sides of the figure-eight ($\dot{tau} > 0$). Ipopt converges within $11$ iterations after $\SI{0.5}{\second} on an M2 MacBook Air for a tolerance requirement of $10^{-6}$.}
    \end{table}"

    println(join(["% ## BEGIN TABLE",
            table_params_str,
            table_opt_str,
            "% ## END TABLE"],
        "\n\n"))

    ## Plot power vs time before and after optimizing

    vbp0 = PMP.build_vbpara(p0)
    vbp_opt = CA(vbp0; (syms .=> solution_[end-length(syms)+1:end])...)
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
    t = 0:3step:tf

    c_0 = :black
    c_opt = :red

    tf_0_str = "\$t_f = $(round(tf_0, sigdigits=2))\$ s"
    tf_opt_str = "\$t_f = $(round(tf_opt, sigdigits=2))\$ s"

    plot(size=plot_size(4 / 3), xticks=(0:0.25:1, ["0", "", "", "", "1"]), xlabel="Normalized time", yticks=0:5:20, ylim=(0, 20), ylabel="Power (kW)", legend=:outerbottom, bottom_margin=-8Plots.mm)

    plot!([], [], label="Optimized, " * tf_opt_str, c=c_opt)
    @_ plot!(_ / tf_0, P(_, sol_0) / 1000, t_0, label="Initial guess, " * tf_0_str, c=c_0)
    @_ plot!(_ / tf_opt, P(_, sol_opt) / 1000, t_opt, c=c_opt)

    # annotate!(1, 12.5, (tf_0_str, ANNOTATIONFONTSIZE, c_0, :right))
    # annotate!(1, 10.2, (tf_opt_str, ANNOTATIONFONTSIZE, c_opt, :right))

    annotate!(0.72, avg_0 - 2, Plots.text("\$$(round(avg_0, sigdigits=2))\$ kW", ANNOTATIONFONTSIZE, c_0, :right, rotation=0, family=FONT))
    annotate!(0.28, avg_opt + 2, Plots.text("\$$(round(avg_opt, sigdigits=2))\$ kW", ANNOTATIONFONTSIZE, c_opt, :left, rotation=0, family=FONT))

    n = 100
    # scatter!(range(0, 1, n), fill(avg_0, n), ms=.5, c=c_0)
    # scatter!(range(0, 1, n), fill(avg_opt, n), ms=.5, c=c_opt)
    plot!([0, 1], fill(avg_0, 2), c=c_0, ls=:dash)
    plot!([0, 1], fill(avg_opt, 2), c=c_opt, ls=:dash)
    display(plot!())

    println("tf_0 = \$$(round(tf_0, sigdigits=3))\$ s")
    println("tf_opt = \$$(round(tf_opt, sigdigits=3))\$ s")

    savefig("test/publications/ECC2026/figs/optimized_power.pdf")
end

main()