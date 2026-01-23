using StaticArrays
using Plots
using LaTeXStrings
using PyFormattedStrings
using Underscores
using Setfield
using SplitApplyCombine

# DO for u0 and -u0
using KEEP: TAU0
import KEEP.LimitCycle as LC
import KEEP.PointMass4 as PM4

include("state_funcs.jl")
include("plots_default.jl")

function main()
    u0 = SA[0, TAU0, 0, 1, 0]
    p = PM4.build_vbpara()
    tf = 25
    cb = LC.build_poincare_callback(eps(1.0))
    sol = PM4.integrate(-u0, tf, p; callback=cb, save_end=false, tol=1e-13)

    x_infty = sol.u[end]
    errors = [LC.distance_on_section(u, x_infty) for u in sol.u[1:end-1]]
    n = length(errors)

    i1, i2 = 2, 6
    convergence_rate = -(log10(errors[i2]) - log10(errors[i1])) / (i2 - i1)

    plot(yscale=:log10, xlim=(0, n + 1), xticks=0:2:100, ylim=(1e-16, 1e2), yticks=10.0 .^ (0:-3:-16), xlabel="lemniscate revolutions", ylabel="distance to limit")
    plot!(i -> errors[i1] * exp10(-convergence_rate * (i - i1)), label="", c=:black)
    annotate!((3.9, 1e-8, (f"slope = {convergence_rate:.2f}", :right, ANNOTATIONFONTSIZE, colorant"black")))
    scatter!(errors .+ 1e-16, label=L"|| x_i - x_\infty  ||")
    display(plot!())
    savefig("test/publications/ECC2026/figs/limit_cycle_convergence.pdf")

    tabular_str = raw"\begin{table}[thpb]
    \centering
    \begin{tabular}{|c|c c|}
        \hline
        $i$ & $t_i$ & $\norm{x_i - x_\infty}$ \\ [0.5ex]
        \midline"
    for i in 1:n
        t_i = sol.t[i]
        x_i = sol.u[i]
        x_infty = sol.u[end]
        norm_x_i_x_infty = LC.distance_on_section(x_i, x_infty)
        tabular_str *= f"
        {i} & {t_i:.2f} & {norm_x_i_x_infty:.2e} \\\\"
    end
    tabular_str *= raw"
    \end{tabular}
    \end{table}"

    println("% ## BEGIN TABLE\n\n", tabular_str, "\n\n% ## END TABLE")

    ## Count the limit cycles
    println("Number of limit cycles: ", length(LC.all_limit_cycles(p)))

    ## Comparing trajectories
    # (+) trajectory
    # (-) trajectory
    # (-) trajectory with necessary symmetry (α, τ) -> (-α, -τ) to show that they are not the "inverse" of one another
    lc_p = LC.compute_limit_cycle(p; sense=+, save_everystep=true)
    lc_m = LC.compute_limit_cycle(p; sense=-, save_everystep=true)


    P_fig = plot(ylabel="P (kW)", xformatter=x -> "", yticks=0:10:20, ylim=(-3, 23))
    @_ plot!(τ(_, lc_m) + 2π, P(_, lc_m) / 1000, lc_m.t, label="", c=2)
    @_ plot!(τ(_, lc_p), P(_, lc_p) / 1000, lc_p.t, label="", c=1)
    @_ plot!(-τ(_, lc_m), P(_, lc_m) / 1000, lc_m.t, label="", c=2, ls=:dot)

    α_fig = plot(ylabel="α (rad)", xlabel="τ (rad)", yticks=-1:1:1, ylim=(-1.1, 1.1))
    @_ plot!(τ(_, lc_m) + 2π, α(_, lc_m), lc_m.t, label="", c=2)
    @_ plot!(τ(_, lc_p), α(_, lc_p), lc_p.t, label="", c=1)
    @_ plot!(-τ(_, lc_m), -α(_, lc_m), lc_m.t, label="", c=2, ls=:dot)

    plot(P_fig, α_fig, layout=(2, 1), size=(400, 300))
    display(plot!())


    ## Eight-trajectory in YZ plane
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
    plot(aspect_ratio=:equal, xlabel="\$y\$ (m)", ylabel="\$z\$ (m)", xtick=-20:5:20, ytick=5:2.5:10, size=plot_size(1.8))
    plot!(invert(yz_split(lc_p, n_splits))..., arrows=true, label="\$\\dot{\\tau} > 0\$")
    plot!(invert(yz_split(lc_m, n_splits))..., arrows=true, label="\$\\dot{\\tau} < 0\$")
    plot!([0], [13])
    lens!([-11, -10], [8, 9], inset=(1, bbox(0.3, 0, 0.3, 0.4)))
    l = plot!().inset_subplots[1]
    xticks!(l, -11:0.5:-10)
    yticks!(l, 8:0.5:9)
    plot!(l, tickfontsize=TICKFONTSIZE - 1)

    display(plot!())

    savefig(plot!(), "test/publications/ECC2026/figs/yz_limit_cycles.pdf")

    #=
    Projection of the kite trajectories onto the YZ plane for the two limit cycles, with equal scaling on both axes.
    Arrows indicate the direction of the kite's motion; zooming in on the trajectory shows that the two paths are distinct.
    =#


    ## Map limit cycles direction in τ

    # more accurate version: do not compute the limit cycle, just integrate for enough time
    # lc_sign(q, p) = first(LC.compute_limit_cycle(SA[q..., 0, 0, 0], p; tol=1e-3))[4] > 0 ? 1 : -1

    lc_sign(q, p) = last(PM4.integrate(SA[q..., 0, 0, 0], 3, p; tol=1e-2))[4] > 0 ? 1 : -1

    n = 10
    qs = [SA[α, τ] for α in range(-π, π, length=100) for τ in range(0, 2π, length=100)]

    @time lc_signs = lc_sign.(qs, Ref(p))

    scatter(invert(qs)..., c=lc_signs, label="", aspect_ratio=:equal, ms=1)
    display(plot!())
end

main()