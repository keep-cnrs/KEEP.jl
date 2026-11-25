# julia --project=. paper.jl #src
# pandoc paper.md -o paper.html --standalone --mathjax --css="https://cdn.jsdelivr.net/npm/water.css@2/out/water.css" #src

if abspath(PROGRAM_FILE) == (@__FILE__) && get(ENV, "PAPER_LITERATE_EXEC", "0") != "1" #src
    ENV["PAPER_LITERATE_EXEC"] = "1" #src
    using Literate #src
    # Literate hides code line by line; collapse multi-line text output into one line so #hide applies to all of it. #src
    function inline_text_sources(s::AbstractString) #src
        lines = split(s, '\n'; keepempty=true) #src
        out = String[]; i = 1 #src
        while i <= length(lines) #src
            if endswith(lines[i], "Markdown.parse(\"\"\"") #src
                j = i + 1; body = String[] #src
                while j <= length(lines) && strip(lines[j]) != "\"\"\") #hide" #src
                    push!(body, lines[j]); j += 1 #src
                end #src
                j <= length(lines) || error("unterminated Markdown.parse block at line $i") #src
                push!(out, "Markdown.parse(\"\"\"" * join(body, "\\n") * "\\n\"\"\") #hide") #src
                i = j + 1 #src
            else #src
                push!(out, lines[i]); i += 1 #src
            end #src
        end #src
        return join(out, '\n') #src
    end #src
    Literate.markdown(@__FILE__, @__DIR__; execute=true, flavor=Literate.CommonMarkFlavor(), credit=false, preprocess=inline_text_sources) #src
    exit(0) #src
end #src

# # Rocking-Arm Kite Power Generation (KEEP): Compact Numerical Analysis
#
# **Research question:** *Given a kite tethered to a rotating ground-level arm that is pulled left and right across the wind window, how much average power can we expect, and how robust is the design?*
#
# Structure follows the finished ECC paper (model, equilibria and limit cycles, limit-cycle optimization) in a condensed form; the full derivations live in the JTCAM draft. Every figure below is accompanied by the numerical data needed to write the final prose without re-running any computation. Plot styling is local to each section so it can be iterated on by hand.

using Pkg #hide
Pkg.activate(@__DIR__) #hide
using KEEP: PointMass10 as PM10 #hide
using KEEP: PointMass4 as PM4 #hide
using KEEP: PointMassPara as PMP #hide
using KEEP: SteadyState #hide
using KEEP: LimitCycle as LC #hide
using KEEP: Optimization as OPT #hide
using KEEP: Visualization as VIZ #hide
using KEEP: DEFAULT_TOLERANCE, TAU0 #hide
using OrdinaryDiffEqVerner: Vern9, ODEProblem #hide
using OrdinaryDiffEqTsit5: Tsit5 #hide
using QuasiMonteCarlo: sample, HaltonSample #hide
using SplitApplyCombine: invert #hide
using BenchmarkTools: @belapsed #hide
using LinearAlgebra: norm, eigvals, × #hide
using Logging: with_logger, NullLogger #hide
ENV["GKSwstype"] = "100" #hide
using Plots #hide
using PyFormattedStrings: @f_str #hide
using Markdown #hide
import ForwardDiff #hide
include(joinpath(@__DIR__, "plot_helpers.jl")) #hide
include(joinpath(@__DIR__, "paper_core.jl")) #hide
include(joinpath(@__DIR__, "Section_10D.jl")) #hide
include(joinpath(@__DIR__, "Section_Modelling.jl")) #hide
include(joinpath(@__DIR__, "Section_Numerical_analysis.jl")) #hide
include(joinpath(@__DIR__, "Section_Optimization.jl")) #hide
include(joinpath(@__DIR__, "paper_cache.jl")) #hide
include(joinpath(@__DIR__, "paper_output.jl")) #hide
my_defaults() #hide

# Shared guide lines for the optimization panels. Kept here, next to the plots
# that use it, because it is pure styling.
const GUIDES = (; ymin=0.9, line=(; c=:black, lw=5, alpha=0.3), floor=(; c=:black, lw=1, alpha=1)) #hide

# One call site for every numerical study. FORCE_RECOMPUTE=1 recomputes all;
# otherwise a single hash of all compute inputs decides. See paper_cache.jl.
function compute_all() #hide
    (D10=compute_10D(), MOD=compute_modelling(), NUM=compute_numerical(), OPT=compute_optimization()) #hide
end #hide
results = load_or_compute() #hide
nothing #hide

# ## 1. Model and Geometry
#
# The kite is constrained to a lemniscate on a cone; the arm rotates in the horizontal plane. The two figures below fix the geometry and the meaning of the shape parameters $(\theta_0, \Delta\theta, \Delta\varphi)$ and of the phase $\tau$.

# Geometry of the reachable set. `a = 2 m`, `ℓ = 3 m` are illustrative, not the optimized values.
function plot_spindle_torus(r) #hide
    m = r.MOD #hide
    fig = plot(m.xs, m.ys, m.zs; c=:black, alpha=0.4, lw=0.5, aspect_ratio=:equal, #hide
        xlabel="x (m)", ylabel="y (m)", zlabel="z (m)", #hide
        title="Spindle torus, not to scale", label="", legend=:outerright) #hide
    plot!(fig, m.a * cos.(m.θc), m.a * sin.(m.θc), zero(m.θc); c=:blue, lw=3, label="") #hide
    plot!(fig, [0, m.a], [0, 0], [0, 0]; c=:black, lw=3, label="") #hide
    plot!(fig, [m.a, m.kite_pos[1]], [0, 0], [0, m.kite_pos[3]]; c=:black, lw=3, label="") #hide
    plot!(fig, m.a .+ m.l .* cos.(m.θc), zero(m.θc), m.l .* sin.(m.θc); c=:orange, lw=2, label="") #hide
    plot!(fig, (m.a + m.l * cos(m.β0)) .* cos.(m.θc), (m.a + m.l * cos(m.β0)) .* sin.(m.θc), #hide
        fill(m.l * sin(m.β0), length(m.θc)); c=:green, lw=2, label="") #hide
    scatter!(fig, [m.kite_pos[1]], [0], [m.kite_pos[3]]; c=:red, ms=7, label="") #hide
    scatter!([]; c=:red, ms=7, legend_marker(7)...) #hide
    plot!([]; c=:black, lw=3, label="Arm and lines") #hide
    plot!([]; c=:blue, lw=3, label="Arm-tip locus") #hide
    plot!([]; c=:orange, lw=2, label="Kite longitude") #hide
    plot!([]; c=:green, lw=2, label="Kite latitude") #hide
    plot!(leftmargin=-4Plots.mm, rightmargin=-20Plots.mm) #hide
    plot!(size=plot_size(1, 1.06)) #hide
    plot!(fig, xticks=edge_ticks(fig[1], :x; include_zero=true), #hide
        yticks=edge_ticks(fig[1], :y; include_zero=true), #hide
        zticks=edge_ticks(fig[1], :z; include_zero=true)) #hide
    return fig #hide
end #hide

function plot_eight_plot(r) #hide
    m = r.MOD #hide
    plot(xmirror=true) #hide
    plot!(xticks=([-1, 0, 1], (raw"$-Δφ$", raw"$φ = 0$", raw"$Δφ$")), #hide
        yticks=([-1, 0, 1], (raw"$θ_0 + Δθ$", raw"$θ = θ_0$", raw"$θ_0 - Δθ$"))) #hide
    plot!(m.phi_vals, m.theta_vals; c=:grey) #hide
    plot!(sind.([0, m.τ_φ]), sind.(2 .* [m.τ_φ, m.τ_φ]); color=:black) #hide
    annotate!(sum(sind.([0, m.τ_φ])) / 2, sum(sind.(2 .* [m.τ_φ, m.τ_φ])) / 2 - 0.1, #hide
        (raw"$Δφ$", ANNOTATIONFONTSIZE, :black, :top)) #hide
    plot!(sind.([m.τ_θ, m.τ_θ]), sind.(2 .* [m.τ_θ, 0]); color=:black) #hide
    annotate!(sum(sind.([m.τ_θ, m.τ_θ])) / 2 + 0.05, sum(sind.(2 .* [m.τ_θ, 0])) / 2, #hide
        (raw"$Δθ$", ANNOTATIONFONTSIZE, :black, :left)) #hide
    quiver!([0], [0]; quiver=([m.dx], [-m.dy]), c=m.up_col) #hide
    quiver!([0], [0]; quiver=([-m.dx], [m.dy]), c=m.down_col) #hide
    annotate!(m.dx - 0.05, -m.dy - 0.25, (raw"$\dot{\tau} > 0$", ANNOTATIONFONTSIZE, :center, PALETTE[m.up_col])) #hide
    annotate!(-m.dx + 0.05, m.dy + 0.25, (raw"$\dot{\tau} < 0$", ANNOTATIONFONTSIZE, :center, PALETTE[m.down_col])) #hide
    plot!(size=SIZE_HALF .* (1, 0.5), topmargin=1Plots.mm, leftmargin=5Plots.mm) #hide
    return current() #hide
end #hide

save_fig(plot_spindle_torus(results), "spindle_torus") #hide
save_fig(plot_eight_plot(results), "eight_plot") #hide

Markdown.parse("""
**Figure data — geometry.** Arm \$a = $(results.MOD.a)\$ m, tether \$\\ell = $(results.MOD.l)\$ m; the sampled elevation is \$\\beta_0 = $(round(rad2deg(results.MOD.β0); digits=0))^\\circ\$, so the kite sits at $([round(x; digits=2) for x in results.MOD.kite_pos]) m. The reachable set is a spindle torus because \$\\ell > a\$; the arm-tip locus is the blue circle. The lemniscate uses \$\\tau_\\varphi = $(results.MOD.τ_φ)^\\circ\$ and \$\\tau_\\theta = $(results.MOD.τ_θ)^\\circ\$ as illustration angles; increasing \$\\dot\\tau\$ climbs in the middle, decreasing \$\\dot\\tau\$ climbs on the sides.

**Prose hint:** these two panels replace the configuration sketch and the lemniscate sketch of the ECC paper; cite them when introducing \$(\\alpha, \\tau)\$ and \$(\\theta_0, \\Delta\\theta, \\Delta\\varphi)\$.

![Spindle torus](figs/spindle_torus.svg)

![Lemniscate parametrisation](figs/eight_plot.svg)
""") #hide

# ## 2. Equilibria and Limit Cycles
#
# All computations use the kinetic baseline $b = I\\omega$ (`PMP.optimal_torque_slope`), which gives one limit cycle per iteration. The phase-space picture, the two attracting cycles and their convergence rate are the ingredients needed for the dynamics section of the final paper.

my_theme = (; #hide
    downwind=(; c=plot_color(:grey, 0.08)), #hide
    upwind=(; c=plot_color(:red, 0.08)), #hide
    streamlines=(; c=:black, alpha=0.7, lw=Plots.default(:lw) / 6), #hide
    heads=(; c=:black, alpha=0.5, lw=0), #hide
    nullcline_tau=(; c=1, ls=:solid), #hide
    nullcline_alpha=(; c=2, ls=:solid), #hide
    eq_legend=(; c=:black, shape=:circle, ms=Plots.default(:ms), legend_marker(Plots.default(:ms))...), #hide
) #hide
theme = my_theme #hide

"Marker style from the equilibrium characteristics: fill = arm wind, shape = zone, stroke = kite side." #hide
ss_style(q, vbp) = ( #hide
    c = is_down(q) ? :white : :black, #hide
    markerstrokecolor = is_down(q) ? :black : :white, #hide
    shape = is_center(q, vbp) ? :circle : :diamond, #hide
    ms = Plots.default(:ms), #hide
    markerstrokewidth = is_down(q) ? 8 : 10, #hide
) #hide

function plot_equilibria_phase_space(r) #hide
    d = r.NUM #hide
    plot(; xlims=(-180, 180), ylims=(-180, 180), aspect_ratio=:equal, #hide
        xticks=([-180, 0, 180], ["-180", "0", "180"]), #hide
        yticks=([-180, 0, 180], ["-180", "0", "180"]), #hide
        xlabel="α (deg)", ylabel="τ (deg)", grid=false, framestyle=:semi, legend=:outerright, colorbar=false) #hide
    vspan!([-180, -90]; theme.upwind..., label="Upwind arm") #hide
    vspan!([90, 180]; theme.upwind...) #hide
    vspan!([-90, 90]; theme.downwind...) #hide
    plot!(rad2deg.(d.line_x), rad2deg.(d.line_y); theme.streamlines...) #hide
    ax, ay = rad2deg.(d.ax), rad2deg.(d.ay) #hide
    L, W = rad2deg(0.2), rad2deg(0.06) #hide
    heads = [Shape([(s[1] + L / 2 * u, s[2] + L / 2 * v), #hide
        (s[1] - L / 2 * u - W * v, s[2] - L / 2 * v + W * u), #hide
        (s[1] - L / 2 * u + W * v, s[2] - L / 2 * v - W * u)]) #hide
        for (s, u, v) in zip(zip(ax, ay), d.au, d.av)] #hide
    plot!(heads; theme.heads...) #hide
    contour!(rad2deg.(d.αs), rad2deg.(d.τs), last.(d.dd); levels=[0], theme.nullcline_tau...) #hide
    contour!(rad2deg.(d.αs), rad2deg.(d.τs), first.(d.dd); levels=[0], theme.nullcline_alpha...) #hide
    plot!([NaN], [NaN]; theme.nullcline_tau..., label="D²(τ) = 0") #hide
    plot!([NaN], [NaN]; theme.nullcline_alpha..., label="D²(α) = 0") #hide
    styles = map(q -> ss_style(q, d.vbp), d.ss) #hide
    scatter!(rad2deg.(first.(d.ss)), rad2deg.(last.(d.ss)); #hide
        c=getproperty.(styles, :c), markerstrokecolor=getproperty.(styles, :markerstrokecolor), #hide
        shape=getproperty.(styles, :shape), ms=getproperty.(styles, :ms), #hide
        markerstrokewidth=getproperty.(styles, :markerstrokewidth)) #hide
    scatter!([NaN], [NaN]; theme.eq_legend..., label="Equilibrium") #hide
    plot!(widen=1.06) #hide
    plot!(leftmargin=6Plots.mm, bottommargin=5Plots.mm, rightmargin=5Plots.mm, topmargin=5Plots.mm, size=plot_size(1, 1.14)) #hide
    return current() #hide
end #hide

function plot_equilibria_eight_circle(r) #hide
    d = r.NUM #hide
    p_downwind = VIZ.plot_eight_circle(d.ss_downwind, d.vbp; plot_kwargs=(; color=:black, lw=2), #hide
        marker_kwargs=map(q -> ss_style(q, d.vbp), d.ss_downwind), project_eight=true) #hide
    plot!(title="Downwind arm", background_color_subplot=theme.downwind.c) #hide
    p_upwind = VIZ.plot_eight_circle(d.ss_upwind, d.vbp; plot_kwargs=(; color=:black, lw=2), #hide
        marker_kwargs=map(q -> ss_style(q, d.vbp), d.ss_upwind), project_eight=true) #hide
    plot!(title="Upwind arm", background_color_subplot=theme.upwind.c) #hide
    return plot(p_downwind, p_upwind; layout=(1, 2), size=plot_size(16 / 15, 2)) #hide
end #hide

function plot_eigvals(r) #hide
    d = r.NUM #hide
    hline([0], c=:black) #hide
    for (k, q) in enumerate(d.ss_) #hide
        local style = is_down(q) ? theme.downwind : theme.upwind #hide
        vspan!([k - 1/4, k + 1/4]; style...) #hide
        scatter!(fill(k, length(d.evs[k])), real.(d.evs[k]); ss_style(q, d.vbp)...) #hide
    end #hide
    annotate!([2.5, 2.5], [7.5, -1.5], ["Repulsive", "Attractive"]) #hide
    plot!(xlabel="Equilibrium", ylabel=raw"max Re(λ${}_i$)", ylim=(-3, 15), grid=false) #hide
    plot!(xticks=(1:4, ["" for _ in 1:4])) #hide
    plot!(leftmargin=5Plots.mm, bottommargin=2Plots.mm) #hide
    return current() #hide
end #hide

function plot_equilibrium_trajectories(r) #hide
    d = r.NUM #hide
    plot(; xlims=(-180 - rad2deg(0.2), 180 + rad2deg(0.2)), #hide
        ylims=(-180 - rad2deg(0.2), 180 + rad2deg(0.2)), aspect_ratio=:equal, #hide
        xticks=([-180, 0, 180], ["-180", "0", "180"]), #hide
        yticks=([-180, 0, 180], ["-180", "0", "180"]), #hide
        xlabel="α (deg)", ylabel="τ (deg)", grid=false, framestyle=:semi, legend=:outerright, colorbar=false) #hide
    vspan!([-180, -90]; theme.upwind..., label="Arm upwind") #hide
    vspan!([90, 180]; theme.upwind...) #hide
    vspan!([-90, 90]; theme.downwind..., label="Arm downwind") #hide
    contour!(rad2deg.(d.αs), rad2deg.(d.τs), last.(d.dd); levels=[0], color=:gray, alpha=0.4, linewidth=1.2, label="") #hide
    contour!(rad2deg.(d.αs), rad2deg.(d.τs), first.(d.dd); levels=[0], color=:gray, alpha=0.4, linewidth=1.2, label="") #hide
    styles = map(q -> ss_style(q, d.vbp), d.ss) #hide
    scatter!(rad2deg.(first.(d.ss)), rad2deg.(last.(d.ss)); #hide
        c=getproperty.(styles, :c), markerstrokecolor=getproperty.(styles, :markerstrokecolor), #hide
        shape=getproperty.(styles, :shape), ms=getproperty.(styles, :ms), #hide
        markerstrokewidth=getproperty.(styles, :markerstrokewidth)) #hide
    scatter!([NaN], [NaN]; theme.eq_legend..., label="Equilibrium") #hide
    for qs in d.traj_qs #hide
        plot!(first.(qs), last.(qs); c=:black, alpha=0.2, lw=1, label="") #hide
    end #hide
    plot!(leftmargin=6Plots.mm, bottommargin=5Plots.mm, rightmargin=5Plots.mm, topmargin=5Plots.mm, size=plot_size(1, 1.22)) #hide
    return current() #hide
end #hide

function plot_convergence_panel(c; legend=true) #hide
    (; dτ, distances, slope, i1) = c #hide
    plt = plot(; yscale=:log10, xlim=(0, 11), xticks=0:2:100, #hide
        ylim=(1e-16, 1e2), yticks=10.0 .^ (0:-3:-16), #hide
        xlabel="""Lemniscate revolutions, D(τ) $(dτ>0 ? ">" : "<") 0""") #hide
    plot!(i -> distances[i1] * exp10(-slope * (i - i1)); c=:black) #hide
    annotate!((5, 1e-11, (f"slope = {slope:.2f}", :right, ANNOTATIONFONTSIZE, colorant"black"))) #hide
    distance_style = (; c=dτ > 0 ? 1 : 2, ls=:dot, marker=:circle, legend_marker()...) #hide
    plot!(distances .+ 1e-16; distance_style...) #hide
    if legend #hide
        plot!([NaN]; distance_style..., label=raw"$|| x_i - x_\infty  ||$") #hide
    end #hide
    return plt #hide
end #hide

function plot_limit_cycle_convergence(r) #hide
    d = r.NUM #hide
    p_conv_p = plot_convergence_panel(d.conv.p; legend=false) #hide
    p_conv_m = plot_convergence_panel(d.conv.m; legend=true) #hide
    plot!(p_conv_p, leftmargin=11Plots.mm) #hide
    plot!(p_conv_m, leftmargin=5Plots.mm) #hide
    plot!(p_conv_p, ylabel="Distance to limit") #hide
    conv_two = plot(p_conv_p, p_conv_m; layout=(1, 2), size=SIZE_FULL) #hide
    plot!(bottommargin=7Plots.mm) #hide
    return conv_two #hide
end #hide

function plot_yz_limit_cycles(r) #hide
    d = r.NUM #hide
    plot(; aspect_ratio=:equal, xlabel="y (m)", ylabel="z (m)", xtick=-40:10:40, ytick=10:5:20) #hide
    plot!(invert(yz_split(d.lc_m, d.n_splits))..., lw=4, c=2, arrows=true) #hide
    plot!(invert(yz_split(d.lc_p, d.n_splits))..., lw=2, c=1, arrows=true) #hide
    plot!([NaN], [NaN], c=1, label="D(τ) > 0") #hide
    plot!([NaN], [NaN], c=2, label="D(τ) < 0") #hide
    plot!([0], [26]) #hide
    c, δ = d.c, d.δ #hide
    lens!([c[1] - δ, c[1] + δ], [c[2] - δ, c[2] + δ], #hide
        inset=(1, bbox(0.3, 0.05, 0.3, 0.4)), lw=1.5, lc=:grey) #hide
    l = plot!().inset_subplots[1] #hide
    xticks_ = (c[1]-δ):δ:(c[1]+δ) #hide
    yticks_ = (c[2]-δ):δ:(c[2]+δ) #hide
    xticks!(l, xticks_, [f"{t:.1f}" for t in xticks_]) #hide
    yticks!(l, yticks_, [f"{t:.1f}" for t in yticks_]) #hide
    plot!(l, tickfontsize=TICKFONTSIZE - 1) #hide
    plot!(; leftmargin=5Plots.mm, bottommargin=5Plots.mm, size=plot_size(2 / 3, 1)) #hide
    return current() #hide
end #hide

save_fig(plot_equilibria_phase_space(results), "equilibria_phase_space") #hide
save_fig(plot_equilibria_eight_circle(results), "equilibria_eight_circle") #hide
save_fig(plot_eigvals(results), "_equilibria_eigvals") #hide
save_fig(plot_equilibrium_trajectories(results), "_equilibria_trajectories") #hide
save_fig(plot_time_series(results.NUM.sol), "time_series_4D") #hide
save_fig(plot_limit_cycle_convergence(results), "limit_cycle_convergence") #hide
save_fig(plot_yz_limit_cycles(results), "yz_limit_cycles") #hide

Markdown.parse("""
**Figure data — equilibria.** There are $(length(results.NUM.ss)) equilibria, all unstable: the largest real eigenvalue lies in \$[$(round(minimum(maximum(real, e) for e in results.NUM.evs); digits=2)), $(round(maximum(maximum(real, e) for e in results.NUM.evs); digits=2))]\\,\\text{s}^{-1}\$. The arm and lines are coplanar at every equilibrium ($(results.NUM.all_aligned)). The two top views separate $(length(results.NUM.ss_downwind)) downwind from $(length(results.NUM.ss_upwind)) upwind equilibria.

**Figure data — limit cycles.** The two attracting cycles carry \$P_+ = $(round(results.NUM.cycles[1].λ; digits=3))\$ (Floquet multiplier) and \$P_- = $(round(results.NUM.cycles[2].λ; digits=3))\$. Their convergence slopes are $(results.NUM.cycle_lines[1]) and $(results.NUM.cycle_lines[2]). A deflated-shooting inventory over $(BK_NSEEDS) Halton seeds finds $(length(results.NUM.found)) additional cycles; quasi-Monte-Carlo routes $(results.NUM.qmc[1]) / $(BK_NSEEDS) seeds to \$D(\\tau)>0\$, $(results.NUM.qmc[2]) / $(BK_NSEEDS) to \$D(\\tau)<0\$, and $(results.NUM.qmc[3]) elsewhere.

**Figure data — trajectories near equilibria.** Each panel seeds $(results.NUM.N_eq) trajectories on a circle of radius \$\\pi/10\$ at rest around one equilibrium and integrates them for $(results.NUM.t_eq) s; the translucent black curves show the local flow directly.

**Figure data — convergence.** Averaged over the two directions, the Poincaré map loses about two digits of distance per revolution until the solver tolerance \$10^{-13}\$ is reached. The slope printed in each panel is the quantity to quote in prose; the dashed perturbed cloud is deliberately absent here (\$n\\text{trace}=0\$).

**Prose hint:** the two cycles are distinct in the \$YZ\$ projection (zoom inset), so the phase of \$\\tau\$ is physically meaningful; the deflation plus QMC inventory supports the claim that only two attracting cycles exist.

**Open questions from the JTCAM draft:** confirm the colour-code wording in the phase-space caption (max \$\\mathrm{Re}(\\lambda)\\approx 3\$ for the circle and \$\\approx 12\$ for the diamonds); decide whether to keep the coplanarity remark in the main text or the appendix.

![Equilibria in phase space](figs/equilibria_phase_space.svg)

![Equilibria in the eight/circle view](figs/equilibria_eight_circle.svg)

![Eigenvalues at representative equilibria](figs/_equilibria_eigvals.svg)

![Trajectories near equilibria](figs/_equilibria_trajectories.svg)

![Transient settling onto the pumping cycle](figs/time_series_4D.svg)

![Convergence on the Poincaré section](figs/limit_cycle_convergence.svg)

![Two limit cycles in the YZ plane](figs/yz_limit_cycles.svg)
""") #hide

# ## 3. Limit-Cycle Optimization and Robustness
#
# We maximize $\\bar P = W(t_f)/t_f$ subject to periodic closure over $\\Theta = (\\ell, I, b, \\Delta\\varphi)$ at $w_0 = 9$ m/s with $\\Delta\\theta = 6^\\circ$ fixed. The reference is the kinetic PTO point ($b = I\\omega$), then $\\Theta$ is optimized. Every panel below comes from the cached `results.OPT` bundle.

function plot_optimization_timeseries(r, which) #hide
    o = r.OPT #hide
    return plot_time_series(o.sols[which == :init ? 1 : 2]) #hide
end #hide

function plot_landscape(r) #hide
    o = r.OPT #hide
    land, guides = o.LANDSCAPE, GUIDES #hide
    (; xs, syms, ys, window, box, center) = land #hide
    lo, hi = 1 - window, 1 + window #hide
    panels = [ #hide
        begin #hide
            x0 = TEX[syms[i]][3](center[i]) #hide
            ticks = ([x0 * lo, x0, x0 * hi], [f"-{100window:.0f}%", f"{x0:.2f}", f"+{100window:.0f}%"]) #hide
            plt = plot(; xlims=(x0 * lo, x0 * hi), ylims=(guides.ymin, 1), xticks=ticks, #hide
                title="$(LABELS[sym])", ylabel=(i == 1 ? "P / P*" : ""), #hide
                yticks=(i == 1 ? ([guides.ymin, 1], [plain_tick(guides.ymin), "1"]) : ([guides.ymin, 1], ["", ""])), #hide
                leftmargin=(i == 1 ? 18 : 3)Plots.mm, rightmargin=4Plots.mm, #hide
                bottommargin=5Plots.mm, topmargin=2.5Plots.mm) #hide
            interp(x) = begin #hide
                k = clamp(searchsortedlast(xs, x), 1, length(xs) - 1) #hide
                ys[i][k] + (x - xs[k]) / (xs[k+1] - xs[k]) * (ys[i][k+1] - ys[i][k]) #hide
            end #hide
            piece(a, b) = vcat(a, xs[(a .< xs) .& (xs .< b)], b) #hide
            xa, xb = max(box[i][1], lo), min(box[i][2], hi) #hide
            xa < xb && plot!(plt, x0 .* piece(xa, xb), interp.(piece(xa, xb)); c=i, ls=:solid) #hide
            for (a, b) in ((lo, min(xa, hi)), (max(xb, lo), hi)) #hide
                a < b && plot!(plt, x0 .* piece(a, b), interp.(piece(a, b)); c=i, ls=:dot) #hide
            end #hide
            for edge in box[i], k in (-1, 1) #hide
                vline!(plt, [x0 * (edge + k * window / 40)]; c=:black, lw=1.2) #hide
            end #hide
            hline!([guides.ymin]; guides.floor...) #hide
            vline!([x0]; guides.line...) #hide
            series = plt[1].series_list #hide
            pushfirst!(series, pop!(series), pop!(series)) #hide
            plt #hide
        end for (i, sym) in enumerate(syms) #hide
    ] #hide
    return plot(panels...; layout=(1, length(syms)), size=plot_size(1/2, 2), framestyle=:grid, grid=false) #hide
end #hide

function plot_dtheta(r) #hide
    o = r.OPT #hide
    ang = [d.Δθ for d in o.dtheta] #hide
    i0 = argmin(abs.(ang .- rad2deg(o.p_ref.Δθ))) #hide
    style = (; c=:black) #hide
    Pn = [d.P for d in o.dtheta] ./ o.dtheta[i0].P #hide
    Pn_opt = [d.P for d in o.dtheta_opt] ./ o.dtheta[i0].P #hide
    Tmin0 = [d.Tmin / 1000 for d in o.dtheta] #hide
    Tmax0 = [d.Tmax / 1000 for d in o.dtheta] #hide
    Tmin1 = [d.Tmin / 1000 for d in o.dtheta_opt] #hide
    Tmax1 = [d.Tmax / 1000 for d in o.dtheta_opt] #hide
    f1 = plot(ang, Pn; style..., seriestype=:scatter, ms=4, legend_marker(4)..., #hide
        label="Θ fixed", ylabel="P / P*", xlabel="Δθ (deg)", #hide
        ylims=(0, 1.03max(1, maximum(Pn_opt))), #hide
        yticks=([0, 0.5, 1], ["0", "0.5", "1"]), legend=:topright) #hide
    plot!(f1, ang, Pn_opt; style..., marker=:none, label="Θ re-optimized") #hide
    f2 = plot(ang, Tmax0; style..., seriestype=:scatter, ms=4, label="", #hide
        ylabel="Line tension (kN)", xlabel="Δθ (deg)") #hide
    scatter!(f2, ang, Tmin0; style..., ms=4, label="") #hide
    plot!(f2, ang, Tmax1; style..., label="", fillrange=Tmin1, fillalpha=0.1) #hide
    plot!(f2, ang, Tmin1; style..., label="") #hide
    foreach(f -> vline!(f, [rad2deg(o.p_ref.Δθ)]; GUIDES.line...), (f1, f2)) #hide
    foreach(f -> plot!(f, xticks=2:2:20, xlim=extrema(ang)), (f1, f2)) #hide
    for f in (f1, f2) #hide
        pushfirst!(f[1].series_list, pop!(f[1].series_list)) #hide
    end #hide
    foreach(f -> vline!(f, [rad2deg(o.p_ref.Δθ)]; c=:black, lw=0.8, alpha=0.7), (f1, f2)) #hide
    plot!(f1, leftmargin=18Plots.mm, bottommargin=7Plots.mm) #hide
    plot!(f2, leftmargin=6Plots.mm, bottommargin=7Plots.mm) #hide
    return plot(f1, f2; layout=(1, 2), size=plot_size(1, 2)) #hide
end #hide

function plot_shape_grid(r) #hide
    o = r.OPT #hide
    function pane(Z, title) #hide
        plt = heatmap(o.x_map, o.shape_dphi, Z; c=:viridis, title=title, titlefontsize=TICKFONTSIZE, #hide
            xticks=-20:10:20, yticks=10:10:90, xlim=extrema(o.x_map), ylim=extrema(o.shape_dphi), #hide
            xlabel="Δθ deviation (%)", ylabel="Δφ (deg)") #hide
        vline!(plt, [0]; c=:white, ls=:dash, lw=1) #hide
        hline!(plt, [rad2deg(o.p_ref.Δφ)]; c=:white, ls=:dash, lw=1) #hide
        scatter!(plt, [0], [rad2deg(o.p_ref.Δφ)]; c=:white, markerstrokecolor=:black, ms=5) #hide
        return plt #hide
    end #hide
    p1 = pane(o.P_map, "P (kW)") #hide
    p2 = pane(o.Tmax_map, "Tmax (kN)") #hide
    p3 = pane(o.amax_map, "max a (m/s²)") #hide
    plot!(p2, ylabel="") #hide
    plot!(p3, ylabel="") #hide
    return plot(p1, p2, p3; layout=(1, 3), size=plot_size(4 / 3, 2), #hide
        leftmargin=6Plots.mm, bottommargin=5Plots.mm, topmargin=4Plots.mm) #hide
end #hide

function plot_wind_sweep(r) #hide
    o = r.OPT #hide
    (; syms, no_opt, opt, theta) = o.sweep #hide
    winds, i_ref = o.winds, o.i_ref #hide
    power_vec = [s.P for s in opt] #hide
    min_tension_vec = [s.Tmin for s in opt] #hide
    max_tension_vec = [s.Tmax for s in opt] #hide
    power0 = [s.P for s in no_opt] #hide
    min_tension0 = [s.Tmin for s in no_opt] #hide
    max_tension0 = [s.Tmax for s in no_opt] #hide
    theta_norm = [maxnormalize(θ) for θ in invert(theta)] #hide
    ref = "Fixed at $(round(Int, WIND_REF)) m/s" #hide
    ref_style = (; ls=:dot, c=:black) #hide
    f1 = plot(ylabel="Power (kW)", legend=:topleft, xformatter=:plain, xlabel="Wind speed (m/s)", #hide
        ylims=(0, 60), yticks=0:20:100, widen=true) #hide
    plot!(f1, winds, power0 ./ 1000; ref_style...) #hide
    plot!(f1, winds, power_vec ./ 1000; c=1) #hide
    f2 = plot(ylabel="Line tension (kN)", legend=:topleft, xformatter=:plain, xlabel="Wind speed (m/s)", #hide
        ylims=(0, 40), yticks=0:20:100, widen=true) #hide
    plot!(f2, winds, [max_tension0 ./ 1000, min_tension0 ./ 1000]; ref_style...) #hide
    plot!(f2, winds, max_tension_vec ./ 1000; fillrange=min_tension_vec ./ 1000, c=1, fillalpha=0.1) #hide
    plot!(f2, winds, min_tension_vec ./ 1000; c=1) #hide
    n = length(syms) #hide
    f3 = plot(yticks=((0, 1), ("0", "max")), ylim=(0, 1), ylabel="Parameters", xformatter=:plain, #hide
        xlabel="Wind speed (m/s)", legend=false, widen=true) #hide
    plot!(f3, winds, theta_norm; c=permutedims(collect(1:n))) #hide
    for f in (f1, f2, f3) #hide
        vline!(f, [WIND_REF]; GUIDES.line...) #hide
        pushfirst!(f[1].series_list, pop!(f[1].series_list)) #hide
    end #hide
    names = Dict(:l => "Line length", :I => "Arm inertia", :b => "Arm braking coefficient", :Δφ => "Lemniscate width") #hide
    ib = findfirst(==(:b), syms) #hide
    annotate!(f3, first(winds) + 0.25, theta_norm[ib][1] - 0.14, #hide
        (names[:b], ANNOTATIONFONTSIZE, :left, PALETTE[ib])) #hide
    y = 0.62 #hide
    for i in sortperm([theta_norm[i][end] for i in eachindex(syms)]; rev=true) #hide
        i == ib && continue #hide
        annotate!(f3, last(winds) - 0.25, y, (names[syms[i]], ANNOTATIONFONTSIZE, :right, PALETTE[i])) #hide
        y -= 0.10 #hide
    end #hide
    fig = plot(f1, f2, f3; layout=(3, 1), size=plot_size(2, 1), leftmargin=10Plots.mm, bottommargin=4Plots.mm) #hide
    foreach(i -> plot!(fig, subplot=i, xformatter=x -> "", xlabel=""), (1, 2)) #hide
    return (; fig, panels=(f1, f2, f3)) #hide
end #hide

save_fig(plot_optimization_timeseries(results, :init), "optimization_timeseries_init") #hide
save_fig(plot_optimization_timeseries(results, :opt), "optimization_timeseries_opt") #hide
save_fig(plot_time_series_grid(results.OPT.sols; titles=("PTO optimum", "Dynamic optimum")), "optimization_timeseries_all") #hide
save_fig(plot_landscape(results), "optimization_landscape") #hide
save_fig(plot_dtheta(results), "optimization_dtheta") #hide
save_fig(plot_shape_grid(results), "optimization_shape_grid") #hide
begin #hide
    (; fig, panels) = plot_wind_sweep(results) #hide
    save_fig(fig, "optimization_wind_sweep_all") #hide
    for (f, panel) in zip(panels, ("power", "tension", "params")) #hide
        plot!(f, size=plot_size(1, 1), leftmargin=10Plots.mm, bottommargin=6Plots.mm) #hide
        save_fig(f, "optimization_wind_sweep_$(panel)") #hide
    end #hide
end #hide

Markdown.parse("""
**Figure data — optimum.** The reference (PTO) point gives $(round(results.OPT.opt_ref.P / 1000; digits=2)) kW; the \$\\Theta\$ optimum gives $(round(maximum(s -> s.P, results.OPT.sweep.opt) / 1000; digits=2)) kW at the reference wind. The achieved optimum is \$\\ell = $(round(results.OPT.p_ref.l; digits=2))\$ m, \$I = $(round(results.OPT.p_ref.I; digits=1))\$ kg m², \$b = $(round(results.OPT.p_ref.b; digits=1))\$ N m s rad⁻¹, \$\\Delta\\varphi = $(round(rad2deg(results.OPT.p_ref.Δφ); digits=2))^\\circ\$. Ipopt converged in $(results.OPT.opt_ref.iter) iterations in $(round(results.OPT.opt_ref.t_wall; digits=2)) s.

**Figure data — landscape.** \$\\ell\$ is the stiffest direction (\$(-11.5\\%, -8.1\\%)\$ at \$\\pm 10\\%\$), \$I\$ is about four times flatter (\$-2.6\\%\$) and \$b\$ more than twenty times flatter (\$-0.4\\%\$). The \$\\Delta\\varphi\$ optimum (\$16.6^\\circ\$) is interior: \$-5.4\\%\$ at \$15.0^\\circ\$, \$-2.9\\%\$ at \$18.3^\\circ\$. The box lower bounds on \$\\ell\$, \$I\$ and \$b\$ are far from the optimum, so none binds.

**Figure data — Δθ sweep.** At fixed \$\\Theta\$ the power falls from \$12.25\$ kW at \$6^\\circ\$ to \$3.96\$ kW at \$2^\\circ\$ and \$0.75\$ kW at \$20^\\circ\$, while \$T_{\\max}\$ rises to \$13.8\$ kN at both ends. Re-optimizing \$\\Theta\$ recovers \$12.82\$ kW at \$2^\\circ\$ and \$8.27\$ kW at \$20^\\circ\$, with flatter loads (\$T_{\\max} = 9.9\$–\$10.8\$ kN). The re-optimized power peaks at \$2\$–\$4^\\circ\$, \$4.6\\%\$ above the nominal point.

**Figure data — shape map.** At fixed optimal \$\\ell, I, b\$ the map peaks at \$\\Delta\\theta = -20\\%\$, \$\\Delta\\varphi = 18.9^\\circ\$ with \$12.46\$ kW and collapses to \$2\$–\$3\$ kW for \$\\Delta\\varphi \\ge 36.7^\\circ\$; \$\\Delta\\varphi = 90^\\circ\$ has no pumping cycle at any \$\\Delta\\theta\$.

**Figure data — wind sweep.** Re-optimizing gains \$+13.5\\%\$ at \$5\$ m/s and \$+8.6\\%\$ at \$15\$ m/s. \$\\ell\$ shortens \$30.3 \\to 26.1\$ m and \$I\$ falls \$451 \\to 326\$ kg m² with wind, while \$b\$ rises \$498 \\to 1151\$ and \$\\Delta\\varphi\$ closes \$17.6^\\circ \\to 16.3^\\circ\$. \$T_{\\min}\$ stays positive throughout, confirming the taut-line assumption a posteriori.

**Prose hint:** these panels support the ECC conclusions — roughly cubic power scaling, quadratic tension scaling and near-linear braking coefficient in wind speed, with little loss when the \$9\$ m/s design is held fixed.

**Open questions from the JTCAM draft:** state explicitly which variables sit on their box bounds in the landscape; decide whether to quote the \$13\\%\$ optimization gain from the finished ECC paper or the value computed here; decide whether the Δθ panel belongs in the main text or the appendix.

![Optimization: initial cycle](figs/optimization_timeseries_init.svg)

![Optimization: optimized cycle](figs/optimization_timeseries_opt.svg)

![PTO vs dynamic optimum](figs/optimization_timeseries_all.svg)

![Objective landscape](figs/optimization_landscape.svg)

![Δθ sensitivity](figs/optimization_dtheta.svg)

![10×10 shape map](figs/optimization_shape_grid.svg)

![Wind sweep](figs/optimization_wind_sweep_all.svg)

![Wind sweep: power](figs/optimization_wind_sweep_power.svg)

![Wind sweep: tension](figs/optimization_wind_sweep_tension.svg)

![Wind sweep: parameters](figs/optimization_wind_sweep_params.svg)
""") #hide

# ## 4. Conclusion
#
# The reduced 4D model reproduces the 10D dynamics without drift and at lower cost; it has two attracting limit cycles; optimizing the design parameters raises average power substantially and the $9$ m/s optimum is robust over $5$–$15$ m/s with positive line tension. The figures above carry the numbers needed to write the corresponding paragraphs.

# ## Appendix A. Validation of the 4D reduction
#
# The original 10D DAE drifts when integrated without a manifold-projection callback. The callback removes the drift at a runtime cost; the 4D model has no algebraic constraint at all.

function plot_trajectory_10D(r) #hide
    d = r.D10 #hide
    fig_traj_short = VIZ.plot_trajectory_10D(d.sol10; tspan=(0, d.tf_short), points_per_second=d.points_per_second) #hide
    fig_traj = VIZ.plot_trajectory_10D(d.sol10; points_per_second=d.points_per_second) #hide
    fig_traj_cb = VIZ.plot_trajectory_10D(d.sol10_cb; points_per_second=d.points_per_second) #hide
    params = (; xlabel="x (m)", ylabel="y (m)", zlabel="z (m)", legend=false, widen=false) #hide
    plot!(fig_traj_short; title="tf = 20 s", params..., legend=:topleft) #hide
    plot!(fig_traj; title="tf = 100 s", params...) #hide
    plot!(fig_traj_cb; title="tf = 100 s + callback", params...) #hide
    tf_fig = plot(fig_traj_short, fig_traj, fig_traj_cb; layout=(1, 3), size=SIZE_FULL) #hide
    for i in 1:3 #hide
        plot!(tf_fig, subplot=i, xticks=edge_ticks(tf_fig[i], :x; include_zero=true), #hide
            yticks=edge_ticks(tf_fig[i], :y; include_zero=true), #hide
            zticks=edge_ticks(tf_fig[i], :z; include_zero=true)) #hide
    end #hide
    return tf_fig #hide
end #hide

function plot_residuals_10D(r) #hide
    d = r.D10 #hide
    plot(; xscale=:log10, xticks=d.xticks_, xformatter=x->"", ylim=d.ylim, yscale=:log10, #hide
        yticks=d.yticks_, ylabel="Residual on position (m)", legend=:topleft) #hide
    plot!(d.ts, norm.(getindex.(d.res, Ref(1:3))), label="Without callback") #hide
    fig_pos = plot!(d.ts, norm.(getindex.(d.res_cb, Ref(1:3))), label="With callback") #hide
    plot(; xscale=:log10, xticks=d.xticks_, xlabel="Time (s)", ylim=d.ylim, yscale=:log10, #hide
        yticks=d.yticks_, ylabel="Residual on speed (m/s)") #hide
    plot!(d.ts, norm.(getindex.(d.res, Ref(4:6)))) #hide
    fig_speed = plot!(d.ts, norm.(getindex.(d.res_cb, Ref(4:6)))) #hide
    fig = plot(fig_pos, fig_speed; layout=(2, 1), size=SIZE_HALF .* (1, 1.5)) #hide
    plot!(bottommargin=3Plots.mm, leftmargin=6Plots.mm) #hide
    return fig #hide
end #hide

function plot_three_models_convergence(r) #hide
    d = r.D10 #hide
    plot(; xlabel="Solver tolerance", ylabel="Norm of error (m)", legend=:bottomright) #hide
    plot!(d.tols, d.err10; marker=:circle, legend_marker()..., label="10D") #hide
    plot!(d.tols, d.err10cb; marker=:square, legend_marker()..., label="10D + callback") #hide
    plot!(d.tols, d.err4; marker=:diamond, legend_marker()..., label="4D") #hide
    plot!(xscale=:log10, xticks=exp10.(-15:3:15), yscale=:log10, yticks=exp10.(-15:3:15)) #hide
    fig = plot!(xflip=true, ymirror=true, legend=:topright, size=SIZE_HALF) #hide
    plot!(bottommargin=4Plots.mm, rightmargin=6Plots.mm) #hide
    return fig #hide
end #hide

save_fig(plot_trajectory_10D(results), "trajectory_10D_short_long_cb") #hide
save_fig(plot_residuals_10D(results), "residuals_10D_interp") #hide
save_fig(plot_three_models_convergence(results), "three_models_convergence") #hide

Markdown.parse("""
**Figure data — 10D drift and callback.** At \$t_f = $(results.D10.tf)\$ s and tolerance \$10^{-3}\$, the raw 10D solution drifts; the projection callback keeps the manifold residual near the callback tolerance \$10^{-10}\$. The callback costs $(round(results.D10.ratio_cb; digits=2))× the runtime of raw 10D.

**Figure data — convergence.** Against a 4D reference at tolerance \$10^{-14}\$, all three formulations converge to the same trajectory as the tolerance decreases from \$10^{-2}\$ to \$10^{-12}\$; the 4D model has no drift by construction. At the production tolerance the 4D run takes $(round(1e3 * results.D10.btime4; digits=2)) ms versus $(round(1e3 * results.D10.btime; digits=2)) ms for raw 10D ($(round(results.D10.speedup; digits=1))× slower) and $(round(1e3 * results.D10.btime_cb; digits=2)) ms for 10D with callback ($(round(results.D10.speedup_cb; digits=1))× slower).

**Prose hint:** the horizontal axis is reversed so that tighter tolerances progress to the right; quote the drift number and the speedup together, as in the ECC paper's model-validation paragraph.

**Open question from the JTCAM draft:** consider dropping the left subplot of the trajectory figure (it duplicates the long-horizon panel); it is kept here for parity with the legacy figure set.

![10D trajectories, short/long with callback](figs/trajectory_10D_short_long_cb.svg)

![Interpolated 10D manifold residuals](figs/residuals_10D_interp.svg)

![Convergence of the three formulations](figs/three_models_convergence.svg)
""") #hide
