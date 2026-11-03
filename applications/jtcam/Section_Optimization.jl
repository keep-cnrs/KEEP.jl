# Optimization section.
#
# Two independent studies, each guarded by a flag at the top of the file:
#  - PARAMETRIC: for each wind speed, maximize the average power over one of
#    the blocks A = (ℓ, I, b), B = (Δθ, Δφ) or (A, B), and compare it with the
#    joint 10 m/s optimum held fixed.
#  - SHAPE: does the optimum of the joint (A, B) optimization depend on what was
#    optimized before it?  A → (A, B), B → (A, B), and (A, B) alone, all started
#    from the same parameters; runtime and ipopt iterations are reported for
#    each path.
#
# Numerical sections and plot sections are separated by `##` markers: run the
# file once, then re-run a plot section to iterate on a figure without
# recomputing the optimizations.
#
# Run from the repo root:
#   julia --project=applications/jtcam applications/jtcam/Section_Optimization.jl
using Pkg
Pkg.activate(@__DIR__)

using KEEP: PointMass4 as PM4
using KEEP: PointMassPara as PMP
using KEEP: LimitCycle as LC
using KEEP: Optimization as OPT

using Plots
using PyFormattedStrings: @f_str
using SplitApplyCombine: invert
using Logging: with_logger, NullLogger
using LinearAlgebra: norm

include(joinpath(@__DIR__, "jtcam_theme.jl"))

## Settings
const PARAMETRIC = true  # wind sweep without and with re-optimization
const SHAPE = true       # do the optimization paths over A and B agree?
const TOL = 1e-5
const MAX_WALL_TIME = 120.0
const WIND_REF = 10.0
const SENSE = +
const SYMS_A = [:r, :I_eq, :torque_slope]  # ℓ, I, b
const SYMS_B = [:Δθ, :Δφ]                  # lemniscate half-height, half-width
const SYMS_AB = vcat(SYMS_A, SYMS_B)       # the joint optimization

# Parameter names and units for the axis titles and the legends
const LABELS = Dict(:r => "ℓ", :I_eq => "I", :torque_slope => "b",
    :Δθ => "Δθ", :Δφ => "Δφ")
const UNITS = Dict(:r => "m", :I_eq => "kg m²", :torque_slope => "kg m² s rad⁻¹",
    :Δθ => "rad", :Δφ => "rad")

## Helpers
"""
Min and max line tension over one lobe of a limit cycle. The tension is
symmetrical over a revolution, so one half-period holds both extrema; a fine
sampling of the (smooth) tension is as accurate as the global search it
replaces, for a fraction of the cost.
"""
function extrema_tension(sol)
    tmin, tmax = sol.prob.tspan
    ts = range(tmin, tmin + (tmax - tmin) / 2, length=2001)
    return extrema(PM4.compute_line_tension.(sol.(ts), Ref(sol.prob.p)))
end

function trajectory_stats(sol)
    min_tension, max_tension = extrema_tension(sol)
    power = sol.u[end][5] / sol.t[end]
    return power, min_tension, max_tension, sol.t[end]
end

# Limit cycle, optimization, and limit cycle of the optimized system, for one
# value of one parameter. `lower`/`upper` are the box constraints of `syms`.
function compute_results(param, value; syms, p0, lower, upper, sense, tol)
    p0 = copy(p0)
    p0[param] = value

    # Compute limit cycle and corresponding stats
    vbp0 = PMP.build_vbpara(p0)
    lc = LC.compute_limit_cycle(vbp0; sense=sense, save_everystep=true, tol=tol)
    stats_no_opt = trajectory_stats(lc)

    # Optimize the selected parameters
    _, stats, model = with_logger(NullLogger()) do
        OPT.optimize(p0, syms, lower, upper; sense, tol, max_wall_time=MAX_WALL_TIME)
    end
    para_dims, _, iterate_dims, _ = OPT.compute_dims(p0, syms)
    solution = stats.solution  # dimensionless
    theta = solution[6:end] .* para_dims

    # Compute limit cycle and corresponding stats for the optimized parameters
    vbp = copy(vbp0)
    vbp[syms] .= solution[(end-length(syms)+1):end]
    lc_opt = LC.compute_limit_cycle(vbp; sense=sense, save_everystep=true, tol=tol)
    stats_opt = trajectory_stats(lc_opt)

    # Error between the limit cycle found by ipopt and by our own method
    (α0, dα0, dτ0), tf, P = lc_opt.u[1][[1, 3, 4]], lc_opt.t[end], lc_opt.u[end][5] / lc_opt.t[end]
    println("$param = $value
    opt error = $(norm(model.c!(model.meta.ucon, solution)))
    LC error =  $(norm((solution .* iterate_dims)[1:5] - [α0, dα0, dτ0, tf, P]))")

    return stats_no_opt, stats_opt, theta
end

function range_optim(param, values; syms, p_ref, lower, upper, sense, tol)
    return invert([
        compute_results(param, v; syms, p0=p_ref, lower, upper, sense, tol)
        for v in values
    ])
end

"""
Wind sweep for one parameter block: hold `p_ref` fixed, or re-optimize `syms`
inside `lower`/`upper` at each wind speed. Returns the block, the statistics
with and without re-optimization, and the optimized parameters.
"""
function wind_sweep(syms, p_ref, winds, lower, upper)
    no_opt, opt, theta = range_optim(:v_ref, winds;
        syms, p_ref, lower, upper, sense=SENSE, tol=TOL)
    return (; syms, no_opt, opt, theta)
end

maxnormalize(x) = x ./ maximum(x)

"""
One optimization step: optimize `syms` from `p` and record the new parameters,
the optimal power and the cost of the solve (`stats.elapsed_time` is ipopt's own
timer, `t_wall` the wall-clock of the whole call including its initial guess).
"""
function opt_step(p, syms, lower, upper)
    solve() = with_logger(NullLogger()) do
        OPT.optimize(p, syms, lower, upper; sense=SENSE, tol=TOL, max_wall_time=MAX_WALL_TIME)
    end
    # The first solve with a given variable layout pays Julia's compilation of the
    # dual-number Jacobian inside ipopt's wall clock (40-245 s here), which eats
    # the time cap and truncates the solve. Warm the layout up and keep the second
    # run. The real cap is required: a short one times out before the Jacobian is
    # compiled, so the second run pays the compilation and is truncated in turn.
    solve()
    t_wall = time()
    solution, stats, model = solve()
    stats.status == :first_order || @warn "opt_step($(join(syms, ", "))) did not converge" stats.status stats.iter
    p_next = copy(p)
    p_next[syms] .= solution.params[syms]  # physical units

    return (; p=p_next, P=solution.P, iter=stats.iter, t_solve=stats.elapsed_time,
        t_wall=time() - t_wall, residual=norm(model.c!(model.meta.ucon, stats.solution)))
end

"""
Objective landscape around the optimal point `p`, one parameter at a time:
scale one parameter by each factor `x ∈ [1 - window, 1 + window]`, keep the
other parameters at their optimum, and take the power of the limit cycle at
those fixed parameters — the objective the optimizer maximizes once the cycle
is solved. Powers are normalized by their value at `x = 1`, so the optimum
plots at `(1, 1)`. `lower` and `upper` are the box constraints of `syms`
(physical units), returned as normalized factors so the caller can mark them.
Returns the factors, the parameter names, one row of normalized powers and one
`(low, up)` pair of normalized box bounds per parameter.
"""
function landscape(p, syms; lower, upper, window=0.1, npoints=21, tol=TOL)
    xs = range(1 - window, 1 + window, length=npoints)
    power(sym, x) = begin
        q = copy(p)
        q[sym] = x * p[sym]
        lc = LC.compute_limit_cycle(PMP.build_vbpara(q); sense=SENSE, tol=tol)
        return lc.u[end][5] / lc.t[end]
    end
    vals = [power(sym, x) for sym in syms, x in xs]
    k0 = findmin(abs.(xs .- 1))[2]  # the optimum, where the landscape is normalized
    return (; xs=collect(xs), syms=collect(syms), window,
        ys=[row ./ row[k0] for row in eachrow(vals)],
        center=collect(p[syms]),  # the optimum in physical units, for the x axis
        box=[(lo / p[sym], up / p[sym]) for (sym, lo, up) in zip(syms, lower, upper)])
end

"Total power, wall-clock time and ipopt iterations over consecutive optimization steps."
function path_stats(steps...)
    return (; P=last(steps).P, t_wall=sum(s -> s.t_wall, steps),
        t_solve=sum(s -> s.t_solve, steps), iter=sum(s -> s.iter, steps))
end

## Numerical: 10 m/s reference, one optimization per parameter block
if PARAMETRIC || SHAPE
    p_start = PMP.build_para(v_ref=WIND_REF)
    lower_A, upper_A = OPT.make_bounds(p_start, SYMS_A)
    lower_B, upper_B = OPT.make_bounds(p_start, SYMS_B; multiplier=1.5)  # ×3 on Δφ leaves the window
    lower_AB, upper_AB = vcat(lower_A, lower_B), vcat(upper_A, upper_B)

    step_A = opt_step(p_start, SYMS_A, lower_A, upper_A)
    step_B = opt_step(p_start, SYMS_B, lower_B, upper_B)
    step_AB = opt_step(p_start, SYMS_AB, lower_AB, upper_AB)  # the joint optimum
end

## Numerical: wind sweep over A, B and (A, B), from the joint 10 m/s optimum
if PARAMETRIC
    winds = range(5, 15, length=11)
    i_ref = findfirst(x -> x >= WIND_REF, winds)
    p_ref = step_AB.p  # the joint optimum at the reference wind speed
    # the boxes of the 10 m/s study, so the sweep moves inside the same region
    sweeps = (; A=wind_sweep(SYMS_A, p_ref, winds, lower_A, upper_A),
        B=wind_sweep(SYMS_B, p_ref, winds, lower_B, upper_B),
        AB=wind_sweep(SYMS_AB, p_ref, winds, lower_AB, upper_AB))
end

## Numerical: parameter paths A = (r, I_eq, torque_slope), B = (Δθ, Δφ)
if SHAPE
    # The second step is always the joint (A, B) optimization: the three paths
    # differ only by the block optimized before it, or by none at all.
    step_A_AB = opt_step(step_A.p, SYMS_AB, lower_AB, upper_AB)  # A, then (A, B)
    step_B_AB = opt_step(step_B.p, SYMS_AB, lower_AB, upper_AB)  # B, then (A, B)

    P_start = trajectory_stats(LC.compute_limit_cycle(PMP.build_vbpara(p_start);
        sense=SENSE, tol=TOL, save_everystep=true))[1]

    # Path totals and final parameter vectors
    x_path(p) = vcat(collect(p[SYMS_A]), collect(p[SYMS_B]))
    xs = (; A_AB=x_path(step_A_AB.p), B_AB=x_path(step_B_AB.p), AB=x_path(step_AB.p))
    # symmetric relative distance 2|x - y| / |x + y| (dimensionless, one per parameter)
    rel(x, y) = maximum(2 .* abs.(x .- y) ./ abs.(x .+ y))
    Δx_AB(path) = rel(xs[path], xs.AB)  # distance to the joint optimum
    paths = (; A_AB=path_stats(step_A, step_A_AB),
        B_AB=path_stats(step_B, step_B_AB), AB=path_stats(step_AB))
    Δmax = max(Δx_AB(:A_AB), Δx_AB(:B_AB), rel(xs.A_AB, xs.B_AB))
    same_optimum = Δmax <= 10TOL

    # parameters sitting on one of their box bounds: the optimum is not interior
    function on_bound(p)
        values = collect(p[SYMS_AB])
        return [s for (s, v, l, u) in zip(SYMS_AB, values,
            vcat(lower_A, lower_B), vcat(upper_A, upper_B))
                      if isapprox(v, l; rtol=1e-3) || isapprox(v, u; rtol=1e-3)]
    end
    bound_hits = (;
        A_AB=on_bound(step_A_AB.p), B_AB=on_bound(step_B_AB.p), AB=on_bound(step_AB.p))

    # Report: total power, wall-clock, ipopt iterations, relative distance to the joint optimum
    println("\nParameter paths (tol = $TOL, boxes ×3 on A, ×1.5 on B)")
    println(rpad("path", 14), rpad("power", 15), rpad("wall-clock", 12), rpad("solver", 10),
        lpad("iterations", 10), "  ", rpad("max 2|x-y|/(x+y)", 18), "bounds")
    for (name, key, path) in (("A → (A, B)", :A_AB, paths.A_AB),
        ("B → (A, B)", :B_AB, paths.B_AB), ("(A, B)", :AB, paths.AB))
        println(rpad(name, 14), rpad("$(round(path.P / 1000; digits=2)) kW", 15),
            rpad("$(round(path.t_wall; digits=1)) s", 12),
            rpad("$(round(path.t_solve; digits=1)) s", 10),
            lpad(path.iter, 10), "  ", rpad(round(Δx_AB(key); sigdigits=3), 18),
            isempty(bound_hits[key]) ? "interior" : "on bound: $(join(bound_hits[key], ", "))")
    end
    println("\nfinal parameters (A → (A, B) | B → (A, B) | (A, B))")
    for (i, sym) in enumerate(SYMS_AB)
        println(rpad(LABELS[sym], 14), join(round.((xs.A_AB[i], xs.B_AB[i], xs.AB[i]);
            sigdigits=8), "  "))
    end
    println(same_optimum ? "→ all paths reach the same optimum" :
            "→ the paths reach DIFFERENT optima")

    # Objective landscape around the joint optimum, one parameter at a time
    LANDSCAPE = landscape(step_AB.p, SYMS_AB; lower=lower_AB, upper=upper_AB)
end

## Plot: wind sweep (needs the PARAMETRIC section)
"""
Three panels of one wind sweep: average power, line tension (min–max envelope),
and the re-optimized parameters normalized by their maximum over the wind
range. The dotted black curve is the 10 m/s reference held fixed.
"""
function plot_wind_sweep(winds, sweep, i_ref)
    (; syms, no_opt, opt, theta) = sweep
    power_vec, min_tension_vec, max_tension_vec, _ = invert(opt)
    power0, min_tension0, max_tension0, _ = invert(no_opt)
    theta_norm = [maxnormalize(θ) for θ in invert(theta)]
    ref = "Fixed at $(round(Int, WIND_REF)) m/s"
    ref_style = (; ls=:dot, c=:black)

    # Average power
    f1 = plot(ylabel="Power (kW)", legend=:topleft)
    plot!(f1, winds, power0 ./ 1000; ref_style..., label=ref)
    plot!(f1, winds, power_vec ./ 1000; m=:circle, c=1, label="Re-optimized")

    # Line tension, min–max envelope (lower is better for wearing the lines out)
    f2 = plot(ylabel="Line tension (kN)", legend=:topleft)
    plot!(f2, winds, [max_tension0 ./ 1000, min_tension0 ./ 1000];
        ref_style..., label=[ref ""])
    plot!(f2, winds, max_tension_vec ./ 1000; fillrange=min_tension_vec ./ 1000,
        c=1, fillalpha=0.1, label="Re-optimized")
    plot!(f2, winds, min_tension_vec ./ 1000; c=1)

    # Optimized parameters, normalized by their maximum over the wind range
    n = length(syms)
    f3 = plot(yticks=((0, 1), ("0", "max")), ylim=(-0.1, 1.1), ylabel="Parameters",
        xformatter=:plain, xlabel="Wind speed (m/s)", legend=:bottomright)
    plot!(f3, winds, [[p for _ in winds] for p in invert(theta_norm)[i_ref]];
        lw=1.5, ref_style...)
    plot!(f3, winds, theta_norm; marker=[:rect :diamond :star6 :utriangle :dtriangle][1:n],
        ls=[:dash :dot :dashdot :dashdotdot :dash][1:n], c=permutedims(collect(1:n)),
        label=permutedims([LABELS[s] for s in syms]))

    # three panels of a third of the text height, full text width; the margins
    # leave room for the axis labels (the theme sets margin = 0)
    fig = plot(f1, f2, f3; layout=(3, 1), size=plot_size(1; height=3 * 0.375),
        leftmargin=8Plots.mm, bottommargin=4Plots.mm)
    return (; fig, panels=(f1, f2, f3))
end

if PARAMETRIC
    my_defaults()
    default(xticks=4:2:16, xlim=(4, 16), xformatter=x -> "")

    for (name, sweep) in pairs(sweeps)
        (; fig, panels) = plot_wind_sweep(winds, sweep, i_ref)
        display(fig)
        savefig(fig, joinpath(OUTPATH, "optimization_wind_sweep_$name.pdf"))

        for (f, panel) in zip(panels, ("power", "tension", "params"))
            plot!(f, size=SIZE_FULL, leftmargin=8Plots.mm, bottommargin=6Plots.mm)
            display(f)
            savefig(f, joinpath(OUTPATH, "optimization_$(panel)_$name.pdf"))
        end
    end

    #= Observation:
    Each wind speed is run twice: with the joint 10 m/s optimum held fixed
    (dotted black), and with one block re-optimized inside the same boxes.
    Re-optimizing the arm A = (ℓ, I, b) is what adapts the system to the wind:
    +19 % power at 5 m/s, nothing at the reference speed, +6 % at 15 m/s, and at
    15 m/s the line tension at the same time falls from 34.9 to 32.0 kN. The
    optimal arm shortens and lightens as the wind picks up (ℓ 40.7 → 32.9 m,
    I 1480 → 1030 kg m²) while the braking coefficient b more than doubles
    (1330 → 3050 kg m² s⁻¹ rad⁻¹): the cycle turns faster and the generator
    takes more out of the same lines.
    Re-optimizing the shape B = (Δθ, Δφ) changes nothing (at most +2 % at
    5 m/s). Its optimum sits on the lower box bound of both angles, so with the
    box fixed the solver has no room to flatten the lemniscate further, and
    (A, B) reproduces the A curve. The shape is already as flat as the box
    allows, so the wind-speed adaptation is carried entirely by the arm.
    The A gain vanishes at the reference speed by construction — 10 m/s is where
    the fixed parameters were optimized — and grows away from it on both sides.
    =#
end

## Plot: parameter paths, numerical results (left) and performance (right)
NODE_XY = (; start=(0.5, 0.95), A=(0.85, 0.5), B=(0.15, 0.5), AB=(0.5, 0.05))
PATH_FONTSIZE = 9  # annotation size for the diagram, smaller than the paper text

"Node box (data coordinates) around (x, y), tall enough for `nlines` text lines."
function node_box(x, y, nlines; w=0.48)
    h = 0.10 + 0.055nlines
    return Shape([x - w/2, x + w/2, x + w/2, x - w/2],
        [y - h/2, y - h/2, y + h/2, y + h/2])
end

"Block of text centered on (x, y), one `annotate!` per line."
function annotate_block!(plt, x, y, txt; color=:black, dy=0.06, fontsize=PATH_FONTSIZE)
    lines = split(txt, '\n')
    for (k, line) in enumerate(lines)
        annotate!(plt, x, y + (length(lines)/2 - k) * dy, (line, fontsize, :center, color))
    end
    return plt
end

"""
Optimization-path diagram: start at the top, A right, B left, both converging to
(A, B) at the bottom, plus the joint (A, B) optimization as the vertical arrow.
`nodes` maps a node name (:start, :A, :B, :AB) to its text; `edges` is a vector of
`(; from, to, label, c, at)` where `at` is the data coordinate of the label (an
empty label draws no text). Arrows stop on the border of the node boxes.
"""
function plot_paths(nodes, edges; title="", accent=:black)
    nlines = Dict(name => length(split(txt, '\n')) for (name, txt) in pairs(nodes))
    plt = plot(; xlims=(-0.15, 1.15), ylims=(-0.18, 1.15), aspect_ratio=:equal,
        framestyle=:none, grid=false, ticks=nothing, legend=false, title,
        titlefontsize=TICKFONTSIZE)

    for name in keys(nodes)
        (x, y) = NODE_XY[name]
        plot!(plt, node_box(x, y, nlines[name]); c=:white,
            lc=name === :AB ? accent : :black, lw=1.5)
    end
    for e in edges
        (x1, y1), (x2, y2) = NODE_XY[e.from], NODE_XY[e.to]
        dx, dy = x2 - x1, y2 - y1
        # stop the arrow on the border of each box
        trim(name) = min(0.24 / abs(dx), (0.05 + 0.0275nlines[name]) / abs(dy))
        t1, t2 = trim(e.from), trim(e.to)
        ax, ay, bx, by = x1 + dx * t1, y1 + dy * t1, x2 - dx * t2, y2 - dy * t2
        plot!(plt, [ax, bx], [ay, by]; arrow=:closed, c=e.c, lw=1.5)
        isempty(e.label) || annotate_block!(plt, e.at[1], e.at[2], e.label; color=e.c)
    end
    for (name, txt) in pairs(nodes)
        (x, y) = NODE_XY[name]
        annotate_block!(plt, x, y, txt)
    end
    return plt
end

if SHAPE
    my_defaults()
    fmt(x) = replace(string(round(x; sigdigits=2)), r"\.0$" => "")  # 2500, 8.0e-6
    secs(t) = t < 10 ? string(round(t; digits=1)) : string(round(Int, t))  # keep the fast steps visible
    Δcolor(Δ) = Δ <= 10TOL ? :green : :red

    # Left panel: power at each node, distance of each path to the joint optimum
    nodes_num = (
        start="start\nP = $(round(P_start / 1000; digits=1)) kW",
        A="A\nP = $(round(step_A.P / 1000; digits=1)) kW",
        B="B\nP = $(round(step_B.P / 1000; digits=1)) kW",
        AB="(A, B)\nP = $(round(step_AB.P / 1000; digits=1)) kW\nΔmax = $(fmt(Δmax))" *
            (isempty(bound_hits.AB) ? "" : "\non bound: $(join(bound_hits.AB, ", "))"),
    )
    edges_num = [
        (; from=:start, to=:A, label="opt A", c=:black, at=(0.80, 0.80)),
        (; from=:start, to=:B, label="opt B", c=:black, at=(0.20, 0.80)),
        (; from=:start, to=:AB, label="", c=:black, at=(0.5, 0.5)),
        (; from=:A, to=:AB, label="opt (A, B)\nΔ = $(fmt(Δx_AB(:A_AB)))",
            c=Δcolor(Δx_AB(:A_AB)), at=(0.90, 0.28)),
        (; from=:B, to=:AB, label="opt (A, B)\nΔ = $(fmt(Δx_AB(:B_AB)))",
            c=Δcolor(Δx_AB(:B_AB)), at=(0.10, 0.28)),
    ]
    panel_num = plot_paths(nodes_num, edges_num; title="Numerical",
        accent=same_optimum ? :green : :red)

    # Right panel: same layout, wall-clock time and ipopt iterations per step,
    # and the totals of the three paths in the bottom node
    nodes_perf = (
        start="start",
        A="A",
        B="B",
        AB="(A, B)\ndirect $(secs(paths.AB.t_wall)) s/$(paths.AB.iter) it\n" *
            "A→(A, B) $(secs(paths.A_AB.t_wall)) s/$(paths.A_AB.iter) it\n" *
            "B→(A, B) $(secs(paths.B_AB.t_wall)) s/$(paths.B_AB.iter) it",
    )
    edges_perf = [
        (; from=:start, to=:A, label="opt A\n$(secs(step_A.t_wall)) s/$(step_A.iter) it",
            c=:black, at=(0.80, 0.80)),
        (; from=:start, to=:B, label="opt B\n$(secs(step_B.t_wall)) s/$(step_B.iter) it",
            c=:black, at=(0.20, 0.80)),
        (; from=:start, to=:AB, label="", c=:black, at=(0.5, 0.5)),
        (; from=:A, to=:AB, label="opt (A, B)\n$(secs(step_A_AB.t_wall)) s/$(step_A_AB.iter) it",
            c=:black, at=(0.90, 0.28)),
        (; from=:B, to=:AB, label="opt (A, B)\n$(secs(step_B_AB.t_wall)) s/$(step_B_AB.iter) it",
            c=:black, at=(0.10, 0.28)),
    ]
    panel_perf = plot_paths(nodes_perf, edges_perf; title="Performance")

    fig_paths = plot(panel_num, panel_perf; layout=(1, 2), size=plot_size(1; height=0.5))
    display(fig_paths)
    savefig(fig_paths, joinpath(OUTPATH, "optimization_paths.pdf"))

    #= Observation:
    All three paths start from the same unoptimized parameters (P = 12.8 kW) and
    end with the same joint (A, B) solve; they differ only by the block optimized
    before it:
      - A → (A, B): 16.1 kW, 20 ipopt iterations
      - B → (A, B): 16.1 kW, 17 iterations
      - (A, B) alone: 16.1 kW, 10 iterations
    They reach the same optimum to solver accuracy: the largest symmetric
    relative difference 2|x - y| / (x + y) between the three final points is
    7.4e-6 (10 tol = 1e-4), i.e. the pre-optimization does not move the joint
    optimum. Blocking the arm first, then the shape, gave a different result
    only because that second step was a *restricted* problem (shape optimized at
    fixed arm); as soon as both blocks are free, the extra start is forgotten.
    Cost: the direct joint solve is the cheapest (10 iterations, 0.6 s); the
    pre-optimization adds a full block solve for nothing (1.2 s and 0.8 s), so
    warm-starting with A or B is not worth it here.
    All three final points sit on both B box bounds (Δθ = 2.5°, Δφ = 33°): the
    shape wants a flatter lemniscate than the ±50 % box allows, so the shape
    optimum is set by the box rather than by the objective.
    =#
end

## Plot: objective landscape around the joint optimum (needs the SHAPE section)
"""
Plot a normalized landscape from `landscape`: one pane per parameter, two
columns (the last cell stays empty for an odd count). The x axis is the
parameter in its own units, spanning ±`window` around the optimum, whose value
is the middle tick (`-10%`, value, `+10%`). The cost function (still
normalized by its optimum) is a line, solid inside the box constraint and
dotted outside it, the crossing interpolated linearly between the two samples
that bracket the bound. Each box constraint is a black double vertical line,
always drawn (clipped when it falls outside the window) and on top of the line.
"""
function plot_landscape(land)
    (; xs, syms, ys, window, box, center) = land
    nrows = cld(length(syms), 2)
    lo, hi = 1 - window, 1 + window

    # the margins are per panel (the theme sets margin = 0) and leave room for
    # the labels the two-column grid would otherwise clip
    panels = [
        begin
            x0 = center[i]  # the optimum, in physical units
            ticks = ([x0 * lo, x0, x0 * hi],
                [f"-{100window:.0f}%", f"{x0:.2f}", f"+{100window:.0f}%"])
            plt = plot(; xlims=(x0 * lo, x0 * hi), widen=true, xticks=ticks,
                title="$(LABELS[sym])\n($(UNITS[sym]))",
                ylabel=(isodd(i) ? "P / P*" : ""),
                leftmargin=(isodd(i) ? 9 : 5)Plots.mm, rightmargin=5Plots.mm,
                bottommargin=7Plots.mm, topmargin=4Plots.mm)
            # value of the sampled curve at a factor x, linearly interpolated
            interp(x) = begin
                k = clamp(searchsortedlast(xs, x), 1, length(xs) - 1)
                ys[i][k] + (x - xs[k]) / (xs[k+1] - xs[k]) * (ys[i][k+1] - ys[i][k])
            end
            # factors from a to b, with the ends inserted (they need not be samples)
            piece(a, b) = vcat(a, xs[(a .< xs) .& (xs .< b)], b)
            xa, xb = max(box[i][1], lo), min(box[i][2], hi)  # box ∩ window
            xa < xb && plot!(plt, x0 .* piece(xa, xb), interp.(piece(xa, xb));
                c=1, ls=:solid)
            for (a, b) in ((lo, min(xa, hi)), (max(xb, lo), hi))
                a < b && plot!(plt, x0 .* piece(a, b), interp.(piece(a, b)); c=1, ls=:dot)
            end
            for edge in box[i], k in (-1, 1)  # box constraint: double line, on top
                vline!(plt, [x0 * (edge + k * window / 40)]; c=:black, lw=1.2)
            end
            plt
        end for (i, sym) in enumerate(syms)
    ]

    return plot(panels...; layout=(nrows, 2), size=plot_size(1; height=nrows * 0.42))
end

if SHAPE
    my_defaults()
    fig_land = plot_landscape(LANDSCAPE)
    display(fig_land)
    savefig(fig_land, joinpath(OUTPATH, "optimization_landscape.pdf"))

    #= Observation:
    Each pane varies one parameter around the joint optimum and re-solves the
    limit cycle, so the curve is the objective as a function of that parameter
    alone. The arm length r is the stiffest direction (≈ 5 % of power over
    ±10 %), the inertia and the braking coefficient are about four and ten times
    flatter. The shape directions are asymmetric: Δθ and Δφ sit on their lower
    box bound, the objective still rises just below it, and Δφ already loses
    4 % at +10 % — the bound, not a flat objective, stops the shape. Each x axis
    spans ±10 % of one parameter in its own units, and the dashed–solid split
    marks the box bound of that parameter.
    =#
end
