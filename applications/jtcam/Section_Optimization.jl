# Optimization section.
#
# Two independent studies, each guarded by a flag at the top of the file:
#  - PARAMETRIC: for each wind speed, maximize the average power over
#    A = (r, I_eq, torque_slope), with the other parameters at their 10 m/s
#    optimum. Reproduction of
#    test/publications/ECC2026/6_fixed_parametrized_optimization.jl.
#  - SHAPE: do different optimization paths over A (arm) and B = (Δθ, Δφ)
#    (lemniscate half-height and half-width) reach the same optimum?  A → B,
#    B → A and both at once (A, B), started from the same parameters; runtime
#    and ipopt iterations are reported for each path.
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
# value of one parameter.
function compute_results(param, value; syms, p0, sense, tol)
    p0 = copy(p0)
    p0[param] = value

    # Compute limit cycle and corresponding stats
    vbp0 = PMP.build_vbpara(p0)
    lc = LC.compute_limit_cycle(vbp0; sense=sense, save_everystep=true, tol=tol)
    stats_no_opt = trajectory_stats(lc)

    # Optimize the selected parameters
    lower, upper = OPT.make_bounds(p0, syms)
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

function range_optim(param, values; syms, p_ref, sense, tol)
    return invert([
        compute_results(param, v; syms, p0=p_ref, sense, tol)
        for v in values
    ])
end

maxnormalize(x) = x ./ maximum(x)

"""
One optimization step: optimize `syms` from `p` and record the new parameters,
the optimal power and the cost of the solve (`stats.elapsed_time` is ipopt's own
timer, `t_wall` the wall-clock of the whole call including its initial guess).
"""
function opt_step(p, syms, lower, upper)
    solve(; compile_only=false) = with_logger(NullLogger()) do
        OPT.optimize(p, syms, lower, upper; sense=SENSE, tol=TOL, max_wall_time=ifelse(compile_only, 1e-1, MAX_WALL_TIME))
    end
    # The first solve with a given variable layout pays Julia's compilation of the
    # dual-number Jacobian inside ipopt's wall clock (100-250 s here), which eats
    # the time cap and truncates the solve. Warm the layout up, keep the second run.
    solve(compile_only=true)
    t_wall = time()
    solution, stats, model = solve()
    stats.status == :first_order || @warn "opt_step($(join(syms, ", "))) did not converge" stats.status stats.iter
    p_next = copy(p)
    p_next[syms] .= solution.params[syms]  # physical units

    return (; p=p_next, P=solution.P, iter=stats.iter, t_solve=stats.elapsed_time,
        t_wall=time() - t_wall, residual=norm(model.c!(model.meta.ucon, stats.solution)))
end

"Total power, wall-clock time and ipopt iterations over consecutive optimization steps."
function path_stats(steps...)
    return (; P=last(steps).P, t_wall=sum(s -> s.t_wall, steps),
        t_solve=sum(s -> s.t_solve, steps), iter=sum(s -> s.iter, steps))
end

## Numerical: wind sweep
if PARAMETRIC
    winds = range(5, 15, length=11)
    i_ref = findfirst(x -> x >= WIND_REF, winds)
    p0 = PMP.build_para(v_ref=WIND_REF)  # Reference parameters at 10 m/s

    # Optimize once at the reference wind speed, then keep those parameters.
    # Warm up the same variable layout first: it is the one whose compilation
    # would otherwise be charged to this solve (see opt_step).
    opt_step(p0, SYMS_A, OPT.make_bounds(p0, SYMS_A)...)
    _, _, theta0 = compute_results(:v_ref, WIND_REF; syms=SYMS_A, p0, sense=SENSE, tol=TOL)
    p_ref = copy(p0)
    p_ref[SYMS_A] .= theta0

    # Sweep the wind, with and without re-optimizing the parameters
    stats_no_opt_vec, stats_opt_vec, theta_vec =
        range_optim(:v_ref, winds; syms=SYMS_A, p_ref, sense=SENSE, tol=TOL)
end

## Numerical: parameter paths A = (r, I_eq, torque_slope), B = (Δθ, Δφ)
if SHAPE
    p_start = PMP.build_para(v_ref=WIND_REF)
    lower_A, upper_A = OPT.make_bounds(p_start, SYMS_A)
    lower_B, upper_B = OPT.make_bounds(p_start, SYMS_B; multiplier=1.5)  # ×3 on Δφ leaves the window
    lower_AB, upper_AB = vcat(lower_A, lower_B), vcat(upper_A, upper_B)

    step_A = opt_step(p_start, SYMS_A, lower_A, upper_A)
    step_B = opt_step(p_start, SYMS_B, lower_B, upper_B)
    step_A_B = opt_step(step_A.p, SYMS_B, lower_B, upper_B)  # A, then B
    step_B_A = opt_step(step_B.p, SYMS_A, lower_A, upper_A)  # B, then A
    step_AB = opt_step(p_start, vcat(SYMS_A, SYMS_B), lower_AB, upper_AB)  # joint

    P_start = trajectory_stats(LC.compute_limit_cycle(PMP.build_vbpara(p_start);
        sense=SENSE, tol=TOL, save_everystep=true))[1]

    # Path totals and final parameter vectors
    x_path(p) = vcat(collect(p[SYMS_A]), collect(p[SYMS_B]))
    xs = (; A_B=x_path(step_A_B.p), B_A=x_path(step_B_A.p), AB=x_path(step_AB.p))
    # symmetric relative distance 2|x - y| / |x + y| (dimensionless, one per parameter)
    rel(x, y) = maximum(2 .* abs.(x .- y) ./ abs.(x .+ y))
    Δx_AB(path) = rel(xs[path], xs.AB)  # distance to the joint optimum
    paths = (; A_B=path_stats(step_A, step_A_B),
        B_A=path_stats(step_B, step_B_A), AB=path_stats(step_AB))
    Δmax = max(Δx_AB(:A_B), Δx_AB(:B_A), rel(xs.A_B, xs.B_A))
    same_optimum = Δmax <= 10TOL

    # parameters sitting on one of their box bounds: the optimum is not interior
    function on_bound(p)
        values = collect(p[vcat(SYMS_A, SYMS_B)])
        return [s for (s, v, l, u) in zip(vcat(SYMS_A, SYMS_B), values,
            vcat(lower_A, lower_B), vcat(upper_A, upper_B))
                      if isapprox(v, l; rtol=1e-3) || isapprox(v, u; rtol=1e-3)]
    end
    bound_hits = (; A_B=on_bound(step_A_B.p), B_A=on_bound(step_B_A.p), AB=on_bound(step_AB.p))

    # Report: total power, wall-clock, ipopt iterations, relative distance to the joint optimum
    println("\nParameter paths (tol = $TOL, boxes ×3 on A, ×1.5 on B)")
    println(rpad("path", 10), rpad("power", 15), rpad("wall-clock", 12), rpad("solver", 10),
        lpad("iterations", 10), "  ", rpad("max 2|x-y|/(x+y)", 18), "bounds")
    for (name, key, path) in (("A → B", :A_B, paths.A_B), ("B → A", :B_A, paths.B_A),
        ("(A, B)", :AB, paths.AB))
        println(rpad(name, 10), rpad("$(round(path.P / 1000; digits=2)) kW", 15),
            rpad("$(round(path.t_wall; digits=1)) s", 12),
            rpad("$(round(path.t_solve; digits=1)) s", 10),
            lpad(path.iter, 10), "  ", rpad(round(Δx_AB(key); sigdigits=3), 18),
            isempty(bound_hits[key]) ? "interior" : "on bound: $(join(bound_hits[key], ", "))")
    end
    println("\nfinal parameters (A → B | B → A | joint)")
    for (i, sym) in enumerate(vcat(SYMS_A, SYMS_B))
        println(rpad(string(sym), 14), join(round.((xs.A_B[i], xs.B_A[i], xs.AB[i]);
            sigdigits=8), "  "))
    end
    println(same_optimum ? "→ all paths reach the same optimum" :
            "→ the paths reach DIFFERENT optima")
end

## Plot: wind sweep (needs the PARAMETRIC section)
if PARAMETRIC
    my_defaults()
    power_vec, min_tension_vec, max_tension_vec, _ = invert(stats_opt_vec)
    power_no_opt_vec, min_tension_no_opt_vec, max_tension_no_opt_vec, _ =
        invert(stats_no_opt_vec)
    theta_norm = [maxnormalize(theta) for theta in invert(theta_vec)]

    ref_style = (; ls=:dot, c=:black)
    default(xticks=4:2:16, xlim=(4, 16), xformatter=x -> "")

    # Average power
    f1 = plot(ylabel="Power (kW)", legend=:topleft)
    plot!(f1, winds, power_no_opt_vec ./ 1000; ref_style..., label="Fixed at $(round(Int, WIND_REF)) m/s")
    plot!(f1, winds, power_vec ./ 1000; m=:circle, c=1, label="Re-optimized")

    # Line tension, min–max envelope (lower is better for wearing the lines out)
    f2 = plot(ylabel="Line tension (kN)", legend=:topleft)
    plot!(f2, winds, [max_tension_no_opt_vec ./ 1000, min_tension_no_opt_vec ./ 1000];
        ref_style..., label=["Fixed at $(round(Int, WIND_REF)) m/s" ""])
    plot!(f2, winds, max_tension_vec ./ 1000; fillrange=min_tension_vec ./ 1000,
        c=1, fillalpha=0.1, label="Re-optimized")
    plot!(f2, winds, min_tension_vec ./ 1000; c=1)

    # Optimized parameters, normalized by their maximum over the wind range
    f3 = plot(yticks=((0, 1), ("0", "max")), ylim=(-0.1, 1.1), ylabel="Parameters",
        xformatter=:plain, xlabel="Wind speed (m/s)", legend=:bottomright)
    plot!(f3, winds, [[p for _ in winds] for p in invert(theta_norm)[i_ref]]; lw=1.5, ref_style...)
    plot!(f3, winds, theta_norm; marker=[:rect :diamond :star6], ls=[:dash :dot :dashdot], c=[1 2 3],
        label=["Line length ℓ" "Arm inertia I" "Braking coeff. b"])

    # three panels of a third of the text height, full text width; the margins
    # leave room for the axis labels (the theme sets margin = 0)
    fig = plot(f1, f2, f3; layout=(3, 1), size=plot_size(1; height=3 * 0.375),
        leftmargin=8Plots.mm, bottommargin=4Plots.mm)
    display(fig)

    for (f, name) in zip((f1, f2, f3), ("power", "tension", "params"))
        plot!(f, size=SIZE_FULL, leftmargin=8Plots.mm, bottommargin=6Plots.mm)
        display(f)
        savefig(f, joinpath(OUTPATH, "optimization_$name.pdf"))
    end
    savefig(fig, joinpath(OUTPATH, "optimization_wind_sweep.pdf"))

    #= Observation:
    For each wind speed we maximize the average power (top) over the parameters
    (bottom), while the line tension (min–max envelope, middle) is kept low to
    wear the lines less. For comparison, power and tension are also computed for
    the parameters optimized at 10 m/s and held fixed (black, dotted).
    The 10 m/s parameters are strongly resilient to a change in wind speed:
    re-optimizing brings close to no improvement in power and tension, but the
    optimized parameters produce more power while putting less strain on the
    lines. Along this parameterized optimization the system behaves very
    linearly, although the arm length is kept constant (and other parameters),
    so that a perfect scaling of the model with wind speed is impossible.
    Approximate correlations: power ∝ v³, tension ∝ v², line length ℓ ∝ 1/v,
    arm inertia I ∝ 1/v, arm braking coefficient b ∝ v, period ∝ 1/v.
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
        (; from=:A, to=:AB, label="opt B\nΔ = $(fmt(Δx_AB(:A_B)))",
            c=Δcolor(Δx_AB(:A_B)), at=(0.90, 0.28)),
        (; from=:B, to=:AB, label="opt A\nΔ = $(fmt(Δx_AB(:B_A)))",
            c=Δcolor(Δx_AB(:B_A)), at=(0.10, 0.28)),
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
            "A→B $(secs(paths.A_B.t_wall)) s/$(paths.A_B.iter) it\n" *
            "B→A $(secs(paths.B_A.t_wall)) s/$(paths.B_A.iter) it",
    )
    edges_perf = [
        (; from=:start, to=:A, label="opt A\n$(secs(step_A.t_wall)) s/$(step_A.iter) it",
            c=:black, at=(0.80, 0.80)),
        (; from=:start, to=:B, label="opt B\n$(secs(step_B.t_wall)) s/$(step_B.iter) it",
            c=:black, at=(0.20, 0.80)),
        (; from=:start, to=:AB, label="", c=:black, at=(0.5, 0.5)),
        (; from=:A, to=:AB, label="opt B\n$(secs(step_A_B.t_wall)) s/$(step_A_B.iter) it",
            c=:black, at=(0.90, 0.28)),
        (; from=:B, to=:AB, label="opt A\n$(secs(step_B_A.t_wall)) s/$(step_B_A.iter) it",
            c=:black, at=(0.10, 0.28)),
    ]
    panel_perf = plot_paths(nodes_perf, edges_perf; title="Performance")

    fig_paths = plot(panel_num, panel_perf; layout=(1, 2), size=plot_size(1; height=0.5))
    display(fig_paths)
    savefig(fig_paths, joinpath(OUTPATH, "optimization_paths.pdf"))

    #= Observation:
    The three paths start from the same unoptimized parameters (P = 12.8 kW) and
    do NOT reach the same optimum, so the square does not commute:
      - A → B: 14.3 kW, 17 ipopt iterations
      - B → A: 15.2 kW, 15 iterations
      - (A, B): 16.1 kW, 10 iterations
    Distances are the symmetric relative difference 2|x - y| / (x + y), which puts
    the five parameters on one scale (verdict threshold 10 tol = 1e-4): the two
    sequential paths end 50 % (A → B: 108 % of I_eq) away from the joint optimum,
    i.e. orders of magnitude above the tolerance.
    The arm (r, I, b) and the lemniscate shape (Δθ, Δφ) are strongly coupled:
    the arm optimum found on its own (14.0 kW) is a different local optimum from
    the arm values the joint solve settles on, so optimizing one block first
    locks the other into a different basin. Re-optimizing the arm after the shape
    helps more (+1.3 kW) than re-optimizing the shape after the arm (+0.3 kW).
    Cost: with the parameter boxes fixed, ipopt needs 7 to 17 iterations and
    about one second per solve, so no path is expensive; the iteration count is
    the robust cost metric.
    Caveat: the joint optimum sits on both B box bounds (Δθ = 2.5°, Δφ = 33°):
    the shape wants a flatter lemniscate than the ±50 % box allows, while the
    same box is not active for the sequential paths. The gap is therefore a lower
    bound on the path dependence, and it would grow with a wider shape box.
    =#
end
