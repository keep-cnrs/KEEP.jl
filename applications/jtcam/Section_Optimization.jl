# Optimization section.
#
# The optimization variables are Θ = (ℓ, I, b, Δφ): the tether length, the arm
# inertia, the braking coefficient and the lemniscate half-width. The
# lemniscate half-height Δθ is held at its nominal value; the last study here
# checks what that choice costs (power, line tension, kite acceleration) and a
# 10×10 map of the same quantities over the shape plane (Δθ, Δφ).
#
# The reference is the default parameters with the torque slope set by the
# kinetic criterion b = I ω (`PMP.optimal_torque_slope`, one limit cycle per
# iteration), not the default torque slope. From there:
#  - OPTIMIZATION: maximize the average power over Θ at 9 m/s, then repeat the
#    optimization at each wind speed and compare with the 9 m/s optimum held
#    fixed. ℓ spans 15–200 m and b 1–1e6; the Δφ box is (10°, 90°), and the
#    inertia floor is set by sizing the aluminium arm (see `I_BOUNDS`).
#  - Δθ: sweep it ±10 % around its nominal value at the optimum with Θ fixed,
#    and plot the power, the line-tension envelope and the maximum kite
#    acceleration.
#  - SHAPE: 10×10 map of the same quantities over (Δθ, Δφ) at the optimal
#    (ℓ, I, b), again without re-optimizing anything.
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
include(joinpath(@__DIR__, "time_series.jl"))

## Settings
const TOL = 1e-5
const MAX_WALL_TIME = 120.0
const SWEEP_WALL_TIME = 30.0  # a converged wind-sweep solve takes about a second
const WIND_REF = 9.0
const SENSE = +
const SYMS = [:l, :I, :b, :Δφ]  # ℓ, I, b and the half-width
const L_BOUNDS = (15.0, 200.0)               # tether length ℓ (m)
const B_BOUNDS = (1.0, 1.0e6)                # braking slope b
# Assumptions for I lower bound :
# Maximum line tension T ≈ 34 kN at w0 = 15 m/s
# Moment is then M = T a = 68 kN m. For 6061-T6 (σy ≈ 276 MPa) with a safety factor 2,
# W = M/σallow ≈ 491 cm³, met by a Ø250 mm tube with an 11.5 mm wall (47 kg).
# Its inertia about the pivot is I = m p.a²/3 ≈ 62 kg m² (solid bar: 165 kg m²).
# No upper bound: extra mass is allowed if it buys power.
const I_BOUNDS = (62.0, 1.0e6)               # arm inertia I (kg m²)
const DPHI_BOUNDS = deg2rad.((10.0, 90.0))   # fixed shape box, in radians
const DTHETA_WINDOW = 0.1                     # ±10 % around the nominal Δθ
const SHAPE_WINDOW = 0.2                      # ±20 % for the shape map
const GRID_N = 10

# Parameter names and units for the axis titles and the legends
const LABELS = Dict(:l => "ℓ", :I => "I", :b => "b",
    :Δθ => "Δθ", :Δφ => "Δφ")
const UNITS = Dict(:l => "m", :I => "kg m²", :b => "kg m² s⁻¹ rad⁻¹",
    :Δθ => "rad", :Δφ => "rad")

## Helpers
"Average power, line-tension extrema, peak kite acceleration and period of a cycle."
function cycle_stats(sol; n=2001)
    ts = range(sol.t[1], sol.t[end], n)
    us = sol.(ts)
    p = sol.prob.p
    tension = PM4.compute_line_tension.(us, Ref(p))
    ok = isfinite.(tension)
    acceleration = norm.(PM4.kite_acceleration.(us, Ref(p)))
    return (; P=sol.u[end][5] / sol.t[end],
        Tmin=minimum(tension[ok]), Tmax=maximum(tension[ok]),
        amax=maximum(acceleration), period=sol.t[end])
end

# Limit cycle, optimization, and limit cycle of the optimized system, for one
# value of one parameter. `lower`/`upper` are the box constraints of `syms`.
# `p_fixed` is the parameter set the no-optimization cycle is taken at, so a
# continuation can start from a neighbour's optimum `p0` without moving the
# reference curve.
function compute_results(param, value; syms, p0, p_fixed=p0, lower, upper, sense, tol,
    wall=MAX_WALL_TIME)
    p_fixed = copy(p_fixed)
    p_fixed[param] = value

    # Limit cycle at the reference parameters, before re-optimizing
    vbp0 = PMP.build_vbpara(p_fixed)
    lc0 = LC.compute_limit_cycle(vbp0; sense=sense, save_everystep=true, tol=tol)
    stats_no_opt = cycle_stats(lc0)

    # Optimize the selected parameters, starting from `p0`
    p0 = copy(p0)
    p0[param] = value
    solution, stats, _ = with_logger(NullLogger()) do
        OPT.optimize(p0, syms, lower, upper; sense, tol, max_wall_time=wall)
    end
    p_opt = copy(p0)
    p_opt[syms] .= solution.params[syms]  # physical units

    # Limit cycle at the optimized parameters
    vbp = PMP.build_vbpara(p_opt)
    lc = LC.compute_limit_cycle(vbp; sense=sense, save_everystep=true, tol=tol)
    stats_opt = cycle_stats(lc)

    return stats_no_opt, stats_opt, collect(p_opt[syms]), stats
end

"""
Wind sweep: hold `p_ref` fixed, or re-optimize `syms` inside `lower`/`upper` at
each wind speed, marching outward from the reference wind `i_ref`. Each point is
first solved from the 9 m/s optimum; if that does not converge, it is retried
from the neighbouring wind's optimum. The wide boxes make the objective
multimodal, so neither start alone succeeds at every wind speed.
Returns the statistics with and without re-optimization and the optimized
parameters.
"""
function wind_sweep(syms, p_ref, winds, lower, upper, i_ref)
    n = length(winds)
    no_opt = Vector{Any}(undef, n)
    opt = Vector{Any}(undef, n)
    theta = Vector{Vector{Float64}}(undef, n)
    ok = falses(n)
    from_theta(θ) = (p=copy(p_ref); p[syms].=θ; p)
    function run!(i, p0)
        r = compute_results(:w0, winds[i];
            syms, p0, p_fixed=p_ref, lower, upper, sense=SENSE, tol=TOL,
            wall=SWEEP_WALL_TIME)
        no_opt[i], opt[i], theta[i] = r[1], r[2], r[3]
        ok[i] = r[4].status == :first_order
    end
    for i in i_ref:-1:1
        run!(i, p_ref)
        !ok[i] && i < i_ref && run!(i, from_theta(theta[i+1]))
    end
    for i in (i_ref+1):n
        run!(i, p_ref)
        !ok[i] && run!(i, from_theta(theta[i-1]))
    end
    for i in eachindex(winds)
        ok[i] || @warn "re-optimization at $(winds[i]) m/s did not converge"
    end
    return (; syms, no_opt, opt, theta)
end

maxnormalize(x) = x ./ maximum(x)

"""
One optimization step with the compilation warm-up: the first solve with a
given variable layout pays Julia's compilation of the dual-number Jacobian
inside ipopt's wall clock, which can eat the time cap and truncate the solve.
Warm the layout up and keep the second run.
"""
function opt_step(p, syms, lower, upper)
    solve() = with_logger(NullLogger()) do
        OPT.optimize(p, syms, lower, upper; sense=SENSE, tol=TOL, max_wall_time=MAX_WALL_TIME)
    end
    solve()
    t_wall = time()
    solution, stats, _ = solve()
    stats.status == :first_order || @warn "opt_step($(join(syms, ", "))) did not converge" stats.status stats.iter
    p_next = copy(p)
    p_next[syms] .= solution.params[syms]  # physical units
    return (; p=p_next, P=solution.P, iter=stats.iter, t_wall=time() - t_wall)
end

"""
Δθ sweep at the optimum: vary Δθ by `factors` around its nominal value with Θ
fixed and take the limit cycle of the perturbed system. Starting from the same
`u0 = (α, τ, dα, dτ) = (0, τ0, 0, 1)` picks the same (dτ > 0) pumping cycle at
every point.
"""
function dtheta_study(p, factors; n=2001)
    u0 = PM4.build_u(; dτ=1)
    return map(factors) do f
        q = copy(p)
        q.Δθ = f * p.Δθ
        sol = LC.compute_limit_cycle(u0, PMP.build_vbpara(q);
            tol=TOL, save_everystep=true, tf=100.0)
        (; f, cycle_stats(sol; n)...)
    end
end

"""
10×10 shape map at fixed (ℓ, I, b): one limit cycle per (Δθ, Δφ), no
re-optimization. Rows are Δφ, columns are Δθ factors. A cell stays NaN when no
dτ > 0 pumping cycle is found.
"""
function shape_grid(p, factors, dphi_deg; n=401)
    u0 = PM4.build_u(; dτ=1)
    P = fill(NaN, length(dphi_deg), length(factors))
    Tmin = fill(NaN, size(P))
    Tmax = fill(NaN, size(P))
    amax = fill(NaN, size(P))
    for (j, dphi) in enumerate(dphi_deg), (i, f) in enumerate(factors)
        q = copy(p)
        q.Δθ = f * p.Δθ
        q.Δφ = deg2rad(dphi)
        vbp = PMP.build_vbpara(q)
        sol = try
            LC.compute_limit_cycle(u0, vbp; tol=TOL, save_everystep=true, tf=100.0)
        catch
            nothing
        end
        sol === nothing && continue
        # reject a trajectory that did not close
        norm(LC.endpoint_residuals(sol; sense=SENSE)) < 1e-3 || continue
        s = cycle_stats(sol; n)
        P[j, i], Tmin[j, i], Tmax[j, i], amax[j, i] =
            s.P / 1000, s.Tmin / 1000, s.Tmax / 1000, s.amax
    end
    return (; P, Tmin, Tmax, amax, factors, dphi_deg)
end

"""
Objective landscape around the optimal point `p`, one parameter at a time:
scale one parameter by each factor `x ∈ [1 - window, 1 + window]`, keep the
other parameters at their optimum, and take the power of the limit cycle at
those fixed parameters — the objective the optimizer maximizes once the cycle
is solved. Powers are normalized by their value at `x = 1`, so the optimum
plots at `(1, 1)`. `lower` and `upper` are the box constraints of `syms`
(physical units), returned as normalized factors so the caller can mark them.
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

## Numerical: 9 m/s reference, one optimization over Θ
p_start = PMP.build_para(w0=WIND_REF)
p_start.b = PMP.optimal_torque_slope(p_start)  # b = I ω
lower = [first(L_BOUNDS), first(I_BOUNDS), first(B_BOUNDS), first(DPHI_BOUNDS)]
upper = [last(L_BOUNDS), last(I_BOUNDS), last(B_BOUNDS), last(DPHI_BOUNDS)]

opt_ref = opt_step(p_start, SYMS, lower, upper)
p_ref = opt_ref.p  # the Θ optimum at the reference wind speed

## Numerical: wind sweep over Θ, from the 9 m/s optimum
winds = range(5, 15, length=11)
i_ref = findfirst(x -> x >= WIND_REF, winds)
sweep = wind_sweep(SYMS, p_ref, winds, lower, upper, i_ref)

## Numerical: Δθ sensitivity at the optimum
dtheta = dtheta_study(p_ref, range(1 - DTHETA_WINDOW, 1 + DTHETA_WINDOW, length=11))

## Numerical: 10×10 shape map at the optimal ℓ, I, b
shape = shape_grid(p_ref,
    range(1 - SHAPE_WINDOW, 1 + SHAPE_WINDOW, length=GRID_N),
    range(rad2deg.(DPHI_BOUNDS)..., length=GRID_N))  # shape_grid takes degrees

## Numerical: objective landscape around the Θ optimum
LANDSCAPE = landscape(p_ref, SYMS; lower, upper)

## Numerical: report
function print_final_table()
    stats = (init=cycle_stats(LC.compute_limit_cycle(PMP.build_vbpara(p_start);
        sense=SENSE, tol=TOL, save_everystep=true)),
        opt=cycle_stats(LC.compute_limit_cycle(PMP.build_vbpara(p_ref);
            sense=SENSE, tol=TOL, save_everystep=true)))
    cols = ((p_start, stats.init), (p_ref, stats.opt))
    fmt(x) = string(round(x; sigdigits=6))
    row(name, f) = println(rpad(name, 30), join([lpad(f(c), 14) for c in cols]))

    println("\nfinal parameters")
    println(rpad("", 30), join(lpad.(["init", "Θ"], 14)))
    for sym in SYMS
        row("$(LABELS[sym]) ($(UNITS[sym]))", c -> fmt(c[1][sym]))
    end
    row("Δθ (rad, fixed)", c -> fmt(c[1].Δθ))
    println()
    row("period (s)", c -> fmt(c[2].period))
    row("objective P (kW)", c -> fmt(c[2].P / 1000))
    row("min tension (kN)", c -> fmt(c[2].Tmin / 1000))
    row("max tension (kN)", c -> fmt(c[2].Tmax / 1000))
    row("max acceleration (m/s²)", c -> fmt(c[2].amax))
    println()

    println("wind sweep over Θ (reference at $(round(Int, WIND_REF)) m/s)")
    println(rpad("v (m/s)", 8), rpad("fixed P", 10), rpad("opt P", 10), rpad("gain", 8),
        rpad("Tmax opt", 10), rpad("ℓ (m)", 9), rpad("I", 9), rpad("b", 9), "Δφ (°)")
    for k in eachindex(winds)
        P0, P = sweep.no_opt[k].P / 1000, sweep.opt[k].P / 1000
        θ = sweep.theta[k]
        println(rpad(round(Int, winds[k]), 8), rpad(round(P0; digits=3), 10),
            rpad(round(P; digits=3), 10), rpad(f"{(P / P0 - 1) * 100:+.1f} %", 8),
            rpad(round(sweep.opt[k].Tmax / 1000; digits=2), 10),
            rpad(round(θ[1]; digits=2), 9), rpad(round(θ[2]; digits=1), 9),
            rpad(round(θ[3]; digits=1), 9), round(rad2deg(θ[4]); digits=2))
    end
    println()

    println("Δθ sensitivity at the Θ optimum (nominal Δθ = $(round(rad2deg(p_ref.Δθ); digits=1))°)")
    println(rpad("Δθ (%)", 8), rpad("P (kW)", 10), rpad("Tmin (kN)", 11),
        rpad("Tmax (kN)", 11), rpad("spread (kN)", 12), "max a (m/s²)")
    for d in dtheta
        println(rpad(f"{(d.f - 1) * 100:+.0f}", 8), rpad(round(d.P / 1000; digits=3), 10),
            rpad(round(d.Tmin / 1000; digits=3), 11), rpad(round(d.Tmax / 1000; digits=3), 11),
            rpad(round((d.Tmax - d.Tmin) / 1000; digits=3), 12), round(d.amax; digits=1))
    end
    println()

    println("shape map (no re-optimization, $(GRID_N)×$(GRID_N) limit cycles)")
    finite(Z) = filter(isfinite, Z)
    println("  valid: ", count(isfinite, shape.P), " of ", length(shape.P), " cells")
    for (name, Z, unit) in (("P", shape.P, "kW"), ("Tmax", shape.Tmax, "kN"),
        ("Tmin", shape.Tmin, "kN"))
        lo, hi = extrema(finite(Z))
        println("  ", rpad(name, 5), round(lo; digits=2), " … ", round(hi; digits=2), " ", unit)
    end
    k = argmax(ifelse.(isfinite.(shape.amax), shape.amax, -Inf))
    println("  max a: ", round(shape.amax[k]; digits=0), " m/s² at Δθ = ",
        round((shape.factors[k[2]] - 1) * 100; digits=0), " %, Δφ = ",
        round(shape.dphi_deg[k[1]]; digits=1), "°")
    println()

    println("landscape at the Θ optimum (power relative to the optimum)")
    println(rpad("parameter", 12), rpad("at -10 %", 12), rpad("at +10 %", 12), "box (÷/×)")
    for (i, sym) in enumerate(SYMS)
        println(rpad(LABELS[sym], 12), rpad(round(LANDSCAPE.ys[i][1]; digits=4), 12),
            rpad(round(LANDSCAPE.ys[i][end]; digits=4), 12),
            round.(LANDSCAPE.box[i]; digits=3))
    end
end

## Plot: wind sweep over Θ
"""
Three panels of the wind sweep: average power, line tension (min–max envelope),
and the re-optimized parameters normalized by their maximum over the wind
range. The dotted black curve is the 9 m/s reference held fixed, and the grey
dashed vertical line marks that reference wind speed in every panel.
"""
function plot_wind_sweep(winds, sweep, i_ref)
    (; syms, no_opt, opt, theta) = sweep
    power_vec = [s.P for s in opt]
    min_tension_vec = [s.Tmin for s in opt]
    max_tension_vec = [s.Tmax for s in opt]
    power0 = [s.P for s in no_opt]
    min_tension0 = [s.Tmin for s in no_opt]
    max_tension0 = [s.Tmax for s in no_opt]
    theta_norm = [maxnormalize(θ) for θ in invert(theta)]
    ref = "Fixed at $(round(Int, WIND_REF)) m/s"
    ref_style = (; ls=:dot, c=:black)

    # Average power
    f1 = plot(ylabel="Power (kW)", legend=:topleft)
    plot!(f1, winds, power0 ./ 1000; ref_style..., label=ref)
    plot!(f1, winds, power_vec ./ 1000; m=:circle, c=1, label="Re-optimized")

    # Line tension, min–max envelope (lower is better for wearing the lines out)
    f2 = plot(ylabel="Line tension (kN)", legend=:topleft)
    plot!(f2, winds, [max_tension0 ./ 1000, min_tension0 ./ 1000]; ref_style..., label=[ref ""])
    plot!(f2, winds, max_tension_vec ./ 1000; fillrange=min_tension_vec ./ 1000, c=1, fillalpha=0.1, label="Re-optimized")
    plot!(f2, winds, min_tension_vec ./ 1000; c=1)

    # Optimized parameters, normalized by their maximum over the wind range
    n = length(syms)
    f3 = plot(yticks=((0, 1), ("0", "max")), ylim=(-0.1, 1.1), ylabel="Parameters", xformatter=:plain, xlabel="Wind speed (m/s)", legend=:bottomright)
    # plot!(f3, winds, [[p for _ in winds] for p in invert(theta_norm)[i_ref]]; lw=1.5, ref_style...)
    plot!(f3, winds, theta_norm; c=permutedims(collect(1:n)), label=permutedims([LABELS[s] for s in syms]))

    # mark the wind speed the reference parameters were optimized at
    foreach(f -> vline!(f, [WIND_REF]; c=:gray50, ls=:dash, label=""), (f1, f2, f3))

    # three panels of a third of the text height, full text width; the margins
    # leave room for the axis labels (the theme sets margin = 0)
    fig = plot(f1, f2, f3; layout=(3, 1), size=plot_size(1; height=3 * 0.375),
        leftmargin=8Plots.mm, bottommargin=4Plots.mm)
    return (; fig, panels=(f1, f2, f3))
end

## Plot: Δθ sensitivity (justifies holding Δθ fixed)
"""
Three panels of the Δθ sweep at the optimum: average power, the line-tension
envelope (min–max) and the peak kite acceleration. The x axis is the relative
deviation from the nominal Δθ, and the dashed line marks the nominal value.
"""
function plot_dtheta(dtheta)
    x = [100 * (d.f - 1) for d in dtheta]
    P = [d.P / 1000 for d in dtheta]
    Tmin = [d.Tmin / 1000 for d in dtheta]
    Tmax = [d.Tmax / 1000 for d in dtheta]
    amax = [d.amax for d in dtheta]
    style = (; c=1, m=:circle)
    vline_style = (; c=:gray50, ls=:dash, label="")

    f1 = plot(x, P; style..., label="", ylabel="Power (kW)")
    f2 = plot(x, Tmax; style..., label="max", ylabel="Line tension (kN)",
        fillrange=Tmin, fillalpha=0.1)
    plot!(f2, x, Tmin; style..., label="min")
    f3 = plot(x, amax; style..., label="", ylabel="Max acceleration (m/s²)",
        xlabel="Δθ deviation from nominal (%)")
    foreach(f -> vline!(f, [0]; vline_style...), (f1, f2, f3))
    foreach(f -> plot!(f, xticks=-10:5:10, xlim=(-11, 11)), (f1, f2, f3))
    foreach(f -> plot!(f, xformatter=x -> ""), (f1, f2))

    return plot(f1, f2, f3; layout=(3, 1), size=plot_size(1; height=3 * 0.375),
        leftmargin=9Plots.mm, bottommargin=4Plots.mm)
end

## Plot: shape map over Δθ and Δφ
"""
Three heatmaps of the 10×10 shape map at the optimal ℓ, I, b: average power
(kW), maximum line tension (kN) and peak kite acceleration (m/s²). The white
dashed lines mark the nominal Δθ and the optimal Δφ.
"""
function plot_shape_grid(g, dphi_opt)
    x = 100 .* (g.factors .- 1)
    y = g.dphi_deg
    function pane(Z, title)
        plt = heatmap(x, y, Z; c=:viridis, title=title, titlefontsize=TICKFONTSIZE,
            xticks=-20:10:20, yticks=10:10:90, xlim=extrema(x), ylim=extrema(y),
            xlabel="Δθ deviation (%)", ylabel="Δφ (°)")
        vline!(plt, [0]; c=:white, ls=:dash, lw=1)
        hline!(plt, [dphi_opt]; c=:white, ls=:dash, lw=1)
        scatter!(plt, [0], [dphi_opt]; c=:white, markerstrokecolor=:black, ms=5)
        return plt
    end

    p1 = pane(g.P, "P (kW)")
    p2 = pane(g.Tmax, "Tmax (kN)")
    p3 = pane(g.amax, "max a (m/s²)")
    plot!(p2, ylabel="")
    plot!(p3, ylabel="")
    return plot(p1, p2, p3; layout=(1, 3), size=plot_size(1; height=0.5),
        leftmargin=6Plots.mm, bottommargin=5Plots.mm, topmargin=4Plots.mm)
end

## Plot: objective landscape around the Θ optimum
"""
Plot a normalized landscape from `landscape`: one pane per parameter, two
columns. The x axis is the parameter in its own units, spanning ±`window`
around the optimum, whose value is the middle tick (`-10%`, value, `+10%`).
The cost function (still normalized by its optimum) is a line, solid inside
the box constraint and dotted outside it, the crossing interpolated linearly
between the two samples that bracket the bound. Each box constraint is a black
double vertical line, always drawn (clipped when it falls outside the window)
and on top of the line.
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

## Plot: wind sweep, Δθ, shape map, landscape, time series
my_defaults()
default(xticks=4:2:16, xlim=(4, 16), xformatter=x -> "")

(; fig, panels) = plot_wind_sweep(winds, sweep, i_ref)
display(fig)
savefig(fig, joinpath(OUTPATH, "optimization_wind_sweep_theta.pdf"))
for (f, panel) in zip(panels, ("power", "tension", "params"))
    plot!(f, size=SIZE_FULL, leftmargin=8Plots.mm, bottommargin=6Plots.mm)
    display(f)
    savefig(f, joinpath(OUTPATH, "optimization_$(panel)_theta.pdf"))
end

#= Observation:
Re-optimizing Θ gains +13.5 % at 5 m/s and +8.6 % at 15 m/s, and nothing at
9 m/s by construction. ℓ shortens (30.3 → 26.1 m) and I falls (451 → 326 kg m²)
with wind, while b rises (498 → 1151) and Δφ closes slightly (17.6° → 16.3°).
At 15 m/s Tmax falls from 33.2 to 29.7 kN while Tmin rises from 12.2 to 16.1 kN;
at 5 m/s Tmax rises from 3.34 to 3.49 kN and Tmin falls from 1.99 to 1.72 kN.
=#

my_defaults()
fig_dtheta = plot_dtheta(dtheta)
display(fig_dtheta)
savefig(fig_dtheta, joinpath(OUTPATH, "optimization_dtheta.pdf"))

#= Observation:
Power is flat in Δθ: 2.4 % over ±10 %, with a maximum of 12.279 kW at −4 %
(+0.2 % over nominal). Smaller Δθ lowers Tmin by 5.4 % and raises Tmax by 3.8 %
at −10 %, so the tension spread widens by 14 %; the minimum stays at about half
the maximum. The peak acceleration is U-shaped in Δθ: 1035 m/s² at −10 %,
a 926 m/s² minimum at +2 %, and 1128 m/s² at +10 %.
=#

my_defaults()
fig_shape = plot_shape_grid(shape, rad2deg(p_ref.Δφ))
display(fig_shape)
savefig(fig_shape, joinpath(OUTPATH, "optimization_shape_grid.pdf"))

#= Observation:
At the fixed optimal ℓ, I, b the shape map peaks at Δθ = −20 %, Δφ = 18.9°
with 12.46 kW and collapses away from it: 6.44 kW at Δθ = +20 %, Δφ = 10°,
1.2 kW at Δθ = −20 %, Δφ = 10°, and 2–3 kW for Δφ ≥ 36.7°. Tmax peaks at
14.96 kN at Δθ = −20 %, Δφ = 10° and has a secondary ridge near 12.6 kN at
Δφ = 36.7°; Tmin ranges from 10.5 kN there down to 0.13 kN at Δθ = +20 %,
Δφ = 81°. The peak acceleration reaches 4731 m/s² at Δθ = +20 %, Δφ = 10° and
falls to 129 m/s² at Δθ = +20 %, Δφ = 81°. Δφ = 90° has no pumping cycle at any
Δθ at these fixed ℓ, I, b, so the last row of the map is empty.
=#

my_defaults()
fig_land = plot_landscape(LANDSCAPE)
display(fig_land)
savefig(fig_land, joinpath(OUTPATH, "optimization_landscape.pdf"))

#= Observation:
ℓ is the stiffest direction (−11.5 % at −10 %, −8.1 % at +10 %), I is about
four times flatter (−2.6 %), and b more than twenty times flatter (−0.4 %). The
Δφ optimum (16.6°) is interior now: the objective loses 5.4 % at 15.0° and 2.9 %
at 18.3 %. The box lower bounds on ℓ, I and b sit far from the optimum, so none
of them binds.
=#

my_defaults()
names = (:init, :opt)
sols = map((p_start, p_ref)) do p
    LC.compute_limit_cycle(PMP.build_vbpara(p); sense=SENSE, tol=TOL, save_everystep=true)
end
# savefig only: a second render of the same figure fails Plots' layout pass
for (name, sol) in zip(names, sols)
    let fig = plot_time_series(sol)
        display(fig)
        savefig(fig, joinpath(OUTPATH, "optimization_timeseries_$name.pdf"))
    end
end
fig_ts = plot_time_series_grid(sols)
display(fig_ts)
savefig(fig_ts, joinpath(OUTPATH, "optimization_timeseries_all.pdf"))

## Console report
print_final_table()
