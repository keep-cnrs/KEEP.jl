# time_series.jl
# One simulation, three panels (states, average power, line tension), shared by
# the numerical-analysis and optimization sections. The legends are drawn by
# hand, one per row, in a dedicated column so every panel keeps the same width
# and the x axes can be linked. In a grid, the y limits of a row are shared by
# all columns.
# Include after jtcam_theme.jl.

using KEEP: PointMass4 as PM4
using KEEP: PointMassPara as PMP
using KEEP: Visualization as VIZ
using KEEP.TorqueFunction: torque_function

using Plots
using PyFormattedStrings: @f_str

"""
Wrapped-τ series: main polyline on the grid with a NaN sample at each 2π
crossing, plus the crossings as NaN-separated vertical segments (crossing
time twice, values π and -π) to be drawn dashed.
"""
function tau_wrap(ts, τm)
    Δt = Base.step(ts)
    tcross = [ts[i] + Δt * (π - τm[i]) / (τm[i+1] + 2π - τm[i])
              for i in findall(diff(τm) .< -π)]
    t_jump, y_jump = Float64[], Float64[]
    for tc in tcross
        append!(t_jump, (tc, tc, NaN))
        append!(y_jump, (π, -π, NaN))
    end
    perm = sortperm(vcat(ts, tcross))
    t_main = vcat(ts, tcross)[perm]
    y_main = vcat(τm, fill(NaN, length(tcross)))[perm]
    return t_main, y_main, t_jump, y_jump
end

"""
Same as `VIZ.plot_avg_power_4D`, but plotting in kW instead of W. Returns the
panel and the y data range, so several simulations can share one power axis.
"""
function plot_avg_power_kW(sol, ts)
    p = sol.prob.p
    L, M, T = PMP.lmt(p)
    t0, tf = first(ts), last(ts)
    dα = sol.(ts, idxs=3)
    power = dα .* torque_function.(dα .* T, (p,)) .* (M * L^2 * T^-2) ./ 1000
    avg_power = sol(tf, idxs=5) / (tf - t0) / 1000
    plt = plot(ts, power; c=:black, label="", legend=false, ylabel="Power (kW)")
    plot!([t0, tf], [avg_power, avg_power]; c=:red, label="")
    plot!(xlabel="t (s)")
    return (; plt, lo=min(0.0, minimum(power)), hi=max(maximum(power), avg_power))
end

"""
Three panels of one simulation `sol` (a transient or a limit cycle): the states
(α, D(α), D(τ), with τ mod 2π), the instantaneous and average power, and the
line tension, whose y ticks are 0 and the two extrema. Returns the panels and
the data ranges used to share the y axes.
"""
function time_series_panels(sol; N=500)
    ts = range(sol.t[1], sol.t[end], N)
    q_ts = reduce(hcat, sol.(ts))'  # interpolated states on the grid
    τm = mod.(q_ts[:, 2], 2π) .- π
    @assert all(-π .<= τm .< π)
    @assert all(diff(τm) .< π)  # τ is monotonically increasing through the wraps
    one_period = abs(abs(q_ts[end, 2] - q_ts[1, 2]) - 2π) <= 1e-3

    t_main, y_main, t_jump, y_jump = tau_wrap(ts, τm)
    q_deg = rad2deg.(q_ts[:, [1, 3, 4]])
    p_wts = plot(ts, q_deg; c=[1 3 4], label="", legend=false,
        ylabel="x (deg, deg/s)", xlabel="", xformatter=x -> "")
    plot!(t_main, rad2deg.(y_main); c=2, label="")
    one_period || plot!(t_jump, rad2deg.(y_jump); c=2, ls=:dot, label=false)

    pw_stats = plot_avg_power_kW(sol, ts)
    pw = pw_stats.plt
    plot!(pw, xlabel="", xformatter=x -> "")

    # Line tension; the ticks are the two extrema and zero. Skip samples where
    # the geometry degenerates (aligned arm and kite), which give a NaN tension
    # and would poison the axis limits.
    tension = PM4.compute_line_tension.(sol.(ts), Ref(sol.prob.p)) ./ 1000
    ok = isfinite.(tension)
    tmin, tmax = extrema(tension[ok])
    p_t = plot(ts[ok], tension[ok]; c=5, label="", legend=false,
        ylabel="Tension (kN)", xlabel="t (s)",
        # yticks=([0, tmin, tmax], ["0", f"{tmin:.1f}", f"{tmax:.1f}"]),
        ylim=(0, tmax), widen=true)

    y_wts = vcat(vec(q_deg), rad2deg.(filter(isfinite, y_main)),
        one_period ? Float64[] : rad2deg.(filter(isfinite, y_jump)))
    return (; wts=p_wts, pw, tension=p_t,
        ylims=(states=extrema(y_wts), power=(pw_stats.lo, pw_stats.hi),
            tension=(0.0, tmax)))
end

"""
Legend tile: one fixed-width box per row (the same axis fraction in every row,
so all three legends match), a handle long enough to read, and a label per
entry. The entry spacing is set in points from `row_px`, the height of one
layout row in pixels, so the entries stay clear of the text at any `plot_size`
(a fixed axis-fraction spacing overlaps at small heights). Plots' GR legend
would instead size the box from its labels and the handle from the subplot
width, which is why the entries are drawn by hand.
"""
function time_series_legend(names, styles, row_px)
    p = plot(; framestyle=:none, showaxis=false, grid=false, legend=false,
        xlims=(0, 1), ylims=(0, 1), xticks=[], yticks=[], margin=0Plots.mm)
    n = length(names)
    dy = 1.25 * LEGENDFONTSIZE / row_px
    ys = 0.5 .+ dy .* ((n + 1) / 2 .- (1:n))
    pad = 0.65dy
    plot!(p, Shape([(0.02, last(ys) - pad), (0.98, last(ys) - pad),
        (0.98, first(ys) + pad), (0.02, first(ys) + pad)]); c=:black, lw=1, fillcolor=:white)
    for (y, style, name) in zip(ys, styles, names)
        plot!(p, [0.08, 0.28], [y, y]; style...)
        annotate!(p, 0.32, y, text(name; family=FONT, pointsize=LEGENDFONTSIZE, halign=:left))
    end
    return p
end

"One hand-drawn legend per row: states, power, tension."
function time_series_legends(row_px)
    return (states=time_series_legend(("α", "τ mod 360°", "D(α)", "D(τ)"),
        ((; c=1), (; c=2), (; c=3), (; c=4)), row_px),
        power=time_series_legend(("P", "Avg(P)"), ((; c=:black), (; c=:red)), row_px),
        tension=time_series_legend(("T",), ((; c=5),), row_px))
end

"Time series of one simulation `sol`, one hand-drawn legend per row."
function plot_time_series(sol; kwargs...)
    return plot_time_series_grid((sol,); kwargs...)  # same layout, so the same legend width
end

"""
Several simulations side by side: one column per solution, one row per panel
(states, power, tension). The hand-drawn legend of each row sits in a dedicated
last column, so every panel keeps the same width; the x axes are linked inside
each column. The y limits of a row are shared by all columns, and every panel
keeps its own y ticks. `legend_width` is the fraction of the figure width given
to the legend column.
"""
function plot_time_series_grid(sols; N=500, legend_width=0.16)
    n = length(sols)
    panels = [time_series_panels(sol; N) for sol in sols]

    # one y range per row, common to every column
    states_lims = (minimum(first(q.ylims.states) for q in panels),
        maximum(last(q.ylims.states) for q in panels))
    power_lims = (0.0, maximum(last(q.ylims.power) for q in panels))
    tension_lims = (0.0, maximum(last(q.ylims.tension) for q in panels))
    for q in panels
        plot!(q.wts, ylims=states_lims, widen=true)
        plot!(q.pw, ylims=power_lims, widen=true)
        plot!(q.tension, ylims=tension_lims, widen=true)
    end

    for q in panels[2:end], p in (q.wts, q.pw, q.tension)
        plot!(p, ylabel="")
    end
    for q in panels
        for p in (q.wts, q.pw, q.tension)
            plot!(p, leftmargin=10Plots.mm)
        end
        plot!(q.tension, bottommargin=6Plots.mm)
    end

    sz = plot_size(4 / 3, 2)
    leg = time_series_legends(sz[2] / 3)
    return plot([q.wts for q in panels]..., leg.states,
        [q.pw for q in panels]..., leg.power,
        [q.tension for q in panels]..., leg.tension;
        layout=grid(3, n + 1; widths=vcat(fill((1 - legend_width) / n, n), legend_width)),
        link=:x, size=sz)
end
