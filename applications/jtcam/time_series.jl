# time_series.jl
# One simulation, three panels (states, average power, line tension), shared by
# the numerical-analysis and optimization sections. Legends are built in, one
# per row, and a grid of simulations shares the y limits of each row across
# columns.
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
function plot_avg_power_kW(sol; legend=:outerright)
    t, (t0, tf) = VIZ.build_t(extrema(sol.t), VIZ.DEFAULT_PPS)
    p = sol.prob.p
    L, M, T = PMP.lmt(p)
    dα = sol.(t, idxs=3)
    power = dα .* torque_function.(dα .* T, (p,)) .* (M * L^2 * T^-2) ./ 1000
    avg_power = sol(tf, idxs=5) / (tf - t0) / 1000
    plt = plot(t, power; c=:black, label="P", ylabel="Power (kW)", legend=legend)
    plot!([t0, tf], [avg_power, avg_power]; c=:red, label="Avg(P)")
    plot!(xlabel="t (s)")
    return (; plt, lo=min(0.0, minimum(power)), hi=max(maximum(power), avg_power))
end

"""
Three panels of one simulation `sol` (a transient or a limit cycle): the states
(α, D(α), D(τ), with τ mod 2π), the instantaneous and average power, and the
line tension, whose y ticks are 0 and the two extrema. `legend` is forwarded to
each panel: the right-most column of a grid draws it, the other columns pass
`false`. Returns the panels and the data ranges used to share the y axes.
"""
function time_series_panels(sol; N=500, legend=:outerright)
    ts = range(sol.t[1], sol.t[end], N)
    q_ts = reduce(hcat, sol.(ts))'  # interpolated states on the grid
    τm = mod.(q_ts[:, 2], 2π) .- π
    @assert all(-π .<= τm .< π)
    @assert all(diff(τm) .< π)  # τ is monotonically increasing through the wraps

    t_main, y_main, t_jump, y_jump = tau_wrap(ts, τm)
    p_wts = plot(ts, q_ts[:, [1, 3, 4]]; c=[1 3 4],
        label=["α" "D(α)" "D(τ)"], legend=legend,
        ylabel="x (rad, rad/s)", xlabel="", xformatter=x -> "")
    plot!(t_main, y_main; c=2, label="τ mod 2π")
    plot!(t_jump, y_jump; c=2, ls=:dot, label=false)

    pw_stats = plot_avg_power_kW(sol; legend)
    pw = pw_stats.plt
    plot!(pw, xlabel="", xformatter=x -> "")

    # Line tension; the ticks are the two extrema and zero. Skip samples where
    # the geometry degenerates (aligned arm and kite), which give a NaN tension
    # and would poison the axis limits.
    tension = PM4.compute_line_tension.(sol.(ts), Ref(sol.prob.p)) ./ 1000
    ok = isfinite.(tension)
    tmin, tmax = extrema(tension[ok])
    p_t = plot(ts[ok], tension[ok]; c=5, label="T", legend=legend,
        ylabel="Tension (kN)", xlabel="t (s)",
        yticks=([0, tmin, tmax], ["0", f"{tmin:.1f}", f"{tmax:.1f}"]),
        ylim=(0, tmax), widen=true)

    y_wts = vcat(vec(q_ts[:, [1, 3, 4]]), filter(isfinite, y_main),
        filter(isfinite, y_jump))
    return (; wts=p_wts, pw, tension=p_t,
        ylims=(states=extrema(y_wts), power=(pw_stats.lo, pw_stats.hi),
            tension=(0.0, tmax)))
end

"Time series of one simulation `sol`, one built-in legend per row on the right."
function plot_time_series(sol; N=500)
    (; wts, pw, tension) = time_series_panels(sol; N, legend=:outerright)
    # margins go on the panels, not on the composed figure: re-computing the
    # layout of a figure that already carries margins is not stable in Plots
    plot!(wts, leftmargin=12Plots.mm, rightmargin=4Plots.mm)
    plot!(pw, leftmargin=12Plots.mm, rightmargin=4Plots.mm, topmargin=-5Plots.mm)
    plot!(tension, leftmargin=12Plots.mm, rightmargin=4Plots.mm, topmargin=-5Plots.mm,
        bottommargin=4Plots.mm)
    return plot(wts, pw, tension; layout=@layout([a{0.42h}; b; c]),
        size=plot_size(1; height=1.5))
end

"""
Several simulations side by side: one column per solution, one row per panel
(states, power, tension). The right-most panel of each row carries the built-in
legend, all columns of a row share the same y limits, and every panel keeps its
own y ticks. Only the first column carries the y axis labels.
"""
function plot_time_series_grid(sols; N=500)
    n = length(sols)
    panels = [time_series_panels(sol; N, legend=(j == n ? :outerright : false))
              for (j, sol) in enumerate(sols)]

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
            plot!(p, leftmargin=12Plots.mm)
        end
        plot!(q.tension, bottommargin=4Plots.mm)
    end
    return plot([q.wts for q in panels]..., [q.pw for q in panels]...,
        [q.tension for q in panels]...;
        layout=grid(3, n), size=plot_size(1; height=1.5))
end
