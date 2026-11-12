# plot_helpers.jl
# Shared plotting helpers for the JTCAM figures: the paper theme (fonts,
# dimensions, palette, margins), the axis-tick helper `edge_ticks`, and the
# time-series panels with their hand-drawn legends. Merged from jtcam_theme.jl
# and time_series.jl so a section only includes one file.
import Plots
using Plots
using Plots.PlotMeasures
using PyFormattedStrings: @f_str

using KEEP: PointMass4 as PM4
using KEEP: PointMassPara as PMP
using KEEP: Visualization as VIZ
using KEEP.TorqueFunction: torque_function

Plots.gr()

# --- Font Setup ---
function setup_font(font_filename="LinLibertine_R.ttf")
    paths = [
        joinpath(@__DIR__, font_filename),
        joinpath(expanduser("~/Library/Fonts"), font_filename),
        joinpath("/usr/share/fonts/truetype/linux-libertine", font_filename)
    ]

    idx = findfirst(isfile, paths)
    if idx !== nothing
        font_id = Plots.GR.loadfont(abspath(paths[idx]))
        if font_id >= 0
            Plots.gr_font_family["linux libertine"] = font_id
            return "linux libertine"
        end
    end

    @warn "$font_filename not found. Drop it in ~/Library/Fonts/ or beside this script. Falling back to 'times'."
    return "times"
end

const FONT = setup_font("LinLibertine_R.ttf")

# --- Dimensions ---
# LaTeX DPI is 72.27 pt/in; JTCAM textwidth is 155 mm (~441 pt), half-width is 75 mm (~213 pt)

# W = \showthe\textwidth = \showthe\columnwidth = \showthe\hsize
# WIDTH_FULL_INCH = W
const PT_PER_INCH = 72.27  # Latex constant
const DPI = 72  # Plots.jl constant

const WIDTH_PT = 415.41025

# --- Typography ---
const FONTSIZE = 12
const LABELFONTSIZE = FONTSIZE
const TITLEFONTSIZE = FONTSIZE
const TICKFONTSIZE = 10
const LEGENDFONTSIZE = 10
const ANNOTATIONFONTSIZE = 10

# --- Sizing Functions & Global Objects ---
"""
One width unit is half the text width and one height unit 3/8 of it, so
`plot_size(1, 1)` is a half-width 4/3 plot and `plot_size(1, 2)` spans the text
width at the same height.
"""
function plot_size(height=1, width=1; dpi=DPI)
    width_in = width * WIDTH_PT / 2 / PT_PER_INCH
    height_in = height * 3 * WIDTH_PT / 8 / PT_PER_INCH
    width_px = round(Int, width_in * dpi)
    height_px = round(Int, height_in * dpi)
    return (width_px, height_px)
end

const SIZE_FULL = plot_size(1, 2)
const SIZE_HALF = plot_size(1, 1)

# --- Palette ---
const PALETTE = Plots.palette(:auto)
# Standard 8-color Okabe-Ito (Wong) palette
# const PALETTE = okabe_ito = [
#     colorant"#E69F00",  # Orange
#     colorant"#56B4E9",  # Sky Blue
#     colorant"#009E73",  # Bluish Green
#     colorant"#F0E442",  # Yellow
#     colorant"#0072B2",  # Blue
#     colorant"#D55E00",  # Vermilion
#     colorant"#CC79A7",  # Reddish Purple
#     colorant"#000000",  # Black
# ]

# const PALETTE = Plots.palette(:Dark2_8)

# JTCAM primary dark blue RGB(10, 40, 80) with harmonious accents
# const PALETTE = [
#     # Plots.RGB(10/255, 40/255, 80/255),  # Main accent color is too dark
#     Plots.RGB(0/255, 140/255, 185/255),  # Opting for the brighter one
#     Plots.RGB(0.85, 0.35, 0.15),
#     Plots.RGB(0.15, 0.60, 0.45),
#     Plots.RGB(0.70, 0.20, 0.50),
#     Plots.RGB(0.40, 0.40, 0.40),
# ]

# const PALETTE = [
#     RGB(0/255, 140/255, 185/255),  # 1. Cerulean Blue (Your anchor)
#     RGB(235/255, 105/255, 65/255),  # 2. Coral Orange
#     RGB(45/255, 170/255, 110/255),  # 3. Emerald Green
#     RGB(165/255, 90/255, 195/255),  # 4. Amethyst Violet
#     RGB(235/255, 175/255, 30/255),  # 5. Marigold Amber
#     RGB(225/255, 70/255, 125/255),  # 6. Raspberry Rose
#     RGB(30/255, 180/255, 175/255),  # 7. Bright Teal
#     RGB(90/255, 110/255, 195/255),  # 8. Cornflower Indigo
# ]


# --- Publication Defaults ---
function my_defaults()
    Plots.default()
    Plots.default(
        thickness_scaling=0.5,
        dpi=DPI,
        fontfamily=FONT,
        size=SIZE_HALF,
        palette=PALETTE,
        formatter=:plain,
        label="",

        # Font sizes
        labelfontsize=LABELFONTSIZE,
        tickfontsize=TICKFONTSIZE,
        legendfontsize=LEGENDFONTSIZE,
        titlefontsize=TITLEFONTSIZE,
        annotationfontsize=ANNOTATIONFONTSIZE,

        # Line & Marker aesthetics
        linewidth=3,
        markersize=7,
        markerstrokewidth=0,

        # Axes & Grid layout
        framestyle=:box,
        grid=true,
        gridalpha=0.1,
        gridlinewidth=1,
        minorgrid=false,
        tick_direction=:none,

        # Margins & Legends
        margin=0Plots.mm, # 4Plots.mm,
        legend=:topright,
        legend_background_color=Plots.RGBA(1.0, 1.0, 1.0, 0.85),
        widen=true
    )
end

const OUTPATH = mkpath(joinpath(@__DIR__, "figs"))
my_defaults()

# --- Axes ---
"""
First and last of the default ticks of the axis `axis` (`:x`, `:y` or `:z`) of
subplot `sp`, plus zero when `include_zero` and zero lies strictly between them.
"""
function edge_ticks(sp, axis; include_zero=false)
    ticks = first(Plots.get_ticks(sp, sp[Symbol(axis, :axis)]))
    lo, hi = first(ticks), last(ticks)
    return include_zero && lo < 0 < hi ? [lo, 0.0, hi] : [lo, hi]
end

"Tick label without trailing zeros: `0.9` instead of `0.900`, `1` instead of `1.0`."
plain_tick(x) = f"{x:g}"

# --- Time series ---
# One simulation, three panels (states, average power, line tension), shared by
# the numerical-analysis and optimization sections. The legends are drawn by
# hand, one per row, in a dedicated column so every panel keeps the same width
# and the x axes can be linked. In a grid, the y limits of a row are shared by
# all columns.

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
    plt = plot(ts, power; c=:black, xlabel="t (s)", ylabel="Power (kW)", yticks=0:10:100)
    plot!([t0, tf], [avg_power, avg_power]; c=:red)
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
    p_wts = plot(ts, q_deg; c=[1 3 4], ylabel="x (deg, deg/s)", xlabel="", xformatter=x -> "", yticks=-720:180:720)
    plot!(t_main, rad2deg.(y_main); c=2)
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
    tmin = min(0.0, tmin)
    p_t = plot(ts[ok], tension[ok]; c=5, ylabel="Tension (kN)", xlabel="t (s)", yticks=-100:5:100, ylim=(tmin, tmax), widen=true)

    y_wts = vcat(vec(q_deg), rad2deg.(filter(isfinite, y_main)),
        one_period ? Float64[] : rad2deg.(filter(isfinite, y_jump)))
    return (; wts=p_wts, pw, tension=p_t,
        ylims=(states=extrema(y_wts), power=(pw_stats.lo, pw_stats.hi),
            tension=(tmin, tmax)))
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
    p = plot(; framestyle=:none, showaxis=false, grid=false,
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
each column. The y limits of a row are shared by all columns, and only the
leftmost column labels the y ticks (the others keep the ticks, so the gridlines
stay). `legend_width` is the fraction of the figure width given to the legend
column. With `titles` (one per solution) a shallow first row carries a title
over each column; the legend column of that row stays empty.
"""
function plot_time_series_grid(sols; N=500, legend_width=0.22, titles=false)
    n = length(sols)
    @assert titles === false || length(titles) == n "one title per solution"
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
        plot!(p, ylabel="", yformatter=x -> "")  # ticks stay, for the gridlines
    end
    for (j, q) in enumerate(panels)
        for p in (q.wts, q.pw, q.tension)
            # the leftmost column carries the y labels, the others only ticks
            plot!(p, leftmargin=(j == 1 ? 20 : 3)Plots.mm)
        end
        plot!(q.tension, bottommargin=6Plots.mm)
    end

    sz = plot_size(4 / 3, 2)

    # one relative height per row: the first entry sizes the title row
    heights = Float64[1, 1, 1]
    titles === false || prepend!(heights, 0.01)
    leg = time_series_legends(sz[2] / sum(heights))  # one panel row, in px
    widths = vcat(fill((1 - legend_width) / n, n), legend_width)
    heights ./= sum(heights)

    # title row: one hidden panel per column and an empty legend cell
    header = titles === false ? () :
             (map(t -> plot(; title=t, framestyle=:none, showaxis=false, margin=-10Plots.mm),
        titles)..., plot(; framestyle=:none, showaxis=false))

    return plot(header...,
        [q.wts for q in panels]..., leg.states,
        [q.pw for q in panels]..., leg.power,
        [q.tension for q in panels]..., leg.tension;
        layout=grid(length(heights), n + 1; heights=heights, widths=widths),
        link=:x, size=sz)
end
