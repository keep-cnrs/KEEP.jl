# Presentation assets for the low-wind fold result (ECC2026 / sprint2026).
#
# Data pipeline (flag below):
#   USE_SERIALIZED = false (default): runs the two generator scripts
#     scratch/_branch_fulltol.jl 40  -> scratch/brS_shoot_M40_fulltol.jls
#     scratch/_three_cycles_248.jl   -> scratch/three_cycles_shooting.jls
#     then loads them. Reproducible from a fresh clone (~2.5 min).
#   USE_SERIALIZED = true: loads those scratch/*.jls caches only (a few seconds).
#
# Produces (applications/sprint2026/), each as .png AND .pdf:
#   fig_full_scale      the whole continuation path, w0 ∈ [0, 25]
#   fig_fold_overview   the two-fold branch, both folds in view
#   fig_fold1           zoom on fold 1 (stable + saddle)
#   fig_fold2           zoom on fold 2 (saddle + saddle)
#   fig_phase_planes    6 pairwise phase-space projections, the 3 cycles
#   fig_timeseries      4 state components vs physical time, the 3 cycles
#
# Dotted convention: a curve is dotted where it leaves the "solid window" set
# for its panel (see dot_end). Windows are chosen so the short branch is dotted
# for w0 > 2.53 (zoom) / > 24 (full scale), the long branch for tf > 27 s on
# the fold-1 zoom, and the long-long branch for tf > 105 s — i.e. where the curve
# continues past the point we chose to show.
#
# Run: julia --project=applications/sprint2026 applications/sprint2026/presentation_figs.jl

using Pkg; Pkg.activate(@__DIR__)

using Plots
using Serialization
using PyFormattedStrings
using LinearAlgebra
import BifurcationKit              # only to resolve RecursiveArrayTools on deserialize

# Presentation canvas: 800x600 for the branch figures (the phase-plane and
# timeseries layouts set their own size). `thickness_scaling` scales lines,
# markers and fonts together; 1.5 keeps an 800x600 panel readable on a slide.
# It is set per branch figure rather than globally so the two layout figures
# keep their previous look.
default(size=(800, 600), widen=true)

const HERE = @__DIR__
const SCRATCH = joinpath(HERE, "scratch")

# Turn ON only while iterating on plotting code: reuse the scratch/*.jls caches
# instead of recomputing the branch and cycles. Off = full reproducible run.
const USE_SERIALIZED = true

const M_FIG = 40
const VREF_FIG = 2.48

const C_SHORT = :royalblue
const C_LONG = :crimson
const C_LL = :darkorange
const C_CONT = :gray25        # the dotted "continues" tail / its legend entry

## ===================================================================== ##
## Data                                                                    ##
## ===================================================================== ##
"Run the two pipeline generators (each writes its scratch/*.jls)."
function generate_caches()
    julia = Base.julia_cmd()
    for (script, args) in (("_branch_fulltol.jl", [f"{M_FIG}"]),
                           ("_three_cycles_248.jl", String[]))
        println("running scratch/", script, " ...")
        run(`$julia --project=$HERE $(joinpath(SCRATCH, script)) $args`)
    end
end
USE_SERIALIZED || generate_caches()

"Full-tol branch for arc count M, preferring the deepest (largest-budget) run."
function load_branch(M)
    for tag in ("_fulltol_max800", "_fulltol_max3000", "_fulltol")
        f = joinpath(SCRATCH, f"brS_shoot_M{M}{tag}.jls")
        isfile(f) && return deserialize(f)
    end
    return deserialize(joinpath(SCRATCH, "brS_shoot.jls"))   # legacy archive
end

"Three coexisting cycles at the anchor; falls back to the older 2.47 pair."
function load_cycles()
    f = joinpath(SCRATCH, isfile(joinpath(SCRATCH, "three_cycles_shooting.jls")) ?
        "three_cycles_shooting.jls" : "two_cycles_shooting.jls")
    cyc = deserialize(f)
    vr = hasproperty(cyc, :w0) ? cyc.w0 : 2.47
    ll = hasproperty(cyc, :longlong) ? cyc.longlong : nothing
    println("cycles: ", f, "  (w0 = ", vr, ", ", ll === nothing ? 2 : 3, " cycles)")
    return cyc.short, cyc.long, ll, vr
end

br = load_branch(M_FIG)
short, long, longlong, TWO_P = load_cycles()
Ushort, Ulong = Array(short.u), Array(long.u)
tshort, tlong = short.t, long.t
if longlong !== nothing
    Ulonglong, tlonglong = Array(longlong.u), Array(longlong.t)
end

## ===================================================================== ##
## Dotted-end helper                                                       ##
## ===================================================================== ##
"""
First boundary crossing of the segment `P0→P1` against the window
`(x_min, y_min, x_max, y_max)`, as a linearly interpolated point. The segment is
assumed to start inside and end outside (or vice versa), so the first edge hit
is the solid/dotted limit.
"""
function boundary_point(P0, P1, window)
    xmin, ymin, xmax, ymax = window
    ts = Float64[]
    for (c0, c1, lo, hi) in ((P0[1], P1[1], xmin, xmax), (P0[2], P1[2], ymin, ymax))
        c1 == c0 && continue
        for e in (lo, hi)
            t = (e - c0) / (c1 - c0)
            0.0 <= t <= 1.0 && push!(ts, t)
        end
    end
    isempty(ts) && return (Float64(P1[1]), Float64(P1[2]))
    t = minimum(ts)
    return (P0[1] + t * (P1[1] - P0[1]), P0[2] + t * (P1[2] - P0[2]))
end

"""
Split a curve at ONE end against a window `(x_min, y_min, x_max, y_max)`:
samples inside the window are solid, samples outside are dotted, stopping at the
first sample that is back inside — so the middle of the series is never dotted.
A linearly interpolated point is inserted **exactly on the boundary** so the
limit is exact even where the samples are sparse, and that crossing point is
shared by both runs (the polyline stays connected; no gap). `where` is `:begin`
(dot the head) or `:end` (dot the tail). Returns `(sx, sy, dx, dy)`.
"""
function dot_end(x, y, where::Symbol, window::NTuple{4,Real})
    xmin, ymin, xmax, ymax = window
    n = length(x)
    n == 0 && return (Float64[], Float64[], Float64[], Float64[])
    inside = (xmin .<= x) .& (x .<= xmax) .& (ymin .<= y) .& (y .<= ymax)
    if where === :begin
        k = findfirst(inside)
        if k === nothing                                   # all dotted
            return (Float64[], Float64[], collect(float.(x)), collect(float.(y)))
        elseif k == 1                                      # nothing dotted
            return (collect(float.(x)), collect(float.(y)), Float64[], Float64[])
        end
        C = boundary_point((x[k-1], y[k-1]), (x[k], y[k]), window)
        sx = vcat(C[1], float.(x[k:n]));   sy = vcat(C[2], float.(y[k:n]))
        dx = vcat(float.(x[1:k-1]), C[1]); dy = vcat(float.(y[1:k-1]), C[2])
        @assert dx[end] == sx[1] && dy[end] == sy[1]       # runs meet at C
    elseif where === :end
        k = findlast(inside)
        if k === nothing
            return (Float64[], Float64[], collect(float.(x)), collect(float.(y)))
        elseif k == n
            return (collect(float.(x)), collect(float.(y)), Float64[], Float64[])
        end
        C = boundary_point((x[k], y[k]), (x[k+1], y[k+1]), window)
        sx = vcat(float.(x[1:k]), C[1]);   sy = vcat(float.(y[1:k]), C[2])
        dx = vcat(C[1], float.(x[k+1:n])); dy = vcat(C[2], float.(y[k+1:n]))
        @assert sx[end] == dx[1] && sy[end] == dy[1]
    else
        throw(ArgumentError("where must be :begin or :end, got $where"))
    end
    return sx, sy, dx, dy
end

"Plot `(x,y)` as a solid interior with a dotted end selected by `where`/`window`."
function plot_end!(sp, x, y, where::Symbol, window::NTuple{4,Real};
        color, lw=2.5, label="", dotted_label="")
    sx, sy, dx, dy = dot_end(x, y, where, window)
    isempty(dx) || plot!(sp, dx, dy; color=color, lw=lw, ls=:dot, label=dotted_label)
    isempty(sx) || plot!(sp, sx, sy; color=color, lw=lw, label=label)
    return sp
end

## ===================================================================== ##
## Fig 1 — the branch: overview, one figure per fold, full scale          ##
## ===================================================================== ##
"Fold indices of the M=40 branch: fold 1 = p MIN, fold 2 = p MAX on the long branch."
function branch_layout(M=M_FIG)
    p, tf = collect(br.p), collect(br.tf)
    i0 = argmin(p)
    i2 = argmax(@view p[1:i0])
    return p, tf, i0, i2
end

# NB on units: the branch's `tf` is PHYSICAL (SI) seconds and the state's dα, dτ in
# rad/s. The LimitCycle/Poincaré helpers work in NORMALIZED units (L=2 m, M=6 kg,
# time T0 = a/w0 s): tf_SI = T_norm * T0, dα_SI = dα_norm / T0.
const X_SHORT_ZOOM = 2.53      # short branch dotted beyond this on the zoom panels
const X_SHORT_FULL = 24.0      # short branch dotted beyond this on the full-scale panel
const TF_LONG_DOT = 27.0      # long branch dotted above this on the fold-1 zoom
const TF_LL_DOT = 105.0        # long-long branch dotted above this period

## ===================================================================== ##
## Fig 1 — the two folds: overview + one figure per fold                   ##
## ===================================================================== ##
function fig_fold_overview(M=M_FIG)
    p, tf, i0, i2 = branch_layout(M)
    p1f, tf1f = p[i0], tf[i0]
    p2f, tf2f = p[i2], tf[i2]

    xl, yl = (2.44, 2.54), (0.0, 120.0)
    sp = plot(; xlabel="Wind speed (m/s)", ylabel="Period (s)", xlims=xl, ylims=yl,
        xticks=[2.44, 2.46, 2.48, 2.50, 2.52], yticks=0:20:120,
        title="S-curve with both folds", legend=:topright, legendfontsize=8,
        thickness_scaling=1.5, left_margin=8Plots.mm)
    plot_end!(sp, p[i0:end], tf[i0:end], :end, (xl[1], yl[1], X_SHORT_ZOOM, yl[2]);
        color=C_SHORT, lw=2.5, label="short, stable")
    plot!(sp, p[i2:i0], tf[i2:i0]; color=C_LONG, lw=2.5, label="long, saddle")
    plot_end!(sp, p[1:i2], tf[1:i2], :begin, (xl[1], yl[1], xl[2], TF_LL_DOT);
        color=C_LL, lw=2.5, label="long-long, saddle")
    plot!(sp, [], []; color=C_CONT, lw=2.5, ls=:dot, label="continues")
    vline!(sp, [p1f]; color=C_SHORT, ls=:dash, lw=1.5,
        label=f"fold 1  w0*₁ = {p1f:.4f}")
    vline!(sp, [p2f]; color=C_LONG, ls=:dash, lw=1.5,
        label=f"fold 2  w0*₂ = {p2f:.4f}")
    scatter!(sp, [p1f, p2f], [tf1f, tf2f]; color=:black, ms=6, msw=1.2,
        markerstrokecolor=:white, label="")
    vline!(sp, [TWO_P]; color=:black, ls=:solid, lw=1.2,
        label=f"anchor  w0 = {TWO_P:.4f}")
    annotate!(sp, 2.492, 20.0, text("short, stable\n(continues →)", :left, 8, C_SHORT))
    annotate!(sp, 2.4405, 0.93 * yl[2],
        text(f"long-long\ncontinues ↑", :left, 8, C_LL))
    return sp
end

"Fold-1 zoom: the stable short cycle meets the long saddle."
function fig_fold1(M=M_FIG)
    p, tf, i0, i2 = branch_layout(M)
    p1f, tf1f = p[i0], tf[i0]

    xl, yl = (2.44, 2.54), (0.0, 30.0)
    sp = plot(; xlabel="Wind speed (m/s)", ylabel="Period (s)", xlims=xl, ylims=yl,
        xticks=[2.44, 2.46, 2.48, 2.50, 2.52], yticks=0:10:30,
        title=f"fold 1: stable + saddle\n(w0, tf) = ({p1f:.3f} m/s, {tf1f:.1f} s)",
        thickness_scaling=1.5, left_margin=8Plots.mm)
    plot_end!(sp, p[i0:end], tf[i0:end], :end, (xl[1], yl[1], X_SHORT_ZOOM, yl[2]);
        color=C_SHORT, lw=2.5)
    plot_end!(sp, p[i2:i0], tf[i2:i0], :begin, (xl[1], yl[1], xl[2], TF_LONG_DOT);
        color=C_LONG, lw=2.5)
    vline!(sp, [p1f]; color=C_SHORT, ls=:dash, lw=1.5, label="")
    scatter!(sp, [p1f], [tf1f]; color=:black, ms=6, msw=1.2, markerstrokecolor=:white, label="")
    vline!(sp, [TWO_P]; color=:black, ls=:solid, lw=1.2, label="")
    return sp
end

"Fold-2 zoom: the long and long-long saddles merge."
function fig_fold2(M=M_FIG)
    p, tf, i0, i2 = branch_layout(M)
    p2f, tf2f = p[i2], tf[i2]

    xl, yl = (2.4755, 2.4835), (40.0, 120.0)
    sp = plot(; xlabel="Wind speed (m/s)", ylabel="Period (s)", xlims=xl, ylims=yl,
        xticks=[2.476, 2.478, 2.480, 2.482], yticks=40:20:120,
        title=f"fold 2: saddle + saddle\n(w0, tf) = ({p2f:.3f} m/s, {tf2f:.1f} s)",
        thickness_scaling=1.5, left_margin=8Plots.mm)
    plot!(sp, p[i0:end], tf[i0:end]; color=C_SHORT, lw=2.5, label="")
    plot!(sp, p[i2:i0], tf[i2:i0]; color=C_LONG, lw=2.5, label="")
    plot_end!(sp, p[1:i2], tf[1:i2], :begin, (xl[1], yl[1], xl[2], TF_LL_DOT);
        color=C_LL, lw=2.5)
    vline!(sp, [p2f]; color=C_LONG, ls=:dash, lw=1.5, label="")
    scatter!(sp, [p2f], [tf2f]; color=:black, ms=6, msw=1.2, markerstrokecolor=:white, label="")
    vline!(sp, [TWO_P]; color=:black, ls=:solid, lw=1.2, label="")
    return sp
end

## ===================================================================== ##
## Fig 1b — full-scale branch, w0 ∈ [0, 25]                              ##
## ===================================================================== ##
function fig_full_scale(M=M_FIG)
    p, tf, i0, i2 = branch_layout(M)
    xl, yl = (0.0, 25.0), (0.0, 120.0)
    plt = plot(legend=:topright, legendfontsize=8,
        xlims=xl, ylims=yl, thickness_scaling=1.5, yticks=0:20:120,
        xlabel="Wind speed (m/s)", ylabel="Period (s)",
        title="Whole continuation path  (w0 in [0, 25])", titlelocation=:center,
        left_margin=8Plots.mm, bottom_margin=7Plots.mm)
    # Only the path: no fold/anchor rules or markers here. Each branch is solid
    # inside the frame and dotted in its own colour past the window; the grey
    # "continues" swatch in the legend is a dummy series (see below).
    win = (xl[1], yl[1], X_SHORT_FULL, yl[2])
    sx, sy, dx, dy = dot_end(p[i0:end], tf[i0:end], :end, win)
    plot!(plt, sx, sy; color=C_SHORT, lw=2.5, label="short, stable")
    plot!(plt, p[i2:i0], tf[i2:i0]; color=C_LONG, lw=2.5, label="long, saddle")
    lx, ly, ldx, ldy = dot_end(p[1:i2], tf[1:i2], :begin,
        (xl[1], yl[1], xl[2], TF_LL_DOT))
    plot!(plt, lx, ly; color=C_LL, lw=2.5, label="long-long, saddle")
    plot!(plt, [], []; color=C_CONT, lw=2.5, ls=:dot, label="continues")
    plot!(plt, dx, dy; color=C_SHORT, lw=2.5, ls=:dot, label="")
    plot!(plt, ldx, ldy; color=C_LL, lw=2.5, ls=:dot, label="")
    return plt
end

## ===================================================================== ##
## Fig 2 — phase portraits: all 6 coordinate planes                        ##
## ===================================================================== ##
const LABELS = ["α (rad)", "τ (rad)", "dα (rad/s)", "dτ (rad/s)"]
const PAIRS = [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]

function fig_phase_planes()
    ncyc = longlong === nothing ? 2 : 3
    plt = plot(layout=(2, 3), size=(1250, 760), legend=:topright,
        plot_title=f"{ncyc} prograde limit cycles at w0 = {TWO_P} m/s — phase projections")
    for (k, (i, j)) in enumerate(PAIRS)
        plot!(plt[k], Ushort[i, :], Ushort[j, :];
            xlabel=LABELS[i], ylabel=LABELS[j],
            title=f"{LABELS[j]} vs {LABELS[i]}",
            color=C_SHORT, lw=2, label="short, stable",
            left_margin=6Plots.mm, bottom_margin=5Plots.mm)
        plot!(plt[k], Ulong[i, :], Ulong[j, :];
            color=C_LONG, lw=2.2, label="long, saddle")
        longlong === nothing || plot!(plt[k], Ulonglong[i, :], Ulonglong[j, :];
            color=C_LL, lw=2.2, label="long-long, saddle")
    end
    return plt
end

## ===================================================================== ##
## Fig 3 — timeseries: 4 components, three lines each                      ##
## ===================================================================== ##
function fig_timeseries()
    ncyc = longlong === nothing ? 2 : 3
    plt = plot(layout=(2, 2), size=(1150, 700), legend=:topright,
        plot_title=f"Prograde BVP solutions vs physical time at w0 = {TWO_P} m/s ({ncyc} cycles)")
    for i in 1:4
        plot!(plt[i], tshort, Ushort[i, :];
            xlabel="t (s)", ylabel=LABELS[i], title=LABELS[i],
            color=C_SHORT, lw=2, label=f"short, stable  (tf={short.tf:.1f} s)",
            left_margin=6Plots.mm, bottom_margin=5Plots.mm)
        plot!(plt[i], tlong, Ulong[i, :];
            color=C_LONG, lw=2.2, label=f"long, saddle  (tf={long.tf:.1f} s)")
        longlong === nothing || plot!(plt[i], tlonglong, Ulonglong[i, :];
            color=C_LL, lw=2.2, label=f"long-long, saddle  (tf={longlong.tf:.1f} s)")
    end
    return plt
end

for (f, fn) in (("fig_full_scale", fig_full_scale),
                ("fig_fold_overview", fig_fold_overview),
                ("fig_fold1", fig_fold1),
                ("fig_fold2", fig_fold2),
                ("fig_phase_planes", fig_phase_planes),
                ("fig_timeseries", fig_timeseries))
    plt = fn()
    for ext in ("png", "pdf")
        savefig(plt, joinpath(HERE, f"{f}.{ext}"))
    end
    println("saved ", f, ".png/.pdf")
end
