# Presentation assets for the low-wind fold result (ECC2026 / sprint2026).
#
# Data pipeline (flag below):
#   USE_SERIALIZED = false (default): runs the two generator scripts
#     scratch/_branch_fulltol.jl 40  -> scratch/brS_shoot_M40_fulltol.jls
#     scratch/_three_cycles_248.jl   -> scratch/three_cycles_shooting.jls
#     then loads them. Reproducible from a fresh clone (~2.5 min).
#   USE_SERIALIZED = true: loads those scratch/*.jls caches only (a few seconds).
#
# Produces (applications/sprint2026/):
#   fig_fold_zoom.png      two-fold branch (overview + fold-1 + fold-2 zooms)
#   fig_full_scale.png     the same branch at full scale, v_ref ∈ [0, 25]
#   fig_phase_planes.png   6 pairwise phase-space projections, the 3 cycles
#   fig_timeseries.png     4 state components vs physical time, the 3 cycles
#
# Dotted convention: a curve is dotted where it leaves the "solid window" set
# for its panel (see dot_end). Windows are chosen so the short branch is dotted
# for v_ref > 2.53 (zoom) / > 24 (full scale), the long branch for tf > 27 s on
# the fold-1 zoom, and the long-long branch for tf > 105 s — i.e. where the curve
# continues past the point we chose to show.
#
# Run: julia --project=applications/sprint2026 applications/sprint2026/presentation_figs.jl

using Plots
using Serialization
using Printf
using LinearAlgebra
using BifurcationKit: BifurcationKit              # only to resolve RecursiveArrayTools on deserialize

const HERE = @__DIR__
const SCRATCH = joinpath(HERE, "scratch")

# Turn ON only while iterating on plotting code: reuse the scratch/*.jls caches
# instead of recomputing the branch and cycles. Off = full reproducible run.
const USE_SERIALIZED = false

const M_FIG = 40
const VREF_FIG = 2.48

const C_SHORT = :royalblue
const C_LONG = :crimson
const C_LL = :darkorange

## ===================================================================== ##
## Data                                                                    ##
## ===================================================================== ##
"Run the two pipeline generators (each writes its scratch/*.jls)."
function generate_caches()
    julia = Base.julia_cmd()
    for (script, args) in
        (("_branch_fulltol.jl", ["$(M_FIG)"]), ("_three_cycles_248.jl", String[]))
        println("running scratch/", script, " ...")
        run(`$julia --project=$HERE $(joinpath(SCRATCH, script)) $args`)
    end
end
USE_SERIALIZED || generate_caches()

"Full-tol branch for arc count M, preferring the deepest (largest-budget) run."
function load_branch(M)
    for tag in ("_fulltol_max800", "_fulltol_max3000", "_fulltol")
        f = joinpath(SCRATCH, "brS_shoot_M$(M)$(tag).jls")
        isfile(f) && return deserialize(f)
    end
    return deserialize(joinpath(SCRATCH, "brS_shoot.jls"))   # legacy archive
end

"Three coexisting cycles at the anchor; falls back to the older 2.47 pair."
function load_cycles()
    f = joinpath(
        SCRATCH,
        if isfile(joinpath(SCRATCH, "three_cycles_shooting.jls"))
            "three_cycles_shooting.jls"
        else
            "two_cycles_shooting.jls"
        end,
    )
    cyc = deserialize(f)
    vr = hasproperty(cyc, :v_ref) ? cyc.v_ref : 2.47
    ll = hasproperty(cyc, :longlong) ? cyc.longlong : nothing
    println("cycles: ", f, "  (v_ref = ", vr, ", ", ll === nothing ? 2 : 3, " cycles)")
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
        C = boundary_point((x[k - 1], y[k - 1]), (x[k], y[k]), window)
        sx = vcat(C[1], float.(x[k:n]))
        sy = vcat(C[2], float.(y[k:n]))
        dx = vcat(float.(x[1:(k - 1)]), C[1])
        dy = vcat(float.(y[1:(k - 1)]), C[2])
        @assert dx[end] == sx[1] && dy[end] == sy[1]       # runs meet at C
    elseif where === :end
        k = findlast(inside)
        if k === nothing
            return (Float64[], Float64[], collect(float.(x)), collect(float.(y)))
        elseif k == n
            return (collect(float.(x)), collect(float.(y)), Float64[], Float64[])
        end
        C = boundary_point((x[k], y[k]), (x[k + 1], y[k + 1]), window)
        sx = vcat(float.(x[1:k]), C[1])
        sy = vcat(float.(y[1:k]), C[2])
        dx = vcat(C[1], float.(x[(k + 1):n]))
        dy = vcat(C[2], float.(y[(k + 1):n]))
        @assert sx[end] == dx[1] && sy[end] == dy[1]
    else
        throw(ArgumentError("where must be :begin or :end, got $where"))
    end
    return sx, sy, dx, dy
end

"Plot `(x,y)` as a solid interior with a dotted end selected by `where`/`window`."
function plot_end!(
    sp,
    x,
    y,
    where::Symbol,
    window::NTuple{4,Real};
    color,
    lw=2.5,
    label="",
    dotted_label="",
)
    sx, sy, dx, dy = dot_end(x, y, where, window)
    isempty(dx) || plot!(sp, dx, dy; color=color, lw=lw, ls=:dot, label=dotted_label)
    isempty(sx) || plot!(sp, sx, sy; color=color, lw=lw, label=label)
    return sp
end

## ===================================================================== ##
## Fig 1 — the branch and the coexisting flow attractor                    ##
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
# time T0 = l/v_ref s): tf_SI = T_norm * T0, dα_SI = dα_norm / T0.
const X_SHORT_ZOOM = 2.53      # short branch dotted beyond this on the zoom panels
const X_SHORT_FULL = 24.0      # short branch dotted beyond this on the full-scale panel
const TF_LONG_DOT = 27.0      # long branch dotted above this on the fold-1 zoom
const TF_LL_DOT = 105.0        # long-long branch dotted above this period

function fig_fold_zoom(M=M_FIG)
    p, tf, i0, i2 = branch_layout(M)
    p1f, tf1f = p[i0], tf[i0]
    p2f, tf2f = p[i2], tf[i2]

    plt = plot(;
        layout=(3, 1),
        size=(900, 1420),
        legend=:topright,
        plot_title="Limit-cycle branch: two folds ⇒ three periods in (v_ref*₁, v_ref*₂)",
    )

    ## panel 1 — overview ---------------------------------------------------
    sp = plt[1]
    xl, yl = (2.44, 2.54), (0.0, 120.0)
    plot!(
        sp;
        xlabel="v_ref (m/s)",
        ylabel="Period (s)",
        xlims=xl,
        ylims=yl,
        title="Overview: S-curve with both folds",
        left_margin=11Plots.mm,
    )
    plot_end!(
        sp,
        p[i0:end],
        tf[i0:end],
        :end,
        (xl[1], yl[1], X_SHORT_ZOOM, yl[2]);
        color=C_SHORT,
        lw=2.5,
        label="short, stable",
        dotted_label="continues",
    )
    plot!(sp, p[i2:i0], tf[i2:i0]; color=C_LONG, lw=2.5, label="long, saddle")
    plot_end!(
        sp,
        p[1:i2],
        tf[1:i2],
        :begin,
        (xl[1], yl[1], xl[2], TF_LL_DOT);
        color=C_LL,
        lw=2.5,
        label="long-long, saddle",
    )
    vline!(
        sp,
        [p1f];
        color=C_SHORT,
        ls=:dash,
        lw=1.5,
        label="fold 1  v_ref*₁ = $(round(p1f, digits=4))",
    )
    vline!(
        sp,
        [p2f];
        color=C_LONG,
        ls=:dash,
        lw=1.5,
        label="fold 2  v_ref*₂ = $(round(p2f, digits=4))",
    )
    scatter!(
        sp,
        [p1f, p2f],
        [tf1f, tf2f];
        color=:black,
        ms=6,
        msw=1.2,
        markerstrokecolor=:white,
        label="",
    )
    vline!(
        sp,
        [TWO_P];
        color=:black,
        ls=:solid,
        lw=1.2,
        label=@sprintf("anchor  v_ref = %.4f", TWO_P)
    )
    annotate!(
        sp,
        p1f + 0.0035,
        0.40 * yl[2],
        text("fold 1\n(tf=$(round(tf1f, digits=1)) s)", :left, 8, :gray30),
    )
    annotate!(
        sp,
        p2f + 0.0010,
        0.66 * yl[2],
        text("fold 2\n(tf=$(round(tf2f, digits=1)) s)", :left, 8, :gray30),
    )
    annotate!(sp, 2.492, 20.0, text("short, stable\n(continues →)", :left, 8, C_SHORT))
    annotate!(
        sp,
        2.4455,
        0.93 * yl[2],
        text(
            "long-long continues ↑\n(past tf=$(round(TF_LL_DOT, digits=0)) s, budget)",
            :left,
            8,
            C_LL,
        ),
    )

    ## panel 2 — fold 1 zoom -------------------------------------------------
    sp = plt[2]
    xl, yl = (2.44, 2.54), (0.0, 30.0)
    plot!(
        sp;
        xlabel="v_ref (m/s)",
        ylabel="Period (s)",
        xlims=xl,
        ylims=yl,
        title=@sprintf("fold 1 — stable + saddle  (v_ref = %.3f, tf = %.1f s)", p1f, tf1f),
        left_margin=11Plots.mm,
    )
    plot_end!(
        sp,
        p[i0:end],
        tf[i0:end],
        :end,
        (xl[1], yl[1], X_SHORT_ZOOM, yl[2]);
        color=C_SHORT,
        lw=2.5,
    )
    plot_end!(
        sp,
        p[i2:i0],
        tf[i2:i0],
        :begin,
        (xl[1], yl[1], xl[2], TF_LONG_DOT);
        color=C_LONG,
        lw=2.5,
    )
    vline!(sp, [p1f]; color=C_SHORT, ls=:dash, lw=1.5, label="")
    scatter!(
        sp, [p1f], [tf1f]; color=:black, ms=6, msw=1.2, markerstrokecolor=:white, label=""
    )
    vline!(sp, [TWO_P]; color=:black, ls=:solid, lw=1.2, label="")

    ## panel 3 — fold 2 zoom -------------------------------------------------
    sp = plt[3]
    xl, yl = (2.4755, 2.4835), (40.0, 120.0)
    plot!(
        sp;
        xlabel="v_ref (m/s)",
        ylabel="Period (s)",
        xlims=xl,
        ylims=yl,
        title=@sprintf("fold 2 — saddle + saddle  (v_ref = %.3f, tf = %.1f s)", p2f, tf2f),
        left_margin=11Plots.mm,
    )
    plot!(sp, p[i0:end], tf[i0:end]; color=C_SHORT, lw=2.5, label="")
    plot!(sp, p[i2:i0], tf[i2:i0]; color=C_LONG, lw=2.5, label="")
    plot_end!(
        sp, p[1:i2], tf[1:i2], :begin, (xl[1], yl[1], xl[2], TF_LL_DOT); color=C_LL, lw=2.5
    )
    vline!(sp, [p2f]; color=C_LONG, ls=:dash, lw=1.5, label="")
    scatter!(
        sp, [p2f], [tf2f]; color=:black, ms=6, msw=1.2, markerstrokecolor=:white, label=""
    )
    vline!(sp, [TWO_P]; color=:black, ls=:solid, lw=1.2, label="")

    return plt
end

## ===================================================================== ##
## Fig 1b — full-scale branch, v_ref ∈ [0, 25]                              ##
## ===================================================================== ##
function fig_full_scale(M=M_FIG)
    p, tf, i0, i2 = branch_layout(M)
    p1f, tf1f = p[i0], tf[i0]
    p2f, tf2f = p[i2], tf[i2]
    xl, yl = (0.0, 25.0), (0.0, 120.0)
    plt = plot(;
        size=(1050, 600),
        legend=:topright,
        xlims=xl,
        ylims=yl,
        xlabel="v_ref (m/s)",
        ylabel="Period (s)",
        title="Full-scale limit-cycle branch  (v_ref in [0, 25])",
        titlelocation=:center,
        left_margin=12Plots.mm,
        bottom_margin=7Plots.mm,
    )
    plot_end!(
        plt,
        p[i0:end],
        tf[i0:end],
        :end,
        (xl[1], yl[1], X_SHORT_FULL, yl[2]);
        color=C_SHORT,
        lw=2.5,
        label="short, stable",
        dotted_label="continues",
    )
    plot!(plt, p[i2:i0], tf[i2:i0]; color=C_LONG, lw=2.5, label="long, saddle")
    plot_end!(
        plt,
        p[1:i2],
        tf[1:i2],
        :begin,
        (xl[1], yl[1], xl[2], TF_LL_DOT);
        color=C_LL,
        lw=2.5,
        label="long-long, saddle",
    )
    vline!(plt, [p1f]; color=C_SHORT, ls=:dash, lw=1.5, label="fold 1")
    vline!(plt, [p2f]; color=C_LONG, ls=:dash, lw=1.5, label="fold 2")
    vline!(
        plt,
        [TWO_P];
        color=:black,
        ls=:solid,
        lw=1.2,
        label=@sprintf("anchor  v_ref = %.4f", TWO_P)
    )
    scatter!(
        plt,
        [p1f, p2f],
        [tf1f, tf2f];
        color=:black,
        ms=6,
        msw=1.2,
        markerstrokecolor=:white,
        label="",
    )
    annotate!(
        plt,
        8.0,
        0.80 * yl[2],
        text("no prograde cycle\nfor v_ref below fold 1", :left, 9, :gray30),
    )
    return plt
end

## ===================================================================== ##
## Fig 2 — phase portraits: all 6 coordinate planes                        ##
## ===================================================================== ##
const LABELS = ["α (rad)", "τ (rad)", "dα (rad/s)", "dτ (rad/s)"]
const PAIRS = [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]

function fig_phase_planes()
    ncyc = longlong === nothing ? 2 : 3
    plt = plot(;
        layout=(2, 3),
        size=(1250, 760),
        legend=:topright,
        plot_title="$(ncyc) prograde limit cycles at v_ref = $TWO_P — phase projections",
    )
    for (k, (i, j)) in enumerate(PAIRS)
        plot!(
            plt[k],
            Ushort[i, :],
            Ushort[j, :];
            xlabel=LABELS[i],
            ylabel=LABELS[j],
            title="$(LABELS[j]) vs $(LABELS[i])",
            color=C_SHORT,
            lw=2,
            label="short, stable",
            left_margin=6Plots.mm,
            bottom_margin=5Plots.mm,
        )
        plot!(plt[k], Ulong[i, :], Ulong[j, :]; color=C_LONG, lw=2.2, label="long, saddle")
        longlong === nothing || plot!(
            plt[k],
            Ulonglong[i, :],
            Ulonglong[j, :];
            color=C_LL,
            lw=2.2,
            label="long-long, saddle",
        )
    end
    return plt
end

## ===================================================================== ##
## Fig 3 — timeseries: 4 components, three lines each                      ##
## ===================================================================== ##
function fig_timeseries()
    ncyc = longlong === nothing ? 2 : 3
    plt = plot(;
        layout=(2, 2),
        size=(1150, 700),
        legend=:topright,
        plot_title="Prograde BVP solutions vs physical time at v_ref = $TWO_P ($ncyc cycles)",
    )
    for i in 1:4
        plot!(
            plt[i],
            tshort,
            Ushort[i, :];
            xlabel="t (s)",
            ylabel=LABELS[i],
            title=LABELS[i],
            color=C_SHORT,
            lw=2,
            label="short, stable  (tf=$(round(short.tf, digits=1)) s)",
            left_margin=6Plots.mm,
            bottom_margin=5Plots.mm,
        )
        plot!(
            plt[i],
            tlong,
            Ulong[i, :];
            color=C_LONG,
            lw=2.2,
            label="long, saddle  (tf=$(round(long.tf, digits=1)) s)",
        )
        longlong === nothing || plot!(
            plt[i],
            tlonglong,
            Ulonglong[i, :];
            color=C_LL,
            lw=2.2,
            label="long-long, saddle  (tf=$(round(longlong.tf, digits=1)) s)",
        )
    end
    return plt
end

for (f, fn) in (
    ("fig_fold_zoom.png", fig_fold_zoom),
    ("fig_full_scale.png", fig_full_scale),
    ("fig_phase_planes.png", fig_phase_planes),
    ("fig_timeseries.png", fig_timeseries),
)
    savefig(fn(), joinpath(HERE, f))
    println("saved ", f)
end
