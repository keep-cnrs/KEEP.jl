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
# for v_ref > 2.55 (zoom) / > 24 (full scale) and the long-long sheet for
# tf > 105 s — i.e. where continuation stopped but the arc continues.
#
# Run: julia --project=applications/sprint2026 applications/sprint2026/presentation_figs.jl

using Plots
using Serialization
using Printf
using LinearAlgebra
import BifurcationKit              # only to resolve RecursiveArrayTools on deserialize

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
    for (script, args) in (("_branch_fulltol.jl", ["$(M_FIG)"]),
                           ("_three_cycles_248.jl", String[]))
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
    f = joinpath(SCRATCH, isfile(joinpath(SCRATCH, "three_cycles_shooting.jls")) ?
        "three_cycles_shooting.jls" : "two_cycles_shooting.jls")
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
Split a curve at ONE end against a window `(x_min, y_min, x_max, y_max)`:
samples inside the window are solid, samples outside are dotted, stopping at
the first sample that is back inside — so the middle of the series is never
dotted. `where` is `:begin` (dot the head) or `:end` (dot the tail).
Returns `(solid::UnitRange, dotted::UnitRange)`.
"""
function dot_end(x, y, where::Symbol, window::NTuple{4,Real})
    xmin, ymin, xmax, ymax = window
    inside = (xmin .<= x) .& (x .<= xmax) .& (ymin .<= y) .& (y .<= ymax)
    n = length(x)
    if where === :begin
        k = findfirst(inside)
        return k === nothing ? (1:0, 1:n) : (k:n, 1:k-1)
    elseif where === :end
        k = findlast(inside)
        return k === nothing ? (1:0, 1:n) : (1:k, k+1:n)
    end
    throw(ArgumentError("where must be :begin or :end, got $where"))
end

"Plot `(x,y)` as a solid interior with a dotted end selected by `where`/`window`."
function plot_end!(sp, x, y, where::Symbol, window::NTuple{4,Real};
        color, lw=2.5, label="", dotted_label="")
    solid, dotted = dot_end(x, y, where, window)
    isempty(dotted) || plot!(sp, x[dotted], y[dotted]; color=color, lw=lw, ls=:dot, label=dotted_label)
    isempty(solid) || plot!(sp, x[solid], y[solid]; color=color, lw=lw, label=label)
    return sp
end

"y-position at fraction `f` of a (possibly log) axis with limits `yl`."
yfrac(yl, f) = 10.0^(log10(yl[1]) + f * (log10(yl[2]) - log10(yl[1])))

## ===================================================================== ##
## Fig 1 — the branch and the coexisting flow attractor                    ##
## ===================================================================== ##
"Fold indices of the M=40 branch: fold 1 = p MIN, fold 2 = p MAX on the long sheet."
function branch_layout(M=M_FIG)
    p, tf = collect(br.p), collect(br.tf)
    i0 = argmin(p)
    i2 = argmax(@view p[1:i0])
    return p, tf, i0, i2
end

# NB on units: the branch's `tf` is PHYSICAL (SI) seconds and the state's dα, dτ in
# rad/s. The LimitCycle/Poincaré helpers work in NORMALIZED units (L=2 m, M=6 kg,
# time T0 = l/v_ref s): tf_SI = T_norm * T0, dα_SI = dα_norm / T0.
const X_SHORT_ZOOM = 2.55      # short branch dotted beyond this on the zoom panels
const TF_LL_DOT = 105.0        # long-long sheet dotted above this period
const X_SHORT_FULL = 24.0      # short branch dotted beyond this on the full-scale panel

function fig_fold_zoom(M=M_FIG)
    p, tf, i0, i2 = branch_layout(M)
    p1f, tf1f = p[i0], tf[i0]
    p2f, tf2f = p[i2], tf[i2]

    PANELS = (
        (xlims=(2.44, 2.56), ylims=(0.05, 200.0),
            title="overview — S-curve with both folds"),
        (xlims=(2.44, 2.56), ylims=(5.0, 40.0),
            title="fold 1 zoom — low-wind saddle–node"),
        (xlims=(2.4755, 2.4835), ylims=(40.0, 120.0),
            title="fold 2 zoom — three coexisting periods"),
    )
    plt = plot(layout=(3, 1), size=(880, 1400), legend=:topright,
        plot_title="Limit-cycle branch: two folds ⇒ three periods in (v_ref*₁, v_ref*₂)")

    for (k, spec) in enumerate(PANELS)
        sp = plt[k]
        xl, yl = spec.xlims, spec.ylims
        plot!(sp; xlabel="v_ref [m/s]", ylabel="period  tf [s] (Physical SI, log)",
            xlims=xl, ylims=yl, yscale=:log10, title=spec.title)
        short_win = (xl[1], yl[1], X_SHORT_ZOOM, yl[2])
        ll_win = (xl[1], yl[1], xl[2], TF_LL_DOT)
        # short branch (stable) — dotted right tail (budget) on panels 1–2
        if k == 3
            plot!(sp, p[i0:end], tf[i0:end]; color=C_SHORT, lw=2.5, label="")
        else
            plot_end!(sp, p[i0:end], tf[i0:end], :end, short_win; color=C_SHORT, lw=2.5,
                label=k == 1 ? "short branch, stable" : "",
                dotted_label=k == 1 ? "  continues → (budget)" : "")
        end
        # long sheet (saddle) between the two folds
        plot!(sp, p[i2:i0], tf[i2:i0]; color=C_LONG, lw=2.5,
            label=k == 1 ? "long sheet, saddle" : "")
        # long-long sheet (beyond fold 2) — dotted head (tf > 105 s, budget)
        plot_end!(sp, p[1:i2], tf[1:i2], :begin, ll_win; color=C_LL, lw=2.5,
            label=k == 1 ? "long-long, beyond fold 2" : "")
        # fold tangents (vertical) coloured by the branch each fold belongs to
        vline!(sp, [p1f]; color=C_SHORT, ls=:dash, lw=1.5,
            label=k == 1 ? "fold 1  v_ref*₁ = $(round(p1f, digits=4))" : "")
        vline!(sp, [p2f]; color=C_LONG, ls=:dash, lw=1.5,
            label=k == 1 ? "fold 2  v_ref*₂ = $(round(p2f, digits=4))" : "")
        scatter!(sp, [p1f, p2f], [tf1f, tf2f]; color=:black, ms=6, msw=1.2,
            markerstrokecolor=:white, label="")
        vline!(sp, [TWO_P]; color=:black, ls=:dashdot, lw=1.0,
            label=k == 1 ? "anchor  v_ref = $(round(TWO_P, digits=2))" : "")
        annotate!(sp, p1f + 0.004, yfrac(yl, 0.60),
            text("fold 1\n(tf=$(round(tf1f, digits=1)) s)", :left, 8, :gray30))
        annotate!(sp, p2f + 0.0008, yfrac(yl, 0.87),
            text("fold 2\n(tf=$(round(tf2f, digits=1)) s)", :left, 8, :gray30))
        if k == 1
            annotate!(sp, 2.505, yfrac(yl, 0.22), text("short branch →\n(budget past 2.55)", :left, 8, C_SHORT))
            annotate!(sp, 2.4705, yfrac(yl, 0.96),
                text("long-long continues ↑\n(past tf=$(round(TF_LL_DOT, digits=0)) s, budget)", :left, 8, C_LL))
        end
        if k == 3
            annotate!(sp, 2.4758, yfrac(yl, 0.75),
                text("anchor = 2.480 cuts\n3 times → 3 periods", :left, 8, :gray30))
        end
    end
    return plt
end

## ===================================================================== ##
## Fig 1b — full-scale branch, v_ref ∈ [0, 25]                              ##
## ===================================================================== ##
function fig_full_scale(M=M_FIG)
    p, tf, i0, i2 = branch_layout(M)
    p1f, tf1f = p[i0], tf[i0]
    p2f, tf2f = p[i2], tf[i2]
    xl, yl = (0.0, 25.0), (0.05, 200.0)
    plt = plot(size=(1050, 600), legend=:topright,
        xlims=xl, ylims=yl, yscale=:log10,
        xlabel="v_ref [m/s]", ylabel="period  tf [s] (Physical SI, log)", title="Full-scale limit-cycle branch  (v_ref in [0, 25])", titlelocation=:center,
        left_margin=12Plots.mm, bottom_margin=7Plots.mm)
    plot_end!(plt, p[i0:end], tf[i0:end], :end, (xl[1], yl[1], X_SHORT_FULL, yl[2]);
        color=C_SHORT, lw=2.5, label="short branch, stable", dotted_label="  continues → (budget)")
    plot!(plt, p[i2:i0], tf[i2:i0]; color=C_LONG, lw=2.5, label="long sheet, saddle")
    plot_end!(plt, p[1:i2], tf[1:i2], :begin, (xl[1], yl[1], xl[2], TF_LL_DOT);
        color=C_LL, lw=2.5, label="long-long, beyond fold 2")
    vline!(plt, [p1f]; color=C_SHORT, ls=:dash, lw=1.5, label="fold 1")
    vline!(plt, [p2f]; color=C_LONG, ls=:dash, lw=1.5, label="fold 2")
    vline!(plt, [TWO_P]; color=:black, ls=:dashdot, lw=1.0, label="anchor")
    scatter!(plt, [p1f, p2f], [tf1f, tf2f]; color=:black, ms=6, msw=1.2,
        markerstrokecolor=:white, label="")
    annotate!(plt, 8.0, yfrac(yl, 0.80),
        text("no prograde cycle\nfor v_ref below fold 1", :left, 9, :gray30))
    return plt
end

## ===================================================================== ##
## Fig 2 — phase portraits: all 6 coordinate planes                        ##
## ===================================================================== ##
const LABELS = ["α [rad]", "τ [rad]", "dα [rad/s]", "dτ [rad/s]"]
const PAIRS = [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]

function fig_phase_planes()
    ncyc = longlong === nothing ? 2 : 3
    plt = plot(layout=(2, 3), size=(1250, 760), legend=:topright,
        plot_title="$(ncyc) prograde limit cycles at v_ref = $TWO_P — phase projections")
    for (k, (i, j)) in enumerate(PAIRS)
        plot!(plt[k], Ushort[i, :], Ushort[j, :];
            xlabel=LABELS[i], ylabel=LABELS[j],
            title="$(LABELS[j]) vs $(LABELS[i])",
            color=C_SHORT, lw=2, label="short (stable)",
            left_margin=6Plots.mm, bottom_margin=5Plots.mm)
        plot!(plt[k], Ulong[i, :], Ulong[j, :];
            color=C_LONG, lw=2.2, label="long (saddle)")
        longlong === nothing || plot!(plt[k], Ulonglong[i, :], Ulonglong[j, :];
            color=C_LL, lw=2.2, label="long-long (saddle)")
    end
    return plt
end

## ===================================================================== ##
## Fig 3 — timeseries: 4 components, three lines each                      ##
## ===================================================================== ##
function fig_timeseries()
    ncyc = longlong === nothing ? 2 : 3
    plt = plot(layout=(2, 2), size=(1150, 700), legend=:topright,
        plot_title="Prograde BVP solutions vs physical time at v_ref = $TWO_P ($ncyc cycles)")
    for i in 1:4
        plot!(plt[i], tshort, Ushort[i, :];
            xlabel="t [s]", ylabel=LABELS[i], title=LABELS[i],
            color=C_SHORT, lw=2, label="short  tf=$(round(short.tf, digits=1)) s",
            left_margin=6Plots.mm, bottom_margin=5Plots.mm)
        plot!(plt[i], tlong, Ulong[i, :];
            color=C_LONG, lw=2.2, label="long   tf=$(round(long.tf, digits=1)) s")
        longlong === nothing || plot!(plt[i], tlonglong, Ulonglong[i, :];
            color=C_LL, lw=2.2, label="long-long  tf=$(round(longlong.tf, digits=1)) s")
    end
    return plt
end

for (f, fn) in (("fig_fold_zoom.png", fig_fold_zoom),
                ("fig_full_scale.png", fig_full_scale),
                ("fig_phase_planes.png", fig_phase_planes),
                ("fig_timeseries.png", fig_timeseries))
    savefig(fn(), joinpath(HERE, f))
    println("saved ", f)
end
