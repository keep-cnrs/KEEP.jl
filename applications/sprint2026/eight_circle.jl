# Subplot-aware re-implementation of the "eight + circle" diagram from
# KEEP.Visualization, for the three-limit-cycle animation.
#
# The src version draws to the CURRENT plot and returns (fig, plot_vbp), so it
# cannot target a subplot (forcing fragile Plot-object combining) and a trail
# costs one series per state. Here the geometry is precomputed and every draw
# takes an explicit subplot handle; a trail is a single NaN-joined polyline.
#
# Eight geometry. θ is measured from the zenith (θ=0 zenith, θ=π/2 horizontal),
# so in a y/z plot altitude DECREASES with θ and the lemniscate's vertical
# modulation is the negative of the naive (φ,θ) one:
#     eight_x(τ) = Δφh·sin(τ)
#     eight_y(τ) = (l + r) − Δθv·sin(2τ)      centre = l + r
# with Δφh = (l+r)·Δφ, Δθv = (l+r)·Δθ.
#
# Modes:
#   :src       true scale (above), arm circle radius arm_scaling·l at the origin.
#   :rescaled  shape normalized to half-width 1, centred at y8 — SCHEMATIC,
#              relative arm/eight scale not preserved.
#
# `arm_scaling` (default 1) scales ONLY the arm (circle radius and kite
# position); the eight is untouched.
#
# State q = (α, τ) in radians (same convention as the BVP / LimitCycle).
module EightCircle

using Plots
import KEEP.PointMass4 as PM4
using KEEP.PointMassPara: build_para

export Geometry,
    eight_circle_geometry, draw_background!, add_state!, add_trail!, kite_xy, eight_xy

struct Geometry
    mode::Symbol
    arm_scaling::Float64
    l::Float64
    θ0::Float64          # = l + r : eight vertical centre and radial scale
    Δφh::Float64         # cross-range half-amplitude
    Δθv::Float64         # vertical half-amplitude
    y8::Float64          # vertical placement of the (:rescaled) eight
    circle_x::Vector{Float64}
    circle_y::Vector{Float64}
    eight_x::Vector{Float64}
    eight_y::Vector{Float64}
end

function eight_circle_geometry(
    vbp; mode::Symbol=:src, arm_scaling::Real=1.0, y8::Real=2.1, n::Int=400
)
    mode in (:src, :rescaled) || error("mode must be :src or :rescaled")
    l, r = vbp.l, vbp.r
    θ0 = l + r
    Δφh = θ0 * vbp.Δφ
    Δθv = vbp.Δθ * θ0

    ls = arm_scaling * (mode === :src ? l : 1.0)
    αs = range(-π, π; length=n)
    circle_x = ls .* sin.(αs)
    circle_y = ls .* cos.(αs)

    τs = range(-π, π; length=n)
    if mode === :src
        eight_x = Δφh .* sin.(τs)
        eight_y = θ0 .- Δθv .* sin.(2 .* τs)       # flipped vertical modulation
    else
        eight_x = sin.(τs)
        eight_y = y8 .- (Δθv / Δφh) .* sin.(2 .* τs)
    end
    return Geometry(
        mode,
        Float64(arm_scaling),
        l,
        θ0,
        Δφh,
        Δθv,
        Float64(y8),
        circle_x,
        circle_y,
        eight_x,
        eight_y,
    )
end

"Kite plot point (xdata, ydata) for arm angle α."
function kite_xy(g::Geometry, α)
    a = PM4.compute_αhat(α)              # SA[cos α, sin α, 0]
    s = g.arm_scaling * (g.mode === :src ? g.l : 1.0)
    return (s * a[2], s * a[1])          # (l sinα, l cosα)
end

"Eight plot point (xdata, ydata) for phase τ."
function eight_xy(g::Geometry, τ)
    if g.mode === :src
        return (g.Δφh * sin(τ), g.θ0 - g.Δθv * sin(2τ))
    else
        return (sin(τ), g.y8 - (g.Δθv / g.Δφh) * sin(2τ))
    end
end

"Draw the static circle + eight background into subplot `sp`."
function draw_background!(sp, g::Geometry; color=:gray55, lw=1.2)
    plot!(sp, g.circle_x, g.circle_y; color=color, lw=lw, label="")
    plot!(sp, g.eight_x, g.eight_y; color=color, lw=lw, label="")
    return sp
end

"Draw one state as a spoke: origin → kite → eight point."
function add_state!(sp, q, g::Geometry; kwargs...)
    k = kite_xy(g, q[1])
    e = eight_xy(g, q[2])
    return plot!(sp, [0.0, k[1], e[1]], [0.0, k[2], e[2]]; kwargs...)
end

"""
Draw a whole trajectory `qs` of states as ONE NaN-joined polyline of spokes.
Constant style only (no per-point alpha); use per-state `add_state!` for a ramp.
"""
function add_trail!(sp, qs, g::Geometry; kwargs...)
    xs = Float64[]
    ys = Float64[]
    for q in qs
        k = kite_xy(g, q[1])
        e = eight_xy(g, q[2])
        append!(xs, (0.0, k[1], e[1], NaN))
        append!(ys, (0.0, k[2], e[2], NaN))
    end
    return plot!(sp, xs, ys; kwargs...)
end

end  # module
