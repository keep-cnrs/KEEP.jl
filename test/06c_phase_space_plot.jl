# Equilibria of the 4D point-mass model, seen in the (α, τ) phase space at rest:
# streamlines of the vector field ddq(q, dq = 0) = (ddα, ddτ), both nullclines
# and the 8 equilibria.
#
# Plots.jl reproduction of the Makie equilibria figure
# (applications/ecc2026/3_equilibria.jl), kept in
# applications/jtcam/Section_Numerical_analysis.jl.

using Plots
using LinearAlgebra: norm

import KEEP.PointMassPara as PMP
using KEEP.SteadyState: ddq_partial, all_steady_states

default(; formatter=:plain, label="")

vbp = PMP.build_vbpara()

"""
Port of Makie's streamplot_impl (src/basic_recipes/streamplot.jl), 2D.

A visited-cell mask covers the domain; seeds are drawn at cell centers from a
quasirandom (plastic-number) low-discrepancy sequence; from each unvisited seed
the NORMALIZED field is RK4-integrated with fixed `stepsize` forward and
backward, stopping on domain exit, on a near-zero field (equilibrium), or on
entering any visited cell (which also terminates closed loops). Repeat until
`density` of the cells are visited.
Seeding and one-arrow-per-streamline follow Makie's streamplot, which steps
explicit Euler on the unit field; the RK4 integration here is a deliberate
deviation, as Euler chatters across the attracting ddτ = 0 nullcline near
equilibria.

Returns (line_x, line_y, ax, ay, au, av): NaN-separated polyline coordinates
and seed positions / unit field directions for the arrows.
"""
function streamlines(
    f, xlims, ylims; gridsize=(20, 20), stepsize=0.1, maxsteps=10, density=1.0
)
    gx, gy = gridsize
    xmin, xmax = xlims
    ymin, ymax = ylims
    rx, ry = LinRange(xmin, xmax, gx + 1), LinRange(ymin, ymax, gy + 1)
    dx, dy = step(rx), step(ry)
    mask = trues(gx, gy)  # unvisited cells
    φ = 1.324717957244746  # plastic number
    acoeff = (φ^-1, φ^-2)
    # Unit field; nothing near an equilibrium (norm < 1e-12), which terminates
    # the integration instead of dividing 0/0.
    unit = x -> begin
        v = f(x)
        n = norm(v)
        n < 1e-12 ? nothing : (v[1] / n, v[2] / n)
    end
    # One RK4 step of size h along (or against, d = -1) the unit field.
    function step_rk4(x, d, h)
        k1 = unit(x)
        k1 === nothing && return nothing
        k2 = unit((x[1] + d * h / 2 * k1[1], x[2] + d * h / 2 * k1[2]))
        k2 === nothing && return nothing
        k3 = unit((x[1] + d * h / 2 * k2[1], x[2] + d * h / 2 * k2[2]))
        k3 === nothing && return nothing
        k4 = unit((x[1] + d * h * k3[1], x[2] + d * h * k3[2]))
        k4 === nothing && return nothing
        return (
            x[1] + d * h / 6 * (k1[1] + 2k2[1] + 2k3[1] + k4[1]),
            x[2] + d * h / 6 * (k1[2] + 2k2[2] + 2k3[2] + k4[2]),
        )
    end
    line_x, line_y = Float32[], Float32[]
    ax, ay, au, av = Float32[], Float32[], Float32[], Float32[]
    visited = 0
    ind = 0
    while visited < gx * gy * min(1.0, density)
        cind = ind
        ind += 1
        ci = clamp(ceil(Int, ((0.5 + acoeff[1] * cind) % 1) * gx), 1, gx)
        cj = clamp(ceil(Int, ((0.5 + acoeff[2] * cind) % 1) * gy), 1, gy)
        mask[ci, cj] || continue
        x0 = (xmin + (ci - 0.5) * dx, ymin + (cj - 0.5) * dy)
        v0 = f(x0)
        n0 = norm(v0)
        push!(ax, x0[1])
        push!(ay, x0[2])
        push!(au, v0[1] / n0)
        push!(av, v0[2] / n0)
        mask[ci, cj] = false
        visited += 1
        for d in (-1, 1)
            x = x0
            ccur = (ci, cj)
            push!(line_x, NaN)
            push!(line_y, NaN)
            push!(line_x, x[1])
            push!(line_y, x[2])
            npts = 1
            while xmin ≤ x[1] ≤ xmax && ymin ≤ x[2] ≤ ymax && npts < maxsteps
                x = step_rk4(x, d, stepsize)
                x === nothing && break
                (xmin ≤ x[1] ≤ xmax && ymin ≤ x[2] ≤ ymax) || break
                idx = (
                    clamp(searchsortedlast(rx, x[1]), 1, gx),
                    clamp(searchsortedlast(ry, x[2]), 1, gy),
                )
                if idx != ccur
                    mask[idx...] || break
                    mask[idx...] = false
                    visited += 1
                    ccur = idx
                end
                push!(line_x, x[1])
                push!(line_y, x[2])
                npts += 1
            end
        end
    end
    return (line_x, line_y, ax, ay, au, av)
end

#! format: off
theme = (;
    downwind        = (; c=plot_color(:grey, 0.08)),
    upwind          = (; c=plot_color(:red, 0.08)),
    streamlines     = (; c=:black, alpha=0.7, lw=0.5),
    heads           = (; c=:black, alpha=0.5, lw=0),
    nullcline_tau   = (; c=1, ls=:solid),
    nullcline_alpha = (; c=2, ls=:solid),
    eq_legend       = (; c=:black, shape=:circle, ms=5, msw=2),
)

# Equilibrium characteristics: the symmetry (α, τ) -> (-α, τ+π) flips the kite
# left↔right (φ = φ0 + Δφ·sin(τ)) and preserves the arm wind region and the
# center/side zone, so the 8 equilibria realize all 2×2×2 combinations.
φkite(τ) = vbp.φ0 + vbp.Δφ * sin(τ)                     # kite azimuth (PM4.τ_to_θφ)
is_down(q) = cos(q[1]) > 0                              # arm downwind
is_center(q) = abs(φkite(q[2]) - vbp.φ0) < vbp.Δφ / 2   # closer to φ0 than to ±Δφ

"Marker style from the equilibrium characteristics: fill = arm wind, shape = zone, stroke = kite side."
ss_style(q) = (
    c = is_down(q) ? :white : :black,
    markerstrokecolor = is_down(q) ? :black : :white,
    shape = is_center(q) ? :circle : :diamond,
    ms = is_down(q) ? 6 : 7,
    markerstrokewidth = is_down(q) ? 2 : 3,
)
#! format: on

f = q -> ddq_partial(q, vbp)

# Canvas & Layout
δ = 0
plot(;
    xlims=(-π - δ, π + δ),
    ylims=(-π - δ, π + δ),
    aspect_ratio=:equal,
    xticks=([-π, 0, π], ["-π", "α = 0", "π"]),
    yticks=([-π, 0, π], ["-π", "τ = 0", "π"]),
    grid=false,
    framestyle=:semi,
    legend=:outerright,
    colorbar=false,
)
# Background wind spans
vspan!([-π, -π / 2]; theme.upwind..., label="Upwind arm")
vspan!([π / 2, π]; theme.upwind...)
vspan!([-π / 2, π / 2]; theme.downwind...) #, label="Downwind arm")

# Streamlines & Arrowheads
line_x, line_y, ax, ay, au, av = streamlines(f, (-π, π), (-π, π))
plot!(line_x, line_y; theme.streamlines...)

L, W = 0.2, 0.06
heads = [
    Shape([
        (s[1] + L / 2 * u, s[2] + L / 2 * v),
        (s[1] - L / 2 * u - W * v, s[2] - L / 2 * v + W * u),
        (s[1] - L / 2 * u + W * v, s[2] - L / 2 * v - W * u),
    ]) for (s, u, v) in zip(zip(ax, ay), au, av)
]
plot!(heads; theme.heads...)

# Nullclines (splatting reuses line properties for both contour and the GR legend fix)
N = 20
αs = range(-π, π, N + 1)
τs = range(-π, π, 2N + 1)
dd = [f((α, τ)) for τ in τs, α in αs]

contour!(αs, τs, last.(dd); levels=[0], theme.nullcline_tau...)
contour!(αs, τs, first.(dd); levels=[0], theme.nullcline_alpha...)
plot!([NaN], [NaN]; theme.nullcline_tau..., label="D²(τ) = 0")
plot!([NaN], [NaN]; theme.nullcline_alpha..., label="D²(α) = 0")

# Equilibria
ss = all_steady_states(vbp)

styles = ss_style.(ss)
scatter!(
    first.(ss),
    last.(ss);
    c=getproperty.(styles, :c),
    markerstrokecolor=getproperty.(styles, :markerstrokecolor),
    shape=getproperty.(styles, :shape),
    ms=getproperty.(styles, :ms),
    markerstrokewidth=getproperty.(styles, :markerstrokewidth),
)
scatter!([NaN], [NaN]; theme.eq_legend..., label="Equilibrium")
plot!(; widen=1.06)

display(current())
