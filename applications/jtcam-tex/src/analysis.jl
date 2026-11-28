# Section 2 figures and Table 1: equilibria, one simulation and the convergence of
# the two limit cycles of the 4D point-mass model, plus the parameter table.
#
# Outputs: figs/{equilibria_phase_space,equilibria_eight_circle,time_series_4D,
# limit_cycle_convergence,yz_limit_cycles}.{pdf,png} and tables/table1.tex.
#
# Not part of the paper build (kept in applications/jtcam/Section_Numerical_analysis.jl):
# the eigenvalue and equilibrium-trajectory figures, the Floquet check of the
# convergence slopes and the limit-cycle inventory (deflated shooting vs
# quasi-Monte-Carlo), which print diagnostics only.
using KEEP: PointMass4 as PM4
using KEEP: PointMassPara as PMP
using KEEP: SteadyState
using KEEP: Visualization as VIZ
using KEEP: LimitCycle as LC
using KEEP: DEFAULT_TOLERANCE, TAU0

using Plots
using LinearAlgebra: norm, eigvals, eigen, ⋅
import ForwardDiff
using SplitApplyCombine: invert
using QuasiMonteCarlo: sample, HaltonSample
using PyFormattedStrings: @f_str
using Roots
const SA = PM4.SA  # StaticArrays, through KEEP

include(joinpath(@__DIR__, "helpers.jl"))

## Parameter table for the paper
# LaTeX table of the parameter set used here, sorted alphabetically by symbol.
# b is the optimized PTO slope (kinetic criterion); angles are in degrees
# and b per degree, so b·dα is a torque when dα is read in deg/s.

struct Parameter
    key::Symbol
    sym::String
    descr::String
    unit::String
    value_override::Union{String,Nothing}
end

const LATEX_TO_UNICODE = Dict(
    raw"\ell" => "l",
    raw"\rho" => "ρ",
    raw"\theta" => "θ",
    raw"\Delta" => "Δ",
    raw"\phi" => "φ"
)

function eor_transform(sym)
    s = sym
    for (k, v) in LATEX_TO_UNICODE
        s = replace(s, k => v)
    end
    s = replace(s, raw"_" => "", raw"{" => "", raw"}" => "", raw"\\" => "")
    return lowercase(s)
end

function param_table(p)
    # symbol, description, unit, optional value override
    params = [
        Parameter(:a, "a", "Arm length", raw"\meter", nothing),
        Parameter(:l, raw"\ell", "Line length", raw"\meter", nothing),
        Parameter(:m, "m", "Kite mass", raw"\kg", nothing),
        Parameter(:I, "I", "Arm inertia", raw"\kg\square\meter", nothing),
        Parameter(:ρ_a, raw"\rho_a", "Air density", raw"\kg\per\cubic\meter", nothing),
        Parameter(:d_l, raw"d_\ell", "Line diameter", raw"\meter", nothing),
        Parameter(:S, "S", "Kite surface", raw"\square\meter", nothing),
        Parameter(:C_L, "C_L", "Kite lift coefficient", "", nothing),
        Parameter(:C_D, "C_D", "Kite drag coefficient", "", "C_L / 7"),
        Parameter(:C_D_l, raw"C_{D,\ell}", "Line drag coefficient", "", nothing),
        Parameter(:ρ_l, raw"\rho_\ell", "Line density", raw"\kg\per\cubic\meter", nothing),
        Parameter(:g, "g", "Gravity", raw"\meter\per\square\second", nothing),
        Parameter(:z0, "z_0", "Wind reference altitude", raw"\meter", nothing),
        Parameter(:w0, "w_0", "Wind reference speed", raw"\meter\per\second", nothing),
        Parameter(:b, "b", "Arm braking coefficient",
            raw"\newton\meter\second\per\degree", f"{p.b / rad2deg(1):.9g}"),
        Parameter(:θ0, raw"\theta_0", "Lemniscate mean polar angle",
            raw"\degree", f"{rad2deg(p.θ0):.0f}"),
        Parameter(:Δθ, raw"\Delta\theta", "Lemniscate polar half-height",
            raw"\degree", f"{rad2deg(p.Δθ):.0f}"),
        Parameter(:Δφ, raw"\Delta\phi", "Lemniscate azimuthal half-width",
            raw"\degree", f"{rad2deg(p.Δφ):.0f}"),
    ]

    sort!(params, by=x -> eor_transform(x.sym))
    table_rows = map(params) do p_item
        val = isnothing(p_item.value_override) ? f"{p[p_item.key]:.9g}" : p_item.value_override
        unit_col = isempty(p_item.unit) ? "" : raw"\si{" * p_item.unit * "}"
        f"    ${p_item.sym}$ & {p_item.descr} & ${val}$ & {unit_col} \\\\ "
    end


    table_str = raw"""
\begin{table}[thpb]
\centering
\caption{Parameter Description and Values}
\label{tab:params}
\begin{tabular}{l p{3.2cm} r l}
    \toprule
    {Symbol} & {Description} & {Value} & {Unit} \\
    \midrule
""" * join(table_rows, "\n") * raw"""

    \bottomrule
\end{tabular}
\end{table}"""
    return table_str
end


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
function streamlines(f, xlims, ylims; gridsize=(20, 20), stepsize=0.1, maxsteps=10, density=1.0)
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
        (x[1] + d * h / 6 * (k1[1] + 2k2[1] + 2k3[1] + k4[1]),
            x[2] + d * h / 6 * (k1[2] + 2k2[2] + 2k3[2] + k4[2]))
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
                idx = (clamp(searchsortedlast(rx, x[1]), 1, gx), clamp(searchsortedlast(ry, x[2]), 1, gy))
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

## Convergence rate of the two cycles (helpers)
const LC_CONV_TOL = 1e-13
const CB = LC.build_poincare_callback(eps(1.0))

function get_slope(distances; tol_start, tol_end)
    i1 = findfirst(<(tol_start), distances)
    i2 = findlast(>(tol_end), distances)
    return -(log10(distances[i2]) - log10(distances[i1])) / (i2 - i1)
end

function plot_convergence_panel(p, dτ, color; tf=50, ntrace=5, legend=true, tol=LC_CONV_TOL)
    sol = PM4.integrate(PM4.build_u(; dτ), tf, p; callback=CB, save_end=false, tol)
    u_infty = sol.u[end]
    n_max = 11
    distances = [LC.distance_on_section(u, u_infty) for u in sol.u[1:(end-1)]][1:min(end, n_max-1)]

    tol_start, tol_end = 1e-2, 10tol
    slope = get_slope(distances; tol_start, tol_end)
    i1 = findfirst(<(tol_start), distances)

    plt = plot(; yscale=:log10, xlim=(0, n_max), xticks=0:2:100,
        ylim=(1e-16, 1e2), yticks=10.0 .^ (0:-3:-16),
        xlabel=f"""Lemniscate revolutions, D(τ) {dτ>0 ? ">" : "<"} 0""",
    )

    perturbed_params = (; ls=:dot, lw=2, alpha=0.4, c=:gray, label="Perturbed trajectories")
    plot!(i -> distances[i1] * exp10(-slope * (i - i1)); c=:black)
    if ntrace > 0
        for (α, τ) in eachcol(sample(ntrace, Float64[-π, -π], Float64[π, π], HaltonSample()))
            s = PM4.integrate(PM4.build_u(; α, τ, dτ), tf, p; callback=CB, save_end=false, tol)
            plot!([LC.distance_on_section(u, u_infty) for u in s.u] .+ 1e-16;
                perturbed_params..., label="")
        end
    end

    annotate!((5, 1e-11, (f"slope = {slope:.2f}", :right, ANNOTATIONFONTSIZE, colorant"black")))
    distance_style = (; c=color, ls=:dot, marker=:circle)
    plot!(distances .+ 1e-16; distance_style...)

    if legend
        plot!([NaN]; distance_style..., label=raw"$|| x_i - x_\infty  ||$")
        ntrace > 0 && plot!([NaN]; perturbed_params...)
    end
    return (; plt, slope)
end

## Spectra of the equilibria and Floquet multipliers of the cycles (numbers of the text)

"Jacobian of the 4D vector field f(x) = (dq, ddq) at the state `x`."
state_jacobian(x, vbp) = ForwardDiff.jacobian(
    y -> PM4.dynamics(PM4.build_u(y[1], y[2], y[3], y[4], zero(y[1])), vbp, 0)[1:4], collect(x))

"Monodromy matrix of the cycle `lc`: derivative of the time-T flow at its initial state."
function monodromy(lc, vbp; tol=DEFAULT_TOLERANCE)
    T = lc.t[end]
    flow = y -> PM4.integrate(PM4.build_u(y[1], y[2], y[3], y[4], zero(y[1])), T, vbp; tol).u[end][1:4]
    return ForwardDiff.jacobian(flow, collect(lc.u[1][1:4]))
end

fmt_complex(z) = abs(imag(z)) < 1e-8 ? f"{real(z):.2f}" :
    f"{real(z):.2f} \\pm {abs(imag(z)):.2f}\\,i"

"""
Spectra of the equilibria, grouped by σ-pairs (σ(α, τ) = (-α, τ + π) conjugates the
vector field, so both members of a pair have the same spectrum). The vector field is in
regular form (`PM4._accelerations_tension`, no division by the lever arm of the lines),
so its AD Jacobian is exact to rounding on the fold curves, where the equilibria lie:
both members of each pair are computed and required to agree.
Returns a vector of named tuples, one per pair, sorted by arm side then zone.
"""
function equilibrium_spectra(ss, vbp; is_down, is_center)
    spec = [eigvals(state_jacobian((q[1], q[2], 0.0, 0.0), vbp)) for q in ss]
    all(s -> all(isfinite, s), spec) || error("non-finite Jacobian at an equilibrium")
    pairs = []
    used = falses(length(ss))
    for i in eachindex(ss)
        used[i] && continue
        qi = ss[i]
        j = findfirst(k -> !used[k] && k != i &&
            abs(rem2pi(ss[k][1] + qi[1], RoundNearest)) < 1e-6 &&
            abs(rem2pi(ss[k][2] - qi[2] - π, RoundNearest)) < 1e-6, eachindex(ss))
        j === nothing && error("equilibrium $(rad2deg.(qi)) deg has no σ-twin")
        used[i] = used[j] = true
        d = maximum(abs, sort(spec[i]; by=reim) - sort(spec[j]; by=reim))
        d < 1e-6 || @warn "σ-twins with different spectra" rad2deg.(qi) d
        # representative with τ in [0, π), as in the figures' description
        r = rem2pi(ss[i][2], RoundDown) < π ? i : j
        λ = sort(spec[r]; by=z -> (real(z), imag(z)))
        # weight of the arm (α, dα) in the unstable eigenvectors: the text says they
        # move the kite along the eight, the arm staying at rest to first order
        E = eigen(state_jacobian((ss[r][1], ss[r][2], 0.0, 0.0), vbp))
        arm = maximum(findall(z -> real(z) > 0, E.values)) do k
            v = E.vectors[:, k]
            max(abs(v[1]), abs(v[3])) / maximum(abs, v)
        end
        push!(pairs, (; q=ss[r], down=is_down(ss[r]), center=is_center(ss[r]), λ,
            nunstable=count(z -> real(z) > 0, λ), remax=maximum(real, λ), arm))
    end
    sort!(pairs; by=p -> (!p.down, !p.center))
    return pairs
end

"""
Equilibria as roots of one scalar function per fold curve. At rest on the fold
α = φ(τ) + s (s = 0: downwind arm, s = π: upwind arm) the lever arm of the lines vanishes,
so α̈ = 0 identically, and the equilibria are the zeros of g_s(τ) = τ̈(φ(τ) + s, τ, 0, 0).
They are bracketed by the sign changes of g_s on a grid of `n` values of τ, then refined.
Returns the roots as [α, τ] (angles in (-π, π]).
"""
function fold_roots(vbp; n=3600)
    roots = Vector{Float64}[]
    for s in (0.0, π)
        q(τ) = (PM4.τ_to_θφ(τ, vbp)[2] + s, τ)
        g(τ) = PM4.dynamics(SA[q(τ)..., 0.0, 0.0, 0.0], vbp)[4]
        τs = range(0, 2π, n + 1)
        gs = g.(τs)
        for i in 1:n
            sign(gs[i]) != sign(gs[i+1]) || continue
            τ0 = find_zero(g, (τs[i], τs[i+1]))
            push!(roots, rem2pi.(collect(q(τ0)), RoundNearest))
        end
    end
    return roots
end

"""
Lower bound, over an `n`×`n` grid of the torus of configurations, of the ratio between
the tangential part (on the plane of the cone) of F_a + F_g at rest and its norm. Off the
fold curves, slack lines at rest require this tangential part to vanish (projections of
Newton's law for the kite on r̂ and r̂'), so a positive minimum excludes equilibria with
slack lines.
"""
function slack_ratio(vbp; n=721)
    minimum(Iterators.product(range(-π, π, n), range(0, 2π, n))) do (α, τ)
        k = PM4._kinematics(α, τ, vbp)
        F_a, F_g, _, _ = PM4._forces(k, SA[0.0, 0.0], vbp)
        F = F_a + F_g
        sqrt((F ⋅ k.rhat)^2 + (F ⋅ k.rhatp)^2 / k.nrhatp2) / norm(F)
    end
end

"""
Gap between the cycle `lc` and its image by σ shifted by the time `th` at which τ has
moved by ±π (max over the cycle, angles mod 2π); returns (; gap, half = th / T).
"""
function sigma_gap(lc; n=200)
    σ(u) = (-u[1], u[2] + π, -u[3], u[4])
    T, u0 = lc.t[end], lc.u[1]
    th = find_zero(t -> lc(t)[2] - (u0[2] + sign(u0[4]) * π), (0.0, T))
    gap = maximum(range(0, T, n)) do t
        y, z = σ(lc(t)), lc(mod(t + th, T))
        maximum(abs, (rem2pi(y[1] - z[1], RoundNearest), rem2pi(y[2] - z[2], RoundNearest),
            y[3] - z[3], y[4] - z[4]))
    end
    return (; gap, half=th / T)
end

function equilibria_table(pairs)
    rows = map(pairs) do p
        λs = String[]
        λ = p.λ
        k = 1
        while k <= length(λ)  # a complex pair once, with ±
            push!(λs, fmt_complex(λ[k]))
            k += abs(imag(λ[k])) < 1e-8 ? 1 : 2
        end
        arm, zone = p.down ? "downwind" : "upwind", p.center ? "centre" : "tip"
        α, τ = rad2deg(p.q[1]), rad2deg(rem2pi(p.q[2], RoundDown))
        ev = join(λs, raw",\ ")
        "    $arm & $zone & \$($(f"{α:.1f}"), $(f"{τ:.1f}"))\$ & \$$ev\$ & $(p.nunstable) \\\\"
    end
    return raw"""
\begin{table}[thpb]
\centering
\caption{Equilibria of the system for the parameters of Table~\ref{tab:params}, by pairs exchanged by the symmetry~\eqref{eq:symmetry}: one representative \((\alpha, \tau)\) of each pair (in degrees, \(0 \leq \tau < 180\)), eigenvalues of the Jacobian of \(f\) (in \si{\per\second}) and number of eigenvalues with positive real part.}
\label{tab:equilibria}
\begin{tabular}{l l c l c}
    \toprule
    {Arm} & {Kite} & {\((\alpha, \tau)\)} & {Eigenvalues} & {Unstable} \\
    \midrule
""" * join(rows, "\n") * raw"""

    \bottomrule
\end{tabular}
\end{table}"""
end

function run_analysis(; fast::Bool=false)
    ## Determine optimal parameters using kinetic model
    p0 = PMP.build_para()
    p0.b = PMP.optimal_torque_slope(p0; sense=(+), tol=1e-3)
    vbp = PMP.build_vbpara(p0)

    write_table("table1", param_table(p0))
    files = String[]

    #! format: off
    my_theme = (;
        downwind        = (; c=plot_color(:grey, 0.08)),
        upwind          = (; c=plot_color(:red, 0.08)),
        streamlines     = (; c=:black, alpha=0.7, lw=default(:lw) / 6),
        heads           = (; c=:black, alpha=0.5, lw=0),
        nullcline_tau   = (; c=1, ls=:solid),
        nullcline_alpha = (; c=2, ls=:solid),
        eq_legend       = (; c=:black, shape=:circle, ms=default(:ms), msw=20),
    )

    # Equilibrium characteristics: the symmetry (α, τ) -> (-α, τ+π) flips the kite
    # left↔right (φ = φ0 + Δφ·sin(τ)) and preserves the arm wind region and the
    # center/side zone, so the 8 equilibria realize all 2×2×2 combinations.
    φkite(τ) = vbp.φ0 + vbp.Δφ * sin(τ)                     # kite azimuth (PM4.τ_to_θφ)
    is_down(q) = cos(q[1]) > 0                              # arm downwind
    is_left(q) = φkite(q[2]) > vbp.φ0                       # kite left of window center
    is_center(q) = abs(φkite(q[2]) - vbp.φ0) < vbp.Δφ / 2   # closer to φ0 than to ±Δφ

    "Marker style from the equilibrium characteristics: fill = arm wind, shape = zone, stroke = kite side."
    ss_style(q) = (
        c = is_down(q) ? :white : :black,
        markerstrokecolor = is_down(q) ? :black : :white,
        shape = is_center(q) ? :circle : :diamond,
        ms = default(:ms),
        markerstrokewidth = is_down(q) ? 8 : 10,
        # markerstrokewidth = 10,
    )
    #! format: on

    theme = my_theme

    f = q -> SteadyState.ddq_partial(q, vbp)

    # Canvas & Layout
    δ = 0
    plot(xlims=(-180 - δ, 180 + δ), ylims=(-180 - δ, 180 + δ), aspect_ratio=:equal,
        xticks=([-180, 0, 180], ["-180", "0", "180"]),
        yticks=([-180, 0, 180], ["-180", "0", "180"]),
        xlabel="α (deg)", ylabel="τ (deg)",
        grid=false, framestyle=:semi, legend=:outerright, colorbar=false)
    # Background wind spans
    vspan!([-180, -90]; theme.upwind..., label="Upwind arm")
    vspan!([90, 180]; theme.upwind...)
    vspan!([-90, 90]; theme.downwind...) #, label="Downwind arm")

    # Streamlines & Arrowheads
    line_x, line_y, ax, ay, au, av = streamlines(f, (-π, π), (-π, π))
    plot!(rad2deg.(line_x), rad2deg.(line_y); theme.streamlines...)

    ax, ay = rad2deg.(ax), rad2deg.(ay)
    L, W = rad2deg(0.2), rad2deg(0.06)
    heads = [Shape([(s[1] + L / 2 * u, s[2] + L / 2 * v),
        (s[1] - L / 2 * u - W * v, s[2] - L / 2 * v + W * u),
        (s[1] - L / 2 * u + W * v, s[2] - L / 2 * v - W * u)])
             for (s, u, v) in zip(zip(ax, ay), au, av)]
    plot!(heads; theme.heads...)

    # Nullclines (splatting reuses line properties for both contour and the GR legend fix)
    N = 20
    αs = range(-π, π, N + 1)
    τs = range(-π, π, 2N + 1)
    dd = [f((α, τ)) for τ in τs, α in αs]

    contour!(rad2deg.(αs), rad2deg.(τs), last.(dd); levels=[0], theme.nullcline_tau...)
    contour!(rad2deg.(αs), rad2deg.(τs), first.(dd); levels=[0], theme.nullcline_alpha...)
    plot!([NaN], [NaN]; theme.nullcline_tau..., label="D²(τ) = 0")
    plot!([NaN], [NaN]; theme.nullcline_alpha..., label="D²(α) = 0")

    # Equilibria
    ss = SteadyState.all_steady_states(vbp)

    styles = ss_style.(ss)
    scatter!(rad2deg.(first.(ss)), rad2deg.(last.(ss));
        c=getproperty.(styles, :c),
        markerstrokecolor=getproperty.(styles, :markerstrokecolor),
        shape=getproperty.(styles, :shape),
        ms=getproperty.(styles, :ms),
        markerstrokewidth=getproperty.(styles, :markerstrokewidth))
    scatter!([NaN], [NaN]; theme.eq_legend..., label="Equilibrium")
    plot!(widen=1.06)
    plot!(leftmargin=6Plots.mm, bottommargin=5Plots.mm, rightmargin=5Plots.mm,
        topmargin=5Plots.mm, size=plot_size(1, 1.14))
    append!(files, save_figure(current(), "equilibria_phase_space"))

    ss_downwind = [_ss for _ss in ss if cos(_ss[1]) > 0]
    ss_upwind = [_ss for _ss in ss if cos(_ss[1]) <= 0]
    p_downwind = VIZ.plot_eight_circle(ss_downwind, vbp; plot_kwargs=(; color=:black, lw=2),
        marker_kwargs=ss_style.(ss_downwind), project_eight=true)
    plot!(title="Downwind arm", background_color_subplot=theme.downwind.c)
    p_upwind = VIZ.plot_eight_circle(ss_upwind, vbp; plot_kwargs=(; color=:black, lw=2),
        marker_kwargs=ss_style.(ss_upwind), project_eight=true)
    plot!(title="Upwind arm", background_color_subplot=theme.upwind.c)

    plot(p_downwind, p_upwind; layout=(1, 2), size=plot_size(16 / 15, 2))
    append!(files, save_figure(current(), "equilibria_eight_circle"))

    u0 = PM4.build_u(; dτ=1)
    sol = PM4.integrate(u0, 10, vbp; save_everystep=true)

    sim_4D_fig = plot_time_series(sol)
    append!(files, save_figure(sim_4D_fig, "time_series_4D"))

    (p_conv_p, slope_p) = plot_convergence_panel(vbp, +1, 1; ntrace=0, legend=false)
    (p_conv_m, slope_m) = plot_convergence_panel(vbp, -1, 2; ntrace=0)
    plot!(p_conv_p, leftmargin=11Plots.mm)
    plot!(p_conv_m, leftmargin=5Plots.mm)
    plot!(p_conv_p, ylabel="Distance to limit")
    conv_two = plot(p_conv_p, p_conv_m; layout=(1, 2), size=SIZE_FULL)
    plot!(bottommargin=7Plots.mm)
    append!(files, save_figure(conv_two, "limit_cycle_convergence"))

    r_vec = (t, sol) -> PM4.compute_OK(PM4.compute_Rτ(sol(t), sol.prob.p), sol.prob.p)
    lc_p = LC.compute_limit_cycle(vbp; sense=(+), save_everystep=true)
    lc_m = LC.compute_limit_cycle(vbp; sense=(-), save_everystep=true)
    function compute_lens(lc_p, lc_m, τ, δ_factor; idxs=2:3)
        # τ(t) is monotone along a cycle and spans TAU0 ± 2π: find the time where
        # the τ state hits τ mod 2π, mapped into the cycle's τ range.
        function phase_time(sol, τ)
            t0, t1 = sol.prob.tspan
            lo, hi = minmax(sol(t0)[2], sol(t1)[2])
            τ_target = clamp(τ + 2π * ceil((lo - τ) / 2π), lo, hi)
            return find_zero(t -> sol(t)[2] - τ_target, (t0, t1))
        end
        yz(sol, t) = r_vec(t, sol)[idxs]
        r_p = yz(lc_p, phase_time(lc_p, τ))
        r_m = yz(lc_m, phase_time(lc_m, τ))
        c = (r_p + r_m) / 2
        δ = δ_factor * norm(r_p - r_m)
        return c, δ
    end

    function yz_split(sol, n_splits)
        t = range(sol.prob.tspan..., 300)
        yz = [r_vec(t, sol)[2:3] for t in t]
        Δt = (t[end] - t[1]) / (n_splits + 1)
        for i in 1:n_splits
            t_split = (i - 0.2) * Δt
            ind = searchsortedfirst(t, t_split)
            insert!(yz, ind + i - 1, [NaN, NaN])
        end
        return yz
    end
    n_splits = 4
    # arrows = [permutedims(repeat([true], n_splits - 1)) false]
    plot(aspect_ratio=:equal, xlabel="y (m)", ylabel="z (m)", xtick=-40:10:40, ytick=10:5:20)
    plot!(invert(yz_split(lc_m, n_splits))..., lw=4, c=2, arrows=true)
    plot!(invert(yz_split(lc_p, n_splits))..., lw=2, c=1, arrows=true)
    plot!([NaN], [NaN], c=1, label="D(τ) > 0")
    plot!([NaN], [NaN], c=2, label="D(τ) < 0")
    plot!([0], [26])
    c, δ = compute_lens(lc_p, lc_m, 2.6, 3)  # lens focus phase τ was set by hand
    lens!([c[1] - δ, c[1] + δ], [c[2] - δ, c[2] + δ], inset=(1, bbox(0.3, 0.05, 0.3, 0.4)), lw=1.5, lc=:grey)
    l = plot!().inset_subplots[1]
    xticks_ = (c[1]-δ):δ:(c[1]+δ)
    yticks_ = (c[2]-δ):δ:(c[2]+δ)
    xticks!(l, xticks_, [f"{t:.1f}" for t in xticks_])
    yticks!(l, yticks_, [f"{t:.1f}" for t in yticks_])
    plot!(l, tickfontsize=TICKFONTSIZE - 1)
    plot!(; leftmargin=5Plots.mm, bottommargin=5Plots.mm, size=plot_size(2 / 3, 1))
    append!(files, save_figure(current(), "yz_limit_cycles"))

    ## Numbers of the text: equilibria, cycles, inventory
    pairs = equilibrium_spectra(ss, vbp; is_down, is_center)
    push!(files, write_table("table5", equilibria_table(pairs)))
    length(ss) == 8 || @warn "the text describes 8 equilibria" n = length(ss)
    all(p -> p.nunstable > 0, pairs) || @warn "the text says that all equilibria are unstable"
    all(p -> p.arm < 1e-2, pairs) || @warn "the text says the unstable modes of the equilibria are along the eight" getproperty.(pairs, :arm)
    pick(down, center) = only(filter(p -> p.down == down && p.center == center, pairs))
    all(c -> pick(false, c).nunstable == pick(true, c).nunstable + 1, (true, false)) ||
        @warn "the text says that an upwind equilibrium has one more unstable eigenvalue than its downwind counterpart"
    cyc = map((lc_p, lc_m)) do lc
        μ = eigvals(monodromy(lc, vbp))
        itriv = argmin(abs.(μ .- 1))
        μnt = sort(μ[setdiff(eachindex(μ), itriv)]; by=abs, rev=true)  # non-trivial multipliers
        μn = abs(μnt[1])                                     # largest non-trivial modulus
        (; T=lc.t[end], P=lc.u[end][end] / lc.t[end], μn, μerr=abs(μ[itriv] - 1),
            μ1real=abs(imag(μnt[1])) < 1e-12, μpair=abs(μnt[2]),  # dominant real; then a complex pair
            αmax=rad2deg(maximum(u -> abs(u[1]), lc.u)), dτ0=lc.u[1][4])
    end
    cyc[1].T > cyc[2].T && cyc[1].P > cyc[2].P ||
        @warn "the text says the cycle with dτ > 0 has the longer period and the larger power" cyc
    for (c, s) in zip(cyc, (slope_p, slope_m))
        abs(-log10(c.μn) - s) < 0.2 || @warn "convergence rate not close to the Floquet multiplier" s -log10(c.μn)
        c.αmax < 90 || @warn "the text says |α| stays below 90 deg" c.αmax
        c.μ1real && c.μpair < c.μn ||
            @warn "the text says the dominant multiplier is real, followed by a complex pair" c
        s > -log10(c.μn) || @warn "the text says the measured rate is above the Floquet rate" s -log10(c.μn)
    end
    # tension at rest, exact on the fold curves (regular form)
    ten = [(; Fl=PM4.compute_line_tension(SA[q[1], q[2], 0.0, 0.0], vbp)) for q in ss]
    roots = fold_roots(vbp; n=fast ? 720 : 3600)
    length(roots) == length(ss) && all(q -> any(r -> maximum(abs, rem2pi.(r - q, RoundNearest)) < 1e-6, roots), ss) ||
        @warn "the text says the equilibria are the roots of g± on the fold curves" roots ss
    slack = slack_ratio(vbp; n=fast ? 181 : 721)
    slack > 0.1 || @warn "the text says no equilibrium has slack lines (tangential force bounded below)" slack
    all(t -> t.Fl > 0, ten) || @warn "the text says the lines are taut at the equilibria"
    sig = sigma_gap.((lc_p, lc_m))
    all(g -> g.gap < 1e-7 && abs(g.half - 0.5) < 1e-7, sig) ||
        @warn "the text says each cycle is invariant under σ with a half-period shift" sig
    Fl_cyc = [[PM4.compute_line_tension(lc(t), vbp) for t in range(lc.t[1], lc.t[end], 400)] for lc in (lc_p, lc_m)]
    all(F -> minimum(F) > 0, Fl_cyc) || @warn "the text says the lines stay taut along the cycles"
    n_inv, vmax_inv = 200, 10  # α sampled in [-π, π]: σ maps Σ onto {τ = π}, not onto Σ
    lcs = LC.all_limit_cycles(vbp; αmin=-π, αmax=π, N=n_inv, vmax=vmax_inv)
    length(lcs) == 2 || @warn "the text says the inventory finds the two cycles only" n = length(lcs)
    sci(x) = (e = floor(Int, log10(x)); f"{x / exp10(e):.1f} \\times 10^{{{e}}}")
    expo(x) = "10^{$(ceil(Int, log10(x)))}"
    # numbers quoted in the text of the analysis section (`tables/numbers.tex` belongs
    # to the optimization section); the rates with the precision of the figure's annotations
    push!(files, write_table("numbers_analysis", join([
        "% Generated by main.jl (analysis section): numbers quoted in the text.",
        "\\newcommand{\\ConvRatePlus}{$(f"{slope_p:.2f}")}",
        "\\newcommand{\\ConvRateMinus}{$(f"{slope_m:.2f}")}",
        "\\newcommand{\\ConvTol}{10^{$(round(Int, log10(LC_CONV_TOL)))}}",
        "\\newcommand{\\IntTol}{10^{$(round(Int, log10(DEFAULT_TOLERANCE)))}}",
        "\\newcommand{\\EqCount}{$(length(ss))}",
        "\\newcommand{\\EqReCenterMin}{$(f"{min(pick(true, true).remax, pick(false, true).remax):.1f}")}",
        "\\newcommand{\\EqReCenterMax}{$(f"{max(pick(true, true).remax, pick(false, true).remax):.1f}")}",
        "\\newcommand{\\EqReTipMin}{$(f"{min(pick(true, false).remax, pick(false, false).remax):.2f}")}",
        "\\newcommand{\\EqReTipMax}{$(f"{max(pick(true, false).remax, pick(false, false).remax):.2f}")}",
        "\\newcommand{\\LcPeriodPlus}{$(f"{cyc[1].T:.2f}")}",
        "\\newcommand{\\LcPeriodMinus}{$(f"{cyc[2].T:.2f}")}",
        "\\newcommand{\\LcPowerPlus}{$(f"{cyc[1].P / 1000:.2f}")}",
        "\\newcommand{\\LcPowerMinus}{$(f"{cyc[2].P / 1000:.2f}")}",
        "\\newcommand{\\LcAlphaMaxPlus}{$(f"{cyc[1].αmax:.0f}")}",
        "\\newcommand{\\LcAlphaMaxMinus}{$(f"{cyc[2].αmax:.0f}")}",
        "\\newcommand{\\LcDtauPlus}{$(f"{cyc[1].dτ0:.2f}")}",
        "\\newcommand{\\LcDtauMinus}{$(f"{cyc[2].dτ0:.2f}")}",
        "\\newcommand{\\FloquetPlus}{$(sci(cyc[1].μn))}",
        "\\newcommand{\\FloquetMinus}{$(sci(cyc[2].μn))}",
        "\\newcommand{\\FloquetRatePlus}{$(f"{-log10(cyc[1].μn):.2f}")}",
        "\\newcommand{\\FloquetRateMinus}{$(f"{-log10(cyc[2].μn):.2f}")}",
        "\\newcommand{\\FloquetPairRatePlus}{$(f"{-log10(cyc[1].μpair):.2f}")}",
        "\\newcommand{\\FloquetPairRateMinus}{$(f"{-log10(cyc[2].μpair):.2f}")}",
        "\\newcommand{\\FloquetTrivialErr}{$(expo(max(cyc[1].μerr, cyc[2].μerr)))}",
        "\\newcommand{\\EqTensionMin}{$(f"{minimum(t -> t.Fl, ten):.0f}")}",
        "\\newcommand{\\SlackRatio}{$(f"{floor(slack, digits=2):.2f}")}",
        "\\newcommand{\\FoldGrid}{$(fast ? 720 : 3600)}",
        "\\newcommand{\\SlackGrid}{$(fast ? 181 : 721)}",
        "\\newcommand{\\EqTensionMax}{$(f"{maximum(t -> t.Fl, ten):.0f}")}",
        "\\newcommand{\\SigmaGap}{$(expo(maximum(g -> g.gap, sig)))}",
        "\\newcommand{\\LcTensionMinPlus}{$(f"{minimum(Fl_cyc[1]) / 1000:.2f}")}",
        "\\newcommand{\\LcTensionMaxPlus}{$(f"{maximum(Fl_cyc[1]) / 1000:.2f}")}",
        "\\newcommand{\\LcTensionMinMinus}{$(f"{minimum(Fl_cyc[2]) / 1000:.2f}")}",
        "\\newcommand{\\LcTensionMaxMinus}{$(f"{maximum(Fl_cyc[2]) / 1000:.2f}")}",
        "\\newcommand{\\InvSeeds}{$n_inv}",
        "\\newcommand{\\InvVmax}{$vmax_inv}",
        "\\newcommand{\\InvCount}{$(length(lcs))}",
    ], "\n")))
    return files
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_analysis()
end
