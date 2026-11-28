# Section 3 figures and Tables 2-4: optimization of the design parameters
# Θ = (ℓ, I, b, Δφ) (tether length, arm inertia, braking coefficient, lemniscate
# half-width) at the reference wind speed, then at each wind speed, and the
# sensitivity to the lemniscate half-height Δθ (held fixed in the optimization).
#
# The reference is the default parameters with the torque slope set by the
# kinetic criterion b = I ω (`PMP.optimal_torque_slope`, one limit cycle per
# iteration), not the default torque slope. The optimizations are limited by an
# iteration count (not by wall-clock time), so the results do not depend on the
# machine; `fast=true` shrinks the limits and the grids to test the workflow.
#
# Outputs: figs/{optimization_wind_sweep_{power,tension,params},optimization_dtheta,
# optimization_landscape,optimization_timeseries_all}.{pdf,png} and
# tables/table{2,3,4}.tex.
#
# Not part of the paper build (kept in applications/jtcam/Section_Optimization.jl):
# the 10×10 shape map and the console report.
using KEEP: PointMass4 as PM4
using KEEP: PointMassPara as PMP
using KEEP: LimitCycle as LC
using KEEP: Optimization as OPT

using Plots
using PyFormattedStrings: @f_str
using SplitApplyCombine: invert
using Logging: with_logger, NullLogger
using LinearAlgebra: norm, eigvals
import ForwardDiff

include(joinpath(@__DIR__, "helpers.jl"))


## Settings
const TOL = 1e-5
# Ipopt iteration limits (a converged solve takes a few dozen iterations); the
# wall-clock time is only a safety net that is not expected to bind.
const MAX_ITER = 300
const SWEEP_MAX_ITER = 150
const DTHETA_MAX_ITER = 300   # per continuation step of `dtheta_opt_sweep` (warm start)
const DTHETA_STEP = 0.5       # largest continuation step in Δθ (deg)
const DTHETA_MIN_STEP = 0.05  # a step that fails is halved down to this (deg)
const WALL_SAFETY = 3600.0
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
const DTHETA_RANGE = (2.0, 20.0)              # Δθ sweep, in degrees
const GLIDE_LOW = 5                           # glide ratio of the sensitivity remark (soft kites: 4-6)
const PTO_TS = 1000.0                         # max wing loading T_s (N/m²) of the PTO estimate

# Marks shared by the wind sweep, the landscape and the Δθ plots: the floor of
# the landscape's y axis, the translucent band that marks the reference
# condition (the same as at the reference wind speed) and the thin line of the
# horizontal floor. Built here and passed to the plot functions.
const GUIDES = (; ymin=0.9, line=(; c=:black, lw=5, alpha=0.3),
    floor=(; c=:black, lw=1, alpha=1))

# Parameter names, units and scales from the stored value to the displayed one.
# An angle is rad2deg; a per-angle unit is the inverse, so that b·dα stays a
# torque when dα is read in deg/s.
const TEX = Dict(:l => (raw"$\ell$", raw"\si{\meter}", identity),
    :I => (raw"$I$", raw"\si{\kilogram\meter\squared}", identity),
    :b => (raw"$b$", raw"\si{\newton\meter\second\per\degree}", x -> x / rad2deg(1)),
    :Δφ => (raw"$\Delta\varphi$", raw"\si{\degree}", rad2deg))

# Parameter names and units for the axis titles and the legends
const LABELS = Dict(:l => "ℓ", :I => "I", :b => "b",
    :Δθ => "Δθ", :Δφ => "Δφ")
const UNITS = Dict(:l => "m", :I => "kg m²", :b => "N m deg s⁻¹",
    :Δθ => "deg", :Δφ => "deg")

## Helpers
## Helpers
"""
Average power, line-tension extrema, peak kite acceleration and period of a cycle."""
function cycle_stats(sol; n=2001)
    # `compute_limit_cycle` returns `nothing` when no periodic orbit was found: report
    # NaN instead of crashing, so the caller treats the point as a failure.
    sol === nothing && return (; P=NaN, Tmin=NaN, Tmax=NaN, amax=NaN, period=NaN)
    ts = range(sol.t[1], sol.t[end], n)
    us = sol.(ts)
    p = sol.prob.p
    tension = PM4.compute_line_tension.(us, Ref(p))
    # the tension is exact everywhere (regular form, PM4._accelerations_tension)
    all(isfinite, tension) || error("cycle_stats: non-finite line tension along the cycle")
    acceleration = norm.(PM4.kite_acceleration.(us, Ref(p)))
    return (; P=sol.u[end][5] / sol.t[end],
        Tmin=minimum(tension), Tmax=maximum(tension),
        amax=maximum(acceleration), period=sol.t[end])
end

# Limit cycle, optimization, and limit cycle of the optimized system, for one
# value of one parameter. `lower`/`upper` are the box constraints of `syms`.
# `p_fixed` is the parameter set the no-optimization cycle is taken at, so a
# continuation can start from a neighbour's optimum `p0` without moving the
# reference curve. `guess` is Ipopt's initial guess (see `clamp_guess`), which
# lets a retry start from the neighbouring cycle states instead of cold.
function compute_results(param, value; syms, p0, p_fixed=p0, lower, upper, sense, tol,
    iters=MAX_ITER, guess=(;))
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
        OPT.optimize(p0, syms, lower, upper; sense, tol, max_iter=iters,
            max_wall_time=WALL_SAFETY, initial_guess=guess)
    end
    p_opt = copy(p0)
    p_opt[syms] .= solution.params[syms]  # physical units

    # Limit cycle at the optimized parameters. A start that drifts to a corner of the
    # box can return parameters whose cycle does not close (CI run 45900062200, 7 m/s:
    # period 8 s): report the start as failed instead of stopping the whole sweep.
    vbp = PMP.build_vbpara(p_opt)
    lc = try
        LC.compute_limit_cycle(vbp; sense=sense, save_everystep=true, tol=tol)
    catch e
        @warn "no limit cycle at the optimized parameters: start discarded" param value status = stats.status exception = e
        return nothing
    end
    stats_opt = cycle_stats(lc)

    return stats_no_opt, stats_opt, collect(p_opt[syms]), stats, solution
end

"""
Wind sweep: hold `p_ref` fixed, or re-optimize `syms` inside `lower`/`upper` at
each wind speed, marching outward from the reference wind `i_ref`.

Keeping `p_ref` is a feasible point of the optimization at every wind speed, so
a re-optimization that converges *below* the fixed curve is not the maximum
there: it is a different branch. Each wind speed therefore tries its starts in
order (the reference, then the neighbouring wind's optimum with its cycle states
as Ipopt's initial guess), keeps the best of them, and stops as soon as one
beats the fixed curve; if none does, the point is kept and a warning is emitted.
The wide boxes make the objective multimodal, so no single start succeeds at
every wind speed.

Returns the statistics with and without re-optimization and the optimized
parameters.
"""
function wind_sweep(syms, p_ref, winds, lower, upper, i_ref; rtol=1e-3, iters=SWEEP_MAX_ITER, strict=true)
    n = length(winds)
    no_opt = Vector{Any}(undef, n)
    opt = Vector{Any}(undef, n)
    theta = Vector{Vector{Float64}}(undef, n)
    sols = Vector{Any}(undef, n)
    ok = falses(n)
    from_theta(θ) = (p=copy(p_ref); p[syms].=θ; p)
    guess_of(j) = isassigned(sols, j) ? clamp_guess(sols[j]) : (;)

    "Starts at wind speed `i`: the reference, then the already-solved neighbour."
    function starts(i)
        i == i_ref && return ((p_ref, (;)),)
        j = i < i_ref ? i + 1 : i - 1
        isassigned(theta, j) || return ((p_ref, (;)),)
        return ((p_ref, (;)), (from_theta(theta[j]), guess_of(j)), (p_ref, guess_of(j)))
    end

    function solve!(i)
        best = nothing
        for (p0, guess) in starts(i)
            r = compute_results(:w0, winds[i];
                syms, p0, p_fixed=p_ref, lower, upper, sense=SENSE, tol=TOL,
                iters=iters, guess=guess)
            r === nothing && continue  # no cycle at the optimized parameters
            (!strict || r[4].status == :first_order) || continue  # strict=false: smoke tests keep unconverged points
            (best === nothing || r[2].P > best[2].P) && (best = r)
            # Keeping `p_ref` is always feasible, so the optimum cannot be below the
            # fixed curve: a converged point below it is another branch, try the next start.
            r[2].P >= r[1].P * (1 - rtol) && break
        end
        best === nothing && return (ok[i] = false)
        no_opt[i], opt[i], theta[i], sols[i] = best[1], best[2], best[3], best[5]
        ok[i] = opt[i].P >= no_opt[i].P * (1 - rtol)
        return nothing
    end

    for i in i_ref:-1:1
        @info "  wind sweep $(i)/$(n): $(winds[i]) m/s"
        solve!(i)
    end
    for i in (i_ref+1):n
        @info "  wind sweep $(i)/$(n): $(winds[i]) m/s"
        solve!(i)
    end
    for i in eachindex(winds)
        ok[i] && continue
        isassigned(opt, i) || (@warn "wind sweep at $(winds[i]) m/s: no converged re-optimization"; continue)
        @warn "wind sweep at $(winds[i]) m/s: best re-optimization ($(round(opt[i].P / 1000; digits=3)) kW) does not beat the fixed parameters ($(round(no_opt[i].P / 1000; digits=3)) kW)"
    end
    return (; syms, no_opt, opt, theta)
end

maxnormalize(x) = x ./ maximum(x)

"""
Neighbour solution as an Ipopt initial guess. The cycle states are clamped just
inside the box `OPT.optimize` asserts on: a cycle that starts at α = ±π puts the
guess on the bound, which `optimize` rejects.
"""
clamp_guess(::Nothing) = (;)
clamp_guess(s) = (; α0=clamp(s.α0, -π + 1e-6, π - 1e-6),
    dα0=clamp(s.dα0, -100 + 1e-6, 100 - 1e-6),
    dτ0=clamp(s.dτ0, 1e-6, 100 - 1e-6),
    tf=clamp(s.tf, 1.1e-3, 99.0), P=clamp(s.P, 1.1e-3, 9.9e4))

"""
One optimization step, limited to `iters` Ipopt iterations.
"""
function opt_step(p, syms, lower, upper; iters=MAX_ITER, strict=false)
    solution, stats, _ = with_logger(NullLogger()) do
        OPT.optimize(p, syms, lower, upper; sense=SENSE, tol=TOL, max_iter=iters,
            max_wall_time=WALL_SAFETY)
    end
    if stats.status != :first_order
        msg = "opt_step($(join(syms, ", "))) did not converge: $(stats.status) after $(stats.iter) iterations"
        strict ? error(msg) : @warn msg
    end
    p_next = copy(p)
    p_next[syms] .= solution.params[syms]  # physical units
    return (; p=p_next, P=solution.P, iter=stats.iter)
end

"""
Δθ sweep at the optimum: take the limit cycle at each absolute angle (in
degrees) in `angles_deg` with Θ fixed. Starting from the same
`u0 = (α, τ, dα, dτ) = (0, τ0, 0, 1)` picks the same (dτ > 0) pumping cycle at
every point.
"""
function dtheta_study(p, angles_deg; n=2001)
    u0 = PM4.build_u(; dτ=1)
    return map(angles_deg) do a
        q = copy(p)
        q.Δθ = deg2rad(a)
        sol = LC.compute_limit_cycle(u0, PMP.build_vbpara(q);
            tol=TOL, save_everystep=true, tf=100.0)
        (; Δθ=a, cycle_stats(sol; n)...)
    end
end

"""
Δθ sweep with re-optimization, by continuation: march outward from the nominal
angle `nominal` (degrees, where `p_ref` is the optimum) to each angle of
`angles_deg`, in steps of at most `step` degrees, and re-optimize `syms` inside
`lower`/`upper` at every step from the previous optimum, with its cycle states as
Ipopt's initial guess. Cold starts away from the nominal drift to the corners of
the loose `Θ` box (I, b up to 1e6) and fail; small warm-started steps stay on
the nominal branch. A failed step is retried with half the step, down to
`min_step`; below it, the remaining angles on that side are reported as `NaN`,
not as a gain, and the angles where it is lost are returned (`NaN` if it is not).
`iters` is the per-step iteration limit. Holding `Θ` at the
nominal is always feasible, so a converged point below the fixed-`Θ` value is
flagged (another branch). Statistics come from `dtheta_study` at the solved
point, i.e. from the same dτ > 0 cycle as the fixed curve.
"""
function dtheta_opt_sweep(p_ref, syms, angles_deg, lower, upper; nominal,
    step=DTHETA_STEP, min_step=DTHETA_MIN_STEP, iters=DTHETA_MAX_ITER, rtol=1e-3)
    function attempt(a, p0, g)
        q = copy(p0)
        q.Δθ = deg2rad(a)
        try
            solution, st, _ = with_logger(NullLogger()) do
                OPT.optimize(q, syms, lower, upper; sense=SENSE, tol=TOL,
                    max_iter=iters, max_wall_time=WALL_SAFETY, initial_guess=g)
            end
            p = copy(q)
            p[syms] .= solution.params[syms]
            return (p=p, ok=st.status == :first_order, sol=solution, P=solution.P, iter=st.iter)
        catch
            return (p=copy(p0), ok=false, sol=nothing, P=NaN, iter=0)
        end
    end
    failed(a) = (; Δθ=a, P=NaN, Tmin=NaN, Tmax=NaN, amax=NaN, period=NaN)
    function record(a, cur)
        fixed = only(dtheta_study(p_ref, [a]))   # Θ held at the nominal
        cur.P < fixed.P * (1 - rtol) &&
            @warn "re-optimization at Δθ = $(a)° is below the fixed Θ ($(round(fixed.P / 1000; digits=3)) kW): another branch" cur.P
        return only(dtheta_study(cur.p, [a]))
    end
    # the nominal point itself: p_ref is the optimum there, solve once for its cycle states
    start = attempt(nominal, p_ref, (;))
    start.ok || @warn "Δθ continuation: the nominal re-optimization at $(nominal)° does not converge"
    stats = Vector{Any}(undef, length(angles_deg))
    lost_at = Dict(-1 => NaN, +1 => NaN)  # where the continuation is lost, per side
    for side in (-1, +1)
        # the nominal angle itself (zero steps) goes with the right side
        idx = findall(a -> side < 0 ? a < nominal - 1e-9 : a >= nominal - 1e-9, angles_deg)
        sort!(idx; by=i -> abs(angles_deg[i] - nominal))
        cur, a, h, lost = start, nominal, step, !start.ok
        for i in idx
            target = angles_deg[i]
            while !lost && abs(target - a) > 1e-9
                next = a + side * min(h, abs(target - a))
                r = attempt(next, cur.p, clamp_guess(cur.sol))
                if r.ok
                    cur, a, h = r, next, min(step, 2h)
                else
                    h /= 2
                    h < min_step && (lost = true; lost_at[side] = next;
                        @warn "Δθ continuation lost at $(round(next; digits=3))° (from $(a)°)")
                end
            end
            stats[i] = lost ? failed(target) : try
                record(target, cur)
            catch e  # no closed cycle at the solved point: a gap, not a crash
                @warn "Δθ = $(target)°: no limit cycle at the re-optimized parameters" exception = e
                failed(target)
            end
            @info "  Δθ = $(target)°: " * (lost ? "failed" : "P = $(round(stats[i].P / 1000; digits=3)) kW")
        end
    end
    return (; stats, lost_below=lost_at[-1], lost_above=lost_at[+1])
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
        lc = try
            LC.compute_limit_cycle(PMP.build_vbpara(q); sense=SENSE, tol=tol)
        catch
            nothing        # no closed cycle with this sign: report a gap, not a crash
        end
        return lc === nothing ? NaN : lc.u[end][5] / lc.t[end]
    end
    vals = [power(sym, x) for sym in syms, x in xs]
    k0 = findmin(abs.(xs .- 1))[2]  # the optimum, where the landscape is normalized
    return (; xs=collect(xs), syms=collect(syms), window,
        ys=[row ./ row[k0] for row in eachrow(vals)],
        center=collect(p[syms]),  # the optimum in physical units, for the x axis
        box=[(lo / p[sym], up / p[sym]) for (sym, lo, up) in zip(syms, lower, upper)])
end

## Tables
"One optimized parameter, in the units shown in the table."
struct OptVar
    sym::String
    unit::String
    low::Float64
    up::Float64
    guess::Float64
    sol::Float64
end

"Limit-cycle statistics of one parameter set, plus the norm of its endpoint residual."
function cycle_report(p)
    sol = LC.compute_limit_cycle(PMP.build_vbpara(p); sense=SENSE, tol=TOL,
        save_everystep=true)
    return (; cycle_stats(sol)..., residual=norm(LC.endpoint_residuals(sol; sense=SENSE)))
end

"Math-mode scientific notation, e.g. `\$3.7\\times10^{-9}\$`."
function latex_sci(x)
    x == 0 && return raw"$0$"
    e = floor(Int, log10(abs(x)))
    return raw"$" * string(round(x / exp10(e); digits=1)) * raw"\times10^{" *
        string(e) * raw"}$"
end

"""
The two siunitx tables of the optimization section: the boxes of Θ = (ℓ, I, b,
Δφ) and the optimization results at `WIND_REF`, initial guess against solution.
Every number comes from `p_init`, `p_opt` and their limit cycles; Δφ is shown in
degrees and b per degree. Δθ is held fixed, so it only enters the first caption.
"""
function latex_tables(p_init, p_opt, lower, upper)
    s0, s1 = cycle_report(p_init), cycle_report(p_opt)
    data = [
        let (sym, unit, scale) = TEX[k]
            OptVar(sym, unit, scale(lower[i]), scale(upper[i]),
                scale(p_init[k]), scale(p_opt[k]))
        end for (i, k) in enumerate(SYMS)
    ]

    # b per degree has a small lower bound (1 N m s/rad -> 0.0175)
    num(x) = abs(x) < 0.1 ? f"{x:.4f}" : f"{x:.1f}"
    bounds_rows = [f"    {v.sym} ({v.unit}) & {num(v.low)} & {num(v.up)} \\\\ " for v in data]
    results_rows = [f"    {v.sym} ({v.unit}) & {v.guess:.1f} & {v.sol:.1f} \\\\ " for v in data]
    push!(results_rows, raw"    \midrule")
    push!(results_rows, f"    Average power (\\si{{\\kilo\\watt}}) & {s0.P / 1000:.2f} & {s1.P / 1000:.2f} \\\\ ")
    push!(results_rows, f"    Period (\\si{{\\second}}) & {s0.period:.2f} & {s1.period:.2f} \\\\ ")
    push!(results_rows, f"    Min line tension (\\si{{\\kilo\\newton}}) & {s0.Tmin / 1000:.2f} & {s1.Tmin / 1000:.2f} \\\\ ")
    push!(results_rows, f"    Max line tension (\\si{{\\kilo\\newton}}) & {s0.Tmax / 1000:.2f} & {s1.Tmax / 1000:.2f} \\\\ ")
    push!(results_rows, "    Norm of residual & {" * latex_sci(s0.residual) * "} & {" *
        latex_sci(s1.residual) * "} \\\\ ")

    table_descr = join([
        raw"\begin{table}[thpb]",
        raw"\centering",
        raw"\caption{Description of the optimized parameters ($\Delta\theta = \SI{" *
            string(round(rad2deg(p_init.Δθ); digits=2)) * raw"}{\degree}$ fixed).}",
        raw"\begin{tabular}{l S S}",
        raw"    \toprule",
        "    {Symbol (Unit)} & {Lower bound} & {Upper bound} \\\\",
        raw"    \midrule",
        bounds_rows...,
        raw"    \bottomrule",
        raw"\end{tabular}",
        raw"\label{tab:opt_descr}",
        raw"\end{table}",
    ], "\n")

    table_opt = join([
        raw"\begin{table}[thpb]",
        raw"\centering",
        raw"\caption{Optimization results at $w_0 = \SI{" *
            string(round(Int, WIND_REF)) * raw"}{\meter\per\second}$.}",
        raw"\begin{tabular}{l S S}",
        raw"    \toprule",
        "    {Symbol (Unit)} & {Initial Guess} & {Solution} \\\\",
        raw"    \midrule",
        results_rows...,
        raw"    \bottomrule",
        raw"\end{tabular}",
        raw"\label{tab:opt_results}",
        raw"\end{table}",
    ], "\n")

    return table_descr, table_opt
end

"""
Curvature of the kite trajectory without arm (a = 0), for the eight of `p`: the
geodesic curvature κ_g of the spherical eight (section "Curvature of the kite
trajectory") and ℓκ = √(1 + κ_g²), at the lobe tip τ = π/2 and maximized over a
fine τ grid.
"""
function curvature_stats(p; n=200_000)
    Δθ, Δφ, θ0 = p.Δθ, p.Δφ, p.θ0
    function κg(τ)
        θ = θ0 + Δθ * sin(2τ)
        v = sqrt(4Δθ^2 * cos(2τ)^2 + Δφ^2 * sin(θ)^2 * cos(τ)^2)
        return (2Δθ * Δφ * sin(θ) * sin(τ) * (2 + cos(2τ)) +
                Δφ * cos(θ) * cos(τ) * (8Δθ^2 * cos(2τ)^2 + Δφ^2 * sin(θ)^2 * cos(τ)^2)) / v^3
    end
    K(τ) = sqrt(1 + κg(τ)^2)
    τs = range(0, 2π, n + 1)[1:end-1]
    Ks = K.(τs)
    i = argmax(Ks)
    # κ(τ + π) = κ(τ) (κ_g changes sign): report the maximizer modulo π, independent of
    # which of the two symmetric grid points argmax picks
    return (; kg_tip=κg(π / 2), K_tip=K(π / 2), K_max=Ks[i], τ_max=mod(τs[i], π))
end

"""
Floquet multipliers of the dτ > 0 cycle of `p`: the largest modulus among the three
non-trivial ones (`μ`) and the error on the trivial multiplier 1 (`err`), from the
monodromy matrix computed by forward-mode AD through the integrator.
"""
function floquet(p)
    vbp = PMP.build_vbpara(p)
    lc = LC.compute_limit_cycle(vbp; sense=SENSE, save_everystep=true)
    T = lc.t[end]
    flow = y -> PM4.integrate(PM4.build_u(y[1], y[2], y[3], y[4], zero(y[1])), T, vbp).u[end][1:4]
    μ = eigvals(ForwardDiff.jacobian(flow, collect(lc.u[1][1:4])))
    itriv = argmin(abs.(μ .- 1))
    return (; μ=maximum(abs, μ[setdiff(eachindex(μ), itriv)]), err=abs(μ[itriv] - 1))
end

"""
A posteriori curvature check on the optimized cycle of `p` (dτ > 0), with the arm: the
curvature κ = ‖ṙ × r̈‖/‖ṙ‖³ of the kite trajectory, maximized over the cycle; the bounds
A₁ = max|α'|, A₂ = max|α''| (prime = d/dτ) of the arm law along the cycle; and the upper
bound (eq:curv-bound) evaluated with the box of the optimization (ℓ ≥ ℓ₋ = L_BOUNDS[1],
Δφ in DPHI_BOUNDS, Δθ and θ0 fixed), and at the design point itself (ℓ₋ = ℓ, Δφ₋ = Δφ₊ =
Δφ), with these A₁, A₂. Returns radii (m): `rmin` the smallest radius of curvature of
the cycle, `rbox` and `rpt` the radii guaranteed by the bound (both ≤ `rmin`, checked).
"""
function curvature_check(p; n=4000)
    vbp = PMP.build_vbpara(p)
    lc = LC.compute_limit_cycle(vbp; sense=SENSE, tol=TOL, save_everystep=true)
    vel(u) = (k = PM4._kinematics(u[1], u[2], vbp);
              (k.R_α * u[3] + k.R_τ * u[4]) * k.rhat + k.R * u[4] * k.rhatp)
    cross3(x, y) = (x[2]y[3] - x[3]y[2], x[3]y[1] - x[1]y[3], x[1]y[2] - x[2]y[1])
    rows = map(range(0, lc.t[end], n + 1)[1:end-1]) do t
        u = collect(lc(t)[1:4])
        du = PM4.dynamics(u, vbp)[1:4]
        v = vel(u)
        acc = ForwardDiff.derivative(s -> vel(u .+ s .* du), 0.0)
        κ = norm(collect(cross3(v, acc))) / norm(v)^3
        α1 = u[3] / u[4]                                   # α' = α̇/τ̇
        α2 = (du[3] * u[4] - u[3] * du[4]) / u[4]^3          # α'' = (α̈τ̇ - α̇τ̈)/τ̇³
        (; κ, α1, α2)
    end
    κmax = maximum(r -> r.κ, rows)
    A1, A2 = maximum(r -> abs(r.α1), rows), maximum(r -> abs(r.α2), rows)
    # bound (eq:curv-bound) on a box: ℓ ≥ ℓm, Δφ in [Δφm, Δφp], Δθ and θ0 fixed
    function κbound(ℓm, Δφm, Δφp)
        a, Δθm, Δθp = p.a, p.Δθ, p.Δθ
        J = (p.θ0 - Δθp, p.θ0 + Δθp)
        0 < J[1] && J[2] < π || error("curvature bound: J not in (0, π)")
        js = range(J..., 10_001)
        m, M, C = minimum(sin, js), maximum(sin, js), maximum(x -> abs(cos(x)), js)
        vstar = m * Δφm >= 4Δθm ? 2Δθm : m * Δφm * sqrt(1 / 2 - m^2 * Δφm^2 / (64 * Δθm^2))
        κg = 2sqrt(2) * M * Δθp * Δφp / vstar^3 + 2C / m
        Γ1 = 2Δθp + Δφp + A1
        Γ2 = Γ1^2 + 4Δθp + Δφp + A2
        V1 = 4Δθp + Δφp + Δφp^2 / 2 + 4Δθp * Δφp
        εp = a / sqrt(ℓm^2 - a^2)
        U = εp * (Γ2 + εp^2 * Γ1^2) / vstar^2 + εp * Γ1 * V1 / vstar^3
        return sqrt(κg^2 + (1 + U)^2) / (ℓm - a)
    end
    κbox = κbound(first(L_BOUNDS), DPHI_BOUNDS...)   # the box of the optimization
    κpt = κbound(p.l, p.Δφ, p.Δφ)                     # the design point itself
    κmax <= min(κbox, κpt) || @warn "curvature: the a posteriori bound is violated" κmax κbox κpt
    return (; rmin=1 / κmax, rbox=1 / κbox, rpt=1 / κpt, A1, A2)
end

"""
Power take-off estimate of section "Estimation of power take-off" (Alain Nême's 1D
model) for the parameters `p` of Table 1, with the maximum wing loading `Ts` (N/m²),
which is not a model parameter: kite speed at the centre of the eight, length of the
eight, pulsation ω, then the optimal design under the amplitude constraint |α| ≤ π/2 and
`I ≥ Imin` (the lower bound of the optimization box): `I = Imin`, `b` on the constraint,
and the power per unit kite area (`pto`; `pto0` is its supremum as `Imin → 0`; `Iseq` the
inertia of the sequential argument "b = Iω, then smallest admissible I").
"""
function pto_estimate(p; Ts=PTO_TS, Imin=I_BOUNDS[1], n=100_000)
    G = p.C_L / p.C_D
    Vk = p.w0 * (p.l * cos(p.θ0) / p.z0)^(1 / p.n_wind) * sqrt((1 + G^2) * sin(p.θ0)^2 - 1)
    τs = range(0, 2π, n + 1)[1:end-1]
    L81 = 2π / n * sum(τ -> sqrt((2p.Δθ * cos(2τ))^2 +
                                 sin(p.θ0 + p.Δθ * sin(2τ))^2 * (p.Δφ * cos(τ))^2), τs)
    ω = 2π * Vk / (p.l * L81)
    CM = Ts * p.S * p.a * sin(p.θ0)
    # max of W = b CM² / (2(I²ω² + b²)) under the amplitude constraint
    # CM / (ω √(I²ω² + b²)) ≤ π/2 and I ≥ Imin: at fixed I the best b is Iω while the
    # constraint allows it (I ≥ √2 CM/(πω²)); below, the constraint is active and W
    # increases as I decreases, so the optimum is I = Imin with the constraint active
    I = Imin
    ε = π * ω^2 * I / (2CM)
    ε < 1 / sqrt(2) || error("PTO estimate: Imin above √2 CM/(πω²), the optimum is not on the constraint")
    b = 2CM / (π * ω) * sqrt(1 - ε^2)
    W = b * CM^2 / (2 * (I^2 * ω^2 + b^2))
    W ≈ π * CM * ω / 4 * sqrt(1 - ε^2) || error("PTO estimate: closed form")
    return (; I, b, ε, pto=W / p.S, pto0=π * CM * ω / 4 / p.S, Iseq=sqrt(2) * CM / (π * ω^2))
end

"""
`tables/numbers.tex`: every number the text of the optimization section quotes, as
LaTeX macros, so that the text cannot drift from the generated tables and figures.
"""
function numbers_tex(; p_start, p_ref, opt_iter, winds, sweep, dtheta, dtheta_opt, dtheta_lost)
    s0, s1 = cycle_report(p_start), cycle_report(p_ref)
    ok = [isassigned(sweep.opt, i) && isassigned(sweep.no_opt, i) && isassigned(sweep.theta, i)
          for i in eachindex(winds)]
    w = collect(winds)[ok]
    P(v) = getproperty.(v, :P)
    o, fx = sweep.opt[ok], sweep.no_opt[ok]
    loss = maximum(100 .* (P(o) .- P(fx)) ./ P(o))
    dT = maximum(100 .* (getproperty.(fx, :Tmax) ./ getproperty.(o, :Tmax) .- 1))
    Tmin = minimum(min.(getproperty.(o, :Tmin), getproperty.(fx, :Tmin)))
    # power-law exponents along the optimized sweep (least squares in log-log)
    slope(y) = ([ones(length(w)) log.(w)] \ log.(y))[2]
    b_opt = [θ[findfirst(==(:b), SYMS)] for θ in sweep.theta[ok]]
    c = curvature_stats(p_ref)
    # Δθ sweep (caption of fig:opt-dtheta): power relative to the nominal with Θ fixed
    ang = [d.Δθ for d in dtheta]
    i0 = argmin(abs.(ang .- rad2deg(p_ref.Δθ)))
    Pfix, Popt = P(dtheta) ./ dtheta[i0].P, P(dtheta_opt) ./ dtheta[i0].P
    # the caption claims: Θ fixed peaks at the nominal; re-optimized decreases
    # steadily above the nominal; the continuation is lost below the nominal and
    # converges up to the largest angle
    argmax(Pfix) == i0 || @warn "Δθ caption: with Θ fixed the power no longer peaks at the nominal angle"
    all(isfinite, Popt[i0:end]) && all(diff(Popt[i0:end]) .< 0) ||
        @warn "Δθ caption: the re-optimized power no longer decreases steadily above the nominal angle"
    # below the nominal angle, the continuation is fragile: depending on the host, it is
    # lost near 4° (CI runs 46139319516, 45912488604) or reaches 2° with a peak at 4°
    # (run 46137278622). The caption fragment `DthetaBelow` states what this run gives.
    ibest = argmax(i -> isfinite(Popt[i]) ? Popt[i] : -Inf, eachindex(Popt))  # an index
    below = if isfinite(dtheta_lost)
        "below the nominal angle, the continuation is lost at \\(" * f"{dtheta_lost:.1f}" *
        "^\\circ\\), so that no re-optimized value is shown at smaller angles"
    else
        all(isfinite, Popt) || @warn "Δθ caption: re-optimized values missing" Popt
        "the continuation reaches the whole range, and the re-optimized power peaks, among " *
        "the sampled angles, at \\(" * f"{ang[ibest]:.1f}" * "^\\circ\\), " *
        f"{100 * (Popt[ibest] - 1):.1f}" * "\\,\\% above its nominal value"
    end
    pto = pto_estimate(p_start)
    cc = curvature_check(p_ref)
    fq0, fq1 = floquet(p_start), floquet(p_ref)
    # sensitivity to the glide ratio G = C_L/C_D: the optimized design and the PTO
    # estimate with the lower G of soft kites (the value the parameters' comments cite)
    pG = copy(p_ref); pG.C_D = pG.C_L / GLIDE_LOW
    sG = cycle_report(pG)
    ptoG = pto_estimate(let q = copy(p_start); q.C_D = q.C_L / GLIDE_LOW; q end)
    isfinite(sG.P) || @warn "glide-ratio sensitivity: no cycle at G = $GLIDE_LOW" 
    fq1.μ < 1 && fq1.err < 1e-6 || @warn "the text says the optimized cycle is attracting" fq1
    fq1.μ > fq0.μ || @warn "the text says the optimized cycle is less strongly attracting than the initial one" fq0 fq1
    r1(x) = f"{x:.1f}"
    sciq(x) = (e = floor(Int, log10(x)); f"{x / exp10(e):.1f} \\times 10^{{{e}}}")  # math mode
    r2(x) = f"{x:.2f}"
    macros = [
        "WindRef" => string(round(Int, WIND_REF)),
        "IpoptTol" => "10^{$(round(Int, log10(TOL)))}",
        "IvpTol" => "10^{$(round(Int, log10(TOL / 10)))}",
        "PostTol" => "10^{$(round(Int, log10(TOL)))}",
        "IpoptIter" => string(opt_iter),
        "PowerInit" => r2(s0.P / 1000), "PowerOpt" => r2(s1.P / 1000),
        "PowerGain" => string(round(Int, 100 * (s1.P / s0.P - 1))),
        "PeriodInit" => r2(s0.period), "PeriodOpt" => r2(s1.period),
        "EllOpt" => r1(p_ref.l), "DphiOpt" => r1(rad2deg(p_ref.Δφ)),
        "DthetaFixed" => r1(rad2deg(p_ref.Δθ)), "ThetaZero" => r1(rad2deg(p_ref.θ0)),
        "SweepWindMin" => string(round(Int, first(w))), "SweepWindMax" => string(round(Int, last(w))),
        "SweepMaxLoss" => r1(loss), "SweepMaxTensionIncrease" => r1(dT),
        "SweepMinTension" => r2(Tmin / 1000),
        "ExpPower" => r1(slope(P(o))), "ExpTension" => r1(slope(getproperty.(o, :Tmax))),
        "ExpBrake" => r1(slope(b_opt)), "ExpPeriod" => r1(slope(getproperty.(o, :period))),
        "CurvTipKg" => r2(c.kg_tip), "CurvTipK" => r2(c.K_tip),
        "CurvMaxK" => r2(c.K_max), "CurvMaxTau" => f"{c.τ_max:.3f}",
        "CurvMaxRadius" => r1(p_ref.l / c.K_max),
        "GlideRef" => string(round(Int, p_ref.C_L / p_ref.C_D)), "GlideLow" => string(GLIDE_LOW),
        "PowerOptGlideLow" => r2(sG.P / 1000), "PtoPowerGlideLow" => r1(ptoG.pto / 1000),
        "FloquetInit" => sciq(fq0.μ), "FloquetOpt" => sciq(fq1.μ),
        "CurvArmRadius" => r2(cc.rmin), "CurvBoxRadius" => f"{1000 * cc.rbox:.1f}",
        "CurvPointRadius" => r2(cc.rpt),
        "CurvAOne" => r2(cc.A1), "CurvATwo" => r2(cc.A2),
        "PtoTs" => string(round(Int, PTO_TS)),
        "PtoI" => string(round(Int, pto.I)), "PtoB" => r1(deg2rad(pto.b)),  # b per degree, as in the tables
        "PtoPower" => r1(pto.pto / 1000), "PtoPowerZero" => r1(pto.pto0 / 1000),
        "PtoEps" => f"{pto.ε:.3f}", "PtoISeq" => string(round(Int, pto.Iseq)),
        "PtoGainSeq" => string(round(Int, 100 * (pto.pto / (pto.pto0 / sqrt(2)) - 1))),
        "DthetaStep" => string(DTHETA_STEP),
        "DthetaMax" => string(round(Int, last(ang))),
        "DthetaFixedLoss" => string(round(Int, 100 * (1 - last(Pfix)))),
        "DthetaBelow" => below,
        "DthetaOptLoss" => r1(100 * (1 - last(Popt))),
    ]
    return join(["% Generated by main.jl (optimization section): numbers quoted in the text.";
        ["\\newcommand{\\$k}{$v}" for (k, v) in macros]], "\n")
end

"""
Table 4: the initial guess (I = 1300, b = I ω: not the PTO optimum of the 1D estimate) against
the dynamic optimum, with the statistics of their limit cycles.
"""
function comparison_table(p_init, p_opt)
    stats(p) = cycle_stats(LC.compute_limit_cycle(PMP.build_vbpara(p); sense=SENSE, tol=TOL, save_everystep=true))
    s0, s1 = stats(p_init), stats(p_opt)
    g(x) = string(round(x; sigdigits=6))
    param_rows = [f"    {TEX[k][1]} ({TEX[k][2]}) & {g(TEX[k][3](p_init[k]))} & {g(TEX[k][3](p_opt[k]))} \\\\ " for k in SYMS]
    rows = [
        param_rows...,
        raw"    \midrule",
        f"    Period (\\si{{\\second}}) & {g(s0.period)} & {g(s1.period)} \\\\ ",
        "    Objective " * raw"$P$" * raw" (\si{\kilo\watt}) & " * g(s0.P / 1000) * " & " * g(s1.P / 1000) * raw" \\ ",
        f"    Min tension (\\si{{\\kilo\\newton}}) & {g(s0.Tmin / 1000)} & {g(s1.Tmin / 1000)} \\\\ ",
        f"    Max tension (\\si{{\\kilo\\newton}}) & {g(s0.Tmax / 1000)} & {g(s1.Tmax / 1000)} \\\\ ",
        f"    Max acceleration (\\si{{\\meter\\per\\square\\second}}) & {g(s0.amax)} & {g(s1.amax)} \\\\ ",
    ]
    return join([
        raw"% WARNING: values: to be checked vs. Alain's computation",
        raw"\begin{table}[htbp]",
        raw"  \centering",
        raw"  \caption{Comparison of initial and optimized parameters}",
        raw"  \label{tab:parameter_comparison}",
        raw"  \begin{tabular}{lrr}",
        raw"    \toprule",
        raw"    Parameter & \multicolumn{1}{c}{Initial guess} & \multicolumn{1}{c}{Dynamic optimization} \\ ",
        raw"    \midrule",
        rows...,
        raw"    \bottomrule",
        raw"  \end{tabular}",
        raw"\end{table}",
    ], "\n")
end


## Plots
## Plot: wind sweep over Θ
"""
Three half-width panels of the wind sweep: average power, line tension (min–max
envelope), and the re-optimized parameters normalized by their maximum over the
wind range. The dotted black curve is the 9 m/s reference held fixed, and a
translucent black band marks that reference wind speed in every panel. The
parameter panel has no legend: the braking coefficient is annotated at the lower
left, where its curve is, and the other three are stacked on the right in their
vertical order at the last wind speed.
"""
function plot_wind_sweep(winds, sweep, i_ref, guides)
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
    f1 = plot(ylabel="Power (kW)", legend=:topleft, xformatter=:plain, xlabel="Wind speed (m/s)", ylims=(0, 60), yticks=0:20:100, widen=true)
    plot!(f1, winds, power0 ./ 1000; ref_style...)
    plot!(f1, winds, power_vec ./ 1000; c=1)

    # Line tension, min–max envelope (lower is better for wearing the lines out)
    f2 = plot(ylabel="Line tension (kN)", legend=:topleft, xformatter=:plain, xlabel="Wind speed (m/s)", ylims=(0, 40), yticks=0:20:100, widen=true)
    plot!(f2, winds, [max_tension0 ./ 1000, min_tension0 ./ 1000]; ref_style...)
    plot!(f2, winds, max_tension_vec ./ 1000; fillrange=min_tension_vec ./ 1000, c=1, fillalpha=0.1)
    plot!(f2, winds, min_tension_vec ./ 1000; c=1)

    # Optimized parameters, normalized by their maximum over the wind range
    n = length(syms)
    f3 = plot(yticks=((0, 1), ("0", "max")), ylim=(0, 1), ylabel="Parameters", xformatter=:plain, xlabel="Wind speed (m/s)", legend=false, widen=true)
    # plot!(f3, winds, [[p for _ in winds] for p in invert(theta_norm)[i_ref]]; lw=1.5, ref_style...)
    plot!(f3, winds, theta_norm; c=permutedims(collect(1:n)))

    # the reference wind speed is marked by the same band as the other guide
    # lines; the backend draws series in list order, so moving it in front of
    # the data keeps the band under the curves
    for f in (f1, f2, f3)
        vline!(f, [WIND_REF]; guides.line...)
        pushfirst!(f[1].series_list, pop!(f[1].series_list))
    end

    # annotations instead of a legend: the braking coefficient at the lower
    # left (its curve is the lowest there), then the other three stacked on the
    # right below the curves, in their vertical order at the last wind speed
    names = Dict(:l => "Line length", :I => "Arm inertia",
        :b => "Arm braking coefficient", :Δφ => "Lemniscate width")
    ib = findfirst(==(:b), syms)
    annotate!(f3, first(winds) + 0.25, theta_norm[ib][1] - 0.14,
        (names[:b], ANNOTATIONFONTSIZE, :left, PALETTE[ib]))
    y = 0.62
    for i in sortperm([theta_norm[i][end] for i in eachindex(syms)]; rev=true)
        i == ib && continue
        annotate!(f3, last(winds) - 0.25, y,
            (names[syms[i]], ANNOTATIONFONTSIZE, :right, PALETTE[i]))
        y -= 0.10
    end

    # three panels of a third of the text height, half text width; the margins
    # leave room for the axis labels (the theme sets margin = 0)
    fig = plot(f1, f2, f3; layout=(3, 1), size=plot_size(2, 1),
        leftmargin=10Plots.mm, bottommargin=4Plots.mm)
    # the three panels share one x axis: only the bottom panel labels it
    foreach(i -> plot!(fig, subplot=i, xformatter=x -> "", xlabel=""), (1, 2))
    return (; fig, panels=(f1, f2, f3))
end

## Plot: Δθ sensitivity (justifies holding Δθ fixed)
"""
The power and line-tension panels of the Δθ sweep, side by side: Δθ in degrees,
2° to 20° in steps of 2°, which includes the nominal angle. The power is
normalized by its value at the nominal with Θ fixed, so the two conditions meet
at 1 there and the gap between them is what holding Θ costs. Θ held at the
optimum is drawn as black points, Θ re-optimized at each angle as a black line
with triangles (no legend: see the caption);
the tension panel draws both min–max envelopes, the fixed one as points and the
re-optimized one filled. A `guides.line` band marks the nominal angle.
"""
function plot_dtheta(dtheta, dopt, guides, mark)
    ang = [d.Δθ for d in dtheta]
    i0 = argmin(abs.(ang .- mark))
    style = (; c=:black)

    Pn = [d.P for d in dtheta] ./ dtheta[i0].P
    Pn_opt = [d.P for d in dopt] ./ dtheta[i0].P
    Tmin0 = [d.Tmin / 1000 for d in dtheta]
    Tmax0 = [d.Tmax / 1000 for d in dtheta]
    Tmin1 = [d.Tmin / 1000 for d in dopt]
    Tmax1 = [d.Tmax / 1000 for d in dopt]

    # power relative to its value at the nominal Δθ with Θ fixed, so the nominal
    # sits at 1 on both curves and the axis runs from 0 to the peak
    # no legend: GR draws legend markers at 0.8 × the legend font size, whatever
    # `markersize`, so with `thickness_scaling` they dwarf the data; the caption
    # gives the encoding instead
    f1 = plot(ang, Pn; style..., seriestype=:scatter, ms=4, label="",
        ylabel="P / P*", xlabel="Δθ (deg)", ylims=(0, 1.03max(1, maximum(filter(isfinite, Pn_opt); init=1))),
        yticks=([0, 0.5, 1], ["0", "0.5", "1"]), legend=false)
    # markers: a converged angle between two failed ones (NaN) still shows
    plot!(f1, ang, Pn_opt; style..., markershape=:utriangle, ms=4, label="")

    f2 = plot(ang, Tmax0; style..., seriestype=:scatter, ms=4, label="",
        ylabel="Line tension (kN)", xlabel="Δθ (deg)")
    scatter!(f2, ang, Tmin0; style..., ms=4, label="")
    plot!(f2, ang, Tmax1; style..., label="", fillrange=Tmin1, fillalpha=0.1)
    plot!(f2, ang, Tmin1; style..., label="")

    foreach(f -> vline!(f, [mark]; guides.line...), (f1, f2))
    foreach(f -> plot!(f, xticks=2:2:20, xlim=extrema(ang)), (f1, f2))
    # the backend draws the series in list order, so moving the band in front of
    # the data keeps it under the curves
    for f in (f1, f2)
        pushfirst!(f[1].series_list, pop!(f[1].series_list))
    end

    # the margins are per panel (the theme sets margin = 0): only the left one
    # carries a y label and y tick labels
    plot!(f1, leftmargin=18Plots.mm, bottommargin=10Plots.mm)
    plot!(f2, leftmargin=6Plots.mm, bottommargin=10Plots.mm)
    return plot(f1, f2; layout=(1, 2), size=plot_size(1, 2))
end

## Plot: objective landscape around the Θ optimum
"""
Plot a normalized landscape from `landscape`: one row of four panels, one per
parameter, with a single `P / P*` label on the left. Each panel is drawn in the
colour of its parameter in the wind-sweep parameter panel. The x axis is the
parameter in its own units, spanning ±`window` around the optimum, whose value
is the middle tick (`-10%`, value, `+10%`). The y axis is the same in every
panel, `(guides.ymin, 1)`, labelled on the leftmost panel only, with ticks at
`guides.ymin` and 1 (the power curve of a parameter that loses more than 10 % is
clipped at the floor). The cost function (still normalized by its optimum) is a
line, solid inside the box constraint and dotted outside it, the crossing
interpolated linearly between the two samples that bracket the bound. Each box
constraint is a black double vertical line, always drawn (clipped when it falls
outside the window) and on top of the line. A `guides.line` band marks the
optimum vertically and a `guides.floor` line the power floor horizontally; both
are drawn under the data.
"""
function plot_landscape(land, guides)
    (; xs, syms, ys, window, box, center) = land
    lo, hi = 1 - window, 1 + window

    # the margins are per panel (the theme sets margin = 0) and leave room for
    # the tick labels, which sit outside the axis in every panel
    panels = [
        begin
            x0 = TEX[syms[i]][3](center[i])  # the optimum, in display units
            ticks = ([x0 * lo, x0, x0 * hi],
                [f"-{100window:.0f}%", f"{x0:.2f}", f"+{100window:.0f}%"])
            plt = plot(; xlims=(x0 * lo, x0 * hi), ylims=(guides.ymin, 1),
                xticks=ticks, title="$(LABELS[sym])",
                ylabel=(i == 1 ? "P / P*" : ""),
                # only the leftmost panel labels its two ticks
                yticks=(i == 1 ? ([guides.ymin, 1], [plain_tick(guides.ymin), "1"]) :
                        ([guides.ymin, 1], ["", ""])),
                leftmargin=(i == 1 ? 18 : 3)Plots.mm, rightmargin=4Plots.mm,
                bottommargin=5Plots.mm, topmargin=2.5Plots.mm)
            # value of the sampled curve at a factor x, linearly interpolated
            interp(x) = begin
                k = clamp(searchsortedlast(xs, x), 1, length(xs) - 1)
                ys[i][k] + (x - xs[k]) / (xs[k+1] - xs[k]) * (ys[i][k+1] - ys[i][k])
            end
            # factors from a to b, with the ends inserted (they need not be samples)
            piece(a, b) = vcat(a, xs[(a .< xs) .& (xs .< b)], b)
            xa, xb = max(box[i][1], lo), min(box[i][2], hi)  # box ∩ window
            xa < xb && plot!(plt, x0 .* piece(xa, xb), interp.(piece(xa, xb));
                c=i, ls=:solid)
            for (a, b) in ((lo, min(xa, hi)), (max(xb, lo), hi))
                a < b && plot!(plt, x0 .* piece(a, b), interp.(piece(a, b)); c=i, ls=:dot)
            end
            for edge in box[i], k in (-1, 1)  # box constraint: double line, on top
                vline!(plt, [x0 * (edge + k * window / 40)]; c=:black, lw=1.2)
            end
            # every panel draws the shared floor and the optimum line
            hline!([guides.ymin]; guides.floor...)
            vline!([x0]; guides.line...)
            # the backend draws the series in list order, so move the two guide
            # lines in front of the data to keep them underneath it
            series = plt[1].series_list
            pushfirst!(series, pop!(series), pop!(series))
            plt
        end for (i, sym) in enumerate(syms)
    ]

    return plot(panels...; layout=(1, length(syms)), size=plot_size(1/2, 2),
        framestyle=:grid, grid=false)
end

function run_optimization(; fast::Bool=false)
    my_defaults()
    # fast: tiny iteration limits and grids, only to test the workflow (results are not the paper's)
    it_ref, it_sweep, it_dtheta = fast ? (40, 25, 25) : (MAX_ITER, SWEEP_MAX_ITER, DTHETA_MAX_ITER)
    files = String[]

    #region opt_setup
    p_start = PMP.build_para(w0=WIND_REF)
    p_start.b = PMP.optimal_torque_slope(p_start; tol=1e-9, rtol=1e-9, maxiter=50)  # b = I ω
    lower = [first(L_BOUNDS), first(I_BOUNDS), first(B_BOUNDS), first(DPHI_BOUNDS)]
    upper = [last(L_BOUNDS), last(I_BOUNDS), last(B_BOUNDS), last(DPHI_BOUNDS)]

    @info "opt: reference at $(WIND_REF) m/s over Θ = ($(join(SYMS, ", ")))"
    opt_ref = opt_step(p_start, SYMS, lower, upper; iters=it_ref, strict=!fast)
    p_ref = opt_ref.p  # the Θ optimum at the reference wind speed
    # the power of the solver's orbit must be that of the cycle recomputed from p_ref
    # (settling + Newton), which every table and figure uses
    let Pc = cycle_report(p_ref).P, gap = abs(opt_ref.P - Pc) / abs(Pc)
        msg = "optimum: solver power $(opt_ref.P) vs recomputed cycle $Pc (rel. gap $gap)"
        gap < 1e-3 || (fast ? @warn(msg) : error(msg))
    end
    #endregion

    #region opt_wind_sweep
    winds = range(5, 15, length=fast ? 3 : 11)
    i_ref = findfirst(x -> x >= WIND_REF, winds)
    @info "opt: wind sweep, $(length(winds)) speeds from $(first(winds)) to $(last(winds)) m/s"
    sweep = wind_sweep(SYMS, p_ref, winds, lower, upper, i_ref; iters=it_sweep, strict=!fast)
    #endregion

    #region opt_dtheta
    dtheta_angles = collect(range(DTHETA_RANGE..., length=fast ? 3 : 10))
    @info "opt: Δθ sweep with Θ fixed, $(length(dtheta_angles)) angles"
    dtheta = dtheta_study(p_ref, dtheta_angles)

    dsweep = dtheta_opt_sweep(p_ref, SYMS, dtheta_angles, lower, upper;
        nominal=rad2deg(p_ref.Δθ), iters=it_dtheta, step=fast ? 3.0 : DTHETA_STEP)
    dtheta_opt, dtheta_lost = dsweep.stats, dsweep.lost_below
    #endregion

    #region opt_landscape
    LANDSCAPE = landscape(p_ref, SYMS; lower, upper, npoints=fast ? 5 : 21)
    #endregion

    #region opt_tables
    ## Tables
    table_descr, table_opt = latex_tables(p_start, p_ref, lower, upper)
    push!(files, write_table("table2", table_descr), write_table("table3", table_opt))
    @warn "table4: values to be checked vs. Alain's computation"
    push!(files, write_table("table4", comparison_table(p_start, p_ref)))
    push!(files, write_table("numbers", numbers_tex(; p_start, p_ref, opt_iter=opt_ref.iter, winds, sweep, dtheta, dtheta_opt, dtheta_lost)))
    #endregion

    ## Figures
    my_defaults()
    #region fig_wind_sweep
    (; fig, panels) = plot_wind_sweep(winds, sweep, i_ref, GUIDES)
    for (f, panel) in zip(panels, ("power", "tension", "params"))
        plot!(f, size=plot_size(1, 1), leftmargin=10Plots.mm, bottommargin=6Plots.mm)
        append!(files, save_figure(f, "optimization_wind_sweep_$(panel)"))
    end
    #endregion

    my_defaults()
    #region fig_dtheta
    fig_dtheta = plot_dtheta(dtheta, dtheta_opt, GUIDES, rad2deg(p_ref.Δθ))
    append!(files, save_figure(fig_dtheta, "optimization_dtheta"))
    #endregion

    my_defaults()
    #region fig_landscape
    fig_land = plot_landscape(LANDSCAPE, GUIDES)
    append!(files, save_figure(fig_land, "optimization_landscape"))
    #endregion

    my_defaults()
    #region fig_timeseries
    sols = map((p_start, p_ref)) do p
        LC.compute_limit_cycle(PMP.build_vbpara(p); sense=SENSE, tol=TOL, save_everystep=true)
    end
    # initial guess (I = 1300, b = Iω; not the PTO optimum of section C) and the same
    # system optimized for average power
    fig_ts = plot_time_series_grid(sols; titles=("Initial guess", "Dynamic optimum"))
    append!(files, save_figure(fig_ts, "optimization_timeseries_all"))
    #endregion
    return files
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_optimization()
end
