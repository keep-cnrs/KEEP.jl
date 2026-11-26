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
using LinearAlgebra: norm

include(joinpath(@__DIR__, "helpers.jl"))


## Settings
const TOL = 1e-5
# Ipopt iteration limits (a converged solve takes a few dozen iterations); the
# wall-clock time is only a safety net that is not expected to bind.
const MAX_ITER = 300
const SWEEP_MAX_ITER = 150
const DTHETA_MAX_ITER = (150, 150, 150, 300)  # one limit per start, see `dtheta_opt_sweep`
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
    ok = isfinite.(tension)
    acceleration = norm.(PM4.kite_acceleration.(us, Ref(p)))
    return (; P=sol.u[end][5] / sol.t[end],
        Tmin=minimum(tension[ok]), Tmax=maximum(tension[ok]),
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

    # Limit cycle at the optimized parameters
    vbp = PMP.build_vbpara(p_opt)
    lc = LC.compute_limit_cycle(vbp; sense=sense, save_everystep=true, tol=tol)
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
function opt_step(p, syms, lower, upper; iters=MAX_ITER)
    solution, stats, _ = with_logger(NullLogger()) do
        OPT.optimize(p, syms, lower, upper; sense=SENSE, tol=TOL, max_iter=iters,
            max_wall_time=WALL_SAFETY)
    end
    stats.status == :first_order || @warn "opt_step($(join(syms, ", "))) did not converge" stats.status stats.iter
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
Δθ sweep with re-optimization: at each absolute angle (degrees) in `angles_deg`,
set Δθ and re-optimize `syms` inside `lower`/`upper`, marching outward from the
nominal angle `i0`. The objective is multimodal and a single start is not
reliable: 4° needs the nominal cycle states as initial guess, and at 2° a plain
start from 4° converges to a box-corner local optimum (7.3 kW against 12.8). So
each angle tries, in order, the neighbour's optimum with the neighbour's cycle
states, the same parameters cold, `p_ref` with the nominal's states, and `p_ref`
cold (the cycle-state guess can be worse than none at a distant angle), keeps
the best converged average power, and stops as soon as one beats the fixed-`Θ`
value at that angle (which is always feasible, so nothing below it can be the
optimum); `iters` are the per-attempt iteration limits. A point where every start
fails is reported as `NaN`, not as a gain. Statistics come from `dtheta_study`
at the solved point, i.e. from the same dτ > 0 cycle as the fixed curve.
"""
function dtheta_opt_sweep(p_ref, syms, angles_deg, lower, upper; i0,
    iters=DTHETA_MAX_ITER, rtol=1e-3)
    n = length(angles_deg)
    stats = Vector{Any}(undef, n)
    p_at = Vector{Any}(undef, n)
    sols = Vector{Any}(undef, n)
    function attempt(i, p0, g, max_iter)
        q = copy(p0)
        q.Δθ = deg2rad(angles_deg[i])
        try
            solution, st, _ = with_logger(NullLogger()) do
                OPT.optimize(q, syms, lower, upper; sense=SENSE, tol=TOL,
                    max_iter=max_iter, max_wall_time=WALL_SAFETY, initial_guess=g)
            end
            p = copy(q)
            p[syms] .= solution.params[syms]
            return (p=p, status=st.status, sol=solution, P=solution.P)
        catch
            return (p=copy(p0), status=:error, sol=nothing, P=NaN)
        end
    end
    function run!(i, prev)
        @info "  Δθ re-optimization $(i)/$(n): $(angles_deg[i])°"
        fixed = only(dtheta_study(p_ref, [angles_deg[i]]))   # Θ held at the nominal
        starts = prev === nothing ? ((p_ref, (;)),) :
            ((prev.p, clamp_guess(prev.sol)), (prev.p, (;)),
             (p_ref, clamp_guess(sols[i0])), (p_ref, (;)))
        best = nothing
        for (k, (start, g)) in enumerate(starts)
            r = attempt(i, start, g, iters[k])
            if r.status == :first_order &&
               (best === nothing || best.status != :first_order || r.P > best.P)
                best = r
            elseif best === nothing
                best = r
            end
            k == length(starts) && break
            # Holding Θ at the nominal is always feasible, so a converged point below
            # the fixed-Θ value at this angle is another branch: try the next start.
            best.status == :first_order && best.P >= fixed.P * (1 - rtol) && break
        end
        if best.status != :first_order || best.P < fixed.P * (1 - rtol)
            @warn "re-optimization at Δθ = $(angles_deg[i])° does not beat the fixed Θ ($(round(fixed.P / 1000; digits=3)) kW)" best.status best.P
        end
        p_at[i], sols[i] = best.p, best.sol
        # A timed-out solve leaves half-converged parameters: report a gap rather than a
        # bogus power for a point that was never solved.
        stats[i] = best.status == :first_order ? only(dtheta_study(best.p, [angles_deg[i]])) :
            (; Δθ=angles_deg[i], P=NaN, Tmin=NaN, Tmax=NaN, amax=NaN, period=NaN)
        return nothing
    end
    run!(i0, nothing)
    for i in (i0-1):-1:1
        run!(i, (p=p_at[i+1], sol=sols[i+1]))       # walk left from the nominal
    end
    for i in (i0+1):n                               # restart at the nominal
        run!(i, (p=(i == i0+1 ? p_ref : p_at[i-1]), sol=(i == i0+1 ? sols[i0] : sols[i-1])))
    end
    return stats
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
Table 4: initial guess ("PTO initialization", the power-take-off sizing point) against
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
        raw"    Parameter & \multicolumn{1}{c}{PTO initialization} & \multicolumn{1}{c}{Dynamic optimization} \\",
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
optimum is drawn as black points, Θ re-optimized at each angle as a black line;
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
    f1 = plot(ang, Pn; style..., seriestype=:scatter, ms=4, label="Θ fixed",
        ylabel="P / P*", xlabel="Δθ (deg)", ylims=(0, 1.03max(1, maximum(Pn_opt))),
        yticks=([0, 0.5, 1], ["0", "0.5", "1"]), legend=:topright)
    plot!(f1, ang, Pn_opt; style..., label="Θ re-optimized")

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
    it_ref, it_sweep, it_dtheta = fast ? (40, 25, (25, 25, 25, 25)) : (MAX_ITER, SWEEP_MAX_ITER, DTHETA_MAX_ITER)
    files = String[]

    p_start = PMP.build_para(w0=WIND_REF)
    p_start.b = PMP.optimal_torque_slope(p_start)  # b = I ω
    lower = [first(L_BOUNDS), first(I_BOUNDS), first(B_BOUNDS), first(DPHI_BOUNDS)]
    upper = [last(L_BOUNDS), last(I_BOUNDS), last(B_BOUNDS), last(DPHI_BOUNDS)]

    @info "opt: reference at $(WIND_REF) m/s over Θ = ($(join(SYMS, ", ")))"
    opt_ref = opt_step(p_start, SYMS, lower, upper; iters=it_ref)
    p_ref = opt_ref.p  # the Θ optimum at the reference wind speed

    winds = range(5, 15, length=fast ? 3 : 11)
    i_ref = findfirst(x -> x >= WIND_REF, winds)
    @info "opt: wind sweep, $(length(winds)) speeds from $(first(winds)) to $(last(winds)) m/s"
    sweep = wind_sweep(SYMS, p_ref, winds, lower, upper, i_ref; iters=it_sweep, strict=!fast)

    dtheta_angles = collect(range(DTHETA_RANGE..., length=fast ? 3 : 10))
    @info "opt: Δθ sweep with Θ fixed, $(length(dtheta_angles)) angles"
    dtheta = dtheta_study(p_ref, dtheta_angles)

    dtheta_opt = dtheta_opt_sweep(p_ref, SYMS, dtheta_angles, lower, upper;
        i0=argmin(abs.(dtheta_angles .- rad2deg(p_ref.Δθ))), iters=it_dtheta)

    LANDSCAPE = landscape(p_ref, SYMS; lower, upper, npoints=fast ? 5 : 21)

    ## Tables
    table_descr, table_opt = latex_tables(p_start, p_ref, lower, upper)
    push!(files, write_table("table2", table_descr), write_table("table3", table_opt))
    @warn "table4: values to be checked vs. Alain's computation"
    push!(files, write_table("table4", comparison_table(p_start, p_ref)))

    ## Figures
    my_defaults()
    (; fig, panels) = plot_wind_sweep(winds, sweep, i_ref, GUIDES)
    for (f, panel) in zip(panels, ("power", "tension", "params"))
        plot!(f, size=plot_size(1, 1), leftmargin=10Plots.mm, bottommargin=6Plots.mm)
        append!(files, save_figure(f, "optimization_wind_sweep_$(panel)"))
    end

    my_defaults()
    fig_dtheta = plot_dtheta(dtheta, dtheta_opt, GUIDES, rad2deg(p_ref.Δθ))
    append!(files, save_figure(fig_dtheta, "optimization_dtheta"))

    my_defaults()
    fig_land = plot_landscape(LANDSCAPE, GUIDES)
    append!(files, save_figure(fig_land, "optimization_landscape"))

    my_defaults()
    sols = map((p_start, p_ref)) do p
        LC.compute_limit_cycle(PMP.build_vbpara(p); sense=SENSE, tol=TOL, save_everystep=true)
    end
    # "PTO optimum" is the initial guess (the power-take-off sizing point);
    # "Dynamic optimum" is the same system optimized for average power.
    fig_ts = plot_time_series_grid(sols; titles=("PTO optimum", "Dynamic optimum"))
    append!(files, save_figure(fig_ts, "optimization_timeseries_all"))
    return files
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_optimization()
end
