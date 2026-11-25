# Optimization section — compute only.
#
# The optimization variables are Θ = (ℓ, I, b, Δφ): the tether length, the arm
# inertia, the braking coefficient and the lemniscate half-width. The
# lemniscate half-height Δθ is held at its nominal value; the last study here
# checks what that choice costs (power, line tension, kite acceleration) and a
# 10×10 map of the same quantities over the shape plane (Δθ, Δφ).
#
# The reference is the default parameters with the torque slope set by the
# kinetic criterion b = I ω (`PMP.optimal_torque_slope`, one limit cycle per
# iteration), not the default torque slope. From there:
#  - OPTIMIZATION: maximize the average power over Θ at 9 m/s, then repeat the
#    optimization at each wind speed and compare with the 9 m/s optimum held
#    fixed. ℓ spans 15–200 m and b 1–1e6; the Δφ box is (10°, 90°), and the
#    inertia floor is set by sizing the aluminium arm (see `I_BOUNDS`).
#  - Δθ: sweep it from 2° to 20° at the optimum with Θ fixed, then re-optimize
#    Θ at each angle, and compare the power and the line-tension envelope with
#    and without re-optimization.
#  - SHAPE: 10×10 map of the same quantities over (Δθ, Δφ) at the optimal
#    (ℓ, I, b), again without re-optimizing anything.
#
# All numerical results (including the arrays needed for the figures) are
# returned by `compute_optimization()`. Plot styling, figure notes and all
# figure writes live in `paper.jl`, which is the only runner.

using KEEP: PointMass4 as PM4
using KEEP: PointMassPara as PMP
using KEEP: LimitCycle as LC
using KEEP: Optimization as OPT

using PyFormattedStrings: @f_str
using SplitApplyCombine: invert
using Logging: with_logger, NullLogger
using LinearAlgebra: norm

## Settings
const TOL = 1e-5
const MAX_WALL_TIME = 120.0
const SWEEP_WALL_TIME = 30.0  # a converged wind-sweep solve takes about a second
# `WIND_REF` is a shared optimization setting defined in `paper_core.jl`.
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
const SHAPE_WINDOW = 0.2                      # ±20 % for the shape map
const GRID_N = 10

# Parameter names, units and scales from the stored value to the displayed one.
# An angle is rad2deg; a per-angle unit is the inverse, so that b·dα stays a
# torque when dα is read in deg/s.
const TEX = Dict(:l => (raw"$\ell$", raw"\si{\meter}", identity),
    :I => (raw"$I$", raw"\si{\kilogram\meter\squared}", identity),
    :b => (raw"$b$", raw"\si{\newton\meter\second\per\degree}", x -> x / rad2deg(1)),
    :Δφ => (raw"$\Delta\varphi$", raw"\si{\degree}", rad2deg))

# Parameter names and units for the console tables
const LABELS = Dict(:l => "ℓ", :I => "I", :b => "b",
    :Δθ => "Δθ", :Δφ => "Δφ")
const UNITS = Dict(:l => "m", :I => "kg m²", :b => "N m deg s⁻¹",
    :Δθ => "deg", :Δφ => "deg")

## Helpers
"""
Average power, line-tension extrema, peak kite speed and acceleration, and period
of a cycle."""
function cycle_stats(sol; n=2001)
    # `compute_limit_cycle` returns `nothing` when no periodic orbit was found: report
    # NaN instead of crashing, so the caller treats the point as a failure.
    sol === nothing && return (; P=NaN, Tmin=NaN, Tmax=NaN, amax=NaN, vmax=NaN, period=NaN)
    ts = range(sol.t[1], sol.t[end], n)
    us = sol.(ts)
    p = sol.prob.p
    tension = PM4.compute_line_tension.(us, Ref(p))
    ok = isfinite.(tension)
    acceleration = norm.(PM4.kite_acceleration.(us, Ref(p)))
    velocity(u) = PM4.jacobian(q -> PM4.compute_OK(PM4.compute_Rτ(q, p), p), u[1:2]) * u[3:4]
    vmax = maximum(norm.(velocity.(us)))
    return (; P=sol.u[end][5] / sol.t[end],
        Tmin=minimum(tension[ok]), Tmax=maximum(tension[ok]),
        amax=maximum(acceleration), vmax=vmax, period=sol.t[end])
end

# Limit cycle, optimization, and limit cycle of the optimized system, for one
# value of one parameter. `lower`/`upper` are the box constraints of `syms`.
# `p_fixed` is the parameter set the no-optimization cycle is taken at, so a
# continuation can start from a neighbour's optimum `p0` without moving the
# reference curve. `guess` is Ipopt's initial guess (see `clamp_guess`), which
# lets a retry start from the neighbouring cycle states instead of cold.
function compute_results(param, value; syms, p0, p_fixed=p0, lower, upper, sense, tol,
    wall=MAX_WALL_TIME, guess=(;))
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
        OPT.optimize(p0, syms, lower, upper; sense, tol, max_wall_time=wall,
            initial_guess=guess)
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
function wind_sweep(syms, p_ref, winds, lower, upper, i_ref; rtol=1e-3)
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
                wall=SWEEP_WALL_TIME, guess=guess)
            r[4].status == :first_order || continue
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
One optimization step with the compilation warm-up: the first solve with a
given variable layout pays Julia's compilation of the dual-number Jacobian
inside ipopt's wall clock, which can eat the time cap and truncate the solve.
Warm the layout up and keep the second run.
"""
function opt_step(p, syms, lower, upper)
    solve() = with_logger(NullLogger()) do
        OPT.optimize(p, syms, lower, upper; sense=SENSE, tol=TOL, max_wall_time=MAX_WALL_TIME)
    end
    solve()
    t_wall = time()
    solution, stats, _ = solve()
    stats.status == :first_order || @warn "opt_step($(join(syms, ", "))) did not converge" stats.status stats.iter
    p_next = copy(p)
    p_next[syms] .= solution.params[syms]  # physical units
    return (; p=p_next, P=solution.P, iter=stats.iter, t_wall=time() - t_wall)
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
optimum); `walls` are the per-attempt wall-clock caps. A point where every start
fails is reported as `NaN`, not as a gain. Statistics come from `dtheta_study`
at the solved point, i.e. from the same dτ > 0 cycle as the fixed curve.
"""
function dtheta_opt_sweep(p_ref, syms, angles_deg, lower, upper; i0,
    walls=(60.0, 60.0, 60.0, 90.0), rtol=1e-3)
    n = length(angles_deg)
    stats = Vector{Any}(undef, n)
    p_at = Vector{Any}(undef, n)
    sols = Vector{Any}(undef, n)
    function attempt(i, p0, g, wall)
        q = copy(p0)
        q.Δθ = deg2rad(angles_deg[i])
        try
            solution, st, _ = with_logger(NullLogger()) do
                OPT.optimize(q, syms, lower, upper; sense=SENSE, tol=TOL,
                    max_wall_time=wall, initial_guess=g)
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
            r = attempt(i, start, g, walls[k])
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
10×10 shape map at fixed (ℓ, I, b): one limit cycle per (Δθ, Δφ), no
re-optimization. Rows are Δφ, columns are Δθ factors. A cell stays NaN when no
dτ > 0 pumping cycle is found.
"""
function shape_grid(p, factors, dphi_deg; n=401)
    u0 = PM4.build_u(; dτ=1)
    P = fill(NaN, length(dphi_deg), length(factors))
    Tmin = fill(NaN, size(P))
    Tmax = fill(NaN, size(P))
    amax = fill(NaN, size(P))
    for (j, dphi) in enumerate(dphi_deg), (i, f) in enumerate(factors)
        i == 1 && @info "  shape map row $(j)/$(length(dphi_deg)): Δφ = $(dphi)°"
        q = copy(p)
        q.Δθ = f * p.Δθ
        q.Δφ = deg2rad(dphi)
        vbp = PMP.build_vbpara(q)
        sol = try
            LC.compute_limit_cycle(u0, vbp; tol=TOL, save_everystep=true, tf=100.0)
        catch
            nothing
        end
        sol === nothing && continue
        # reject a trajectory that did not close
        norm(LC.endpoint_residuals(sol; sense=SENSE)) < 1e-3 || continue
        s = cycle_stats(sol; n)
        P[j, i], Tmin[j, i], Tmax[j, i], amax[j, i] =
            s.P / 1000, s.Tmin / 1000, s.Tmax / 1000, s.amax
    end
    return (; P, Tmin, Tmax, amax, factors, dphi_deg)
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
function print_latex_tables(p_init, p_opt, lower, upper)
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
    push!(results_rows, f"    Max kite speed (\\si{{\\meter\\per\\second}}) & {s0.vmax:.2f} & {s1.vmax:.2f} \\\\ ")
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

    println(join(["%% BEGIN TABLE 2", table_descr, table_opt, "%% END TABLE 2"], "\n\n"))
end

"""
Run every numerical study of the optimization section and print the console
and LaTeX tables. Every array needed by the figures in `paper.jl` is returned.

Returned fields:
- `p_start`: initial parameter set (the 9 m/s PTO-sizing point, b = I ω).
- `lower`: lower box bounds of Θ = (ℓ, I, b, Δφ), in physical units.
- `upper`: upper box bounds of Θ, in physical units.
- `opt_ref`: reference optimization result `(p, P, iter, t_wall)` at 9 m/s.
- `p_ref`: optimized parameter set at the reference wind speed (the Θ optimum).
- `winds`: wind speeds (m/s) of the wind sweep, 5 to 15 in 11 steps.
- `i_ref`: index in `winds` of the reference wind speed.
- `sweep`: wind-sweep data `(; syms, no_opt, opt, theta)`; each entry is a cycle
  statistics NamedTuple and `theta[k]` is the optimized Θ at that wind speed.
- `dtheta_angles`: absolute Δθ angles (deg) of the Δθ sweeps, 2 to 20 in 10 steps.
- `dtheta`: cycle statistics with Θ held at the optimum, one per Δθ angle.
- `dtheta_opt`: cycle statistics with Θ re-optimized at each Δθ angle.
- `shape`: shape-map data `(; P, Tmin, Tmax, amax, factors, dphi_deg)`, 10×10,
  rows Δφ, columns Δθ factors; powers in kW, tensions in kN, acceleration m/s².
- `P_map`: average power (kW) of the 10×10 shape map (`shape.P`).
- `Tmin_map`: minimum line tension (kN) of the shape map (`shape.Tmin`).
- `Tmax_map`: maximum line tension (kN) of the shape map (`shape.Tmax`).
- `amax_map`: peak kite acceleration (m/s²) of the shape map (`shape.amax`).
- `x_map`: Δθ deviation (%) on the shape heatmaps, `100 .* (factors .- 1)`.
- `shape_dphi`: Δφ axis (deg) of the shape heatmaps (`shape.dphi_deg`).
- `LANDSCAPE`: objective landscape `(; xs, syms, ys, window, center, box)`, one
  row of normalized powers per optimized parameter around the Θ optimum.
- `sols`: the two limit-cycle solutions at (`p_start`, `p_ref`) for time series.
- `names`: the series names `(:init, :opt)` matching `sols`.
"""
function compute_optimization()
    ## Numerical: 9 m/s reference, one optimization over Θ
    p_start = PMP.build_para(w0=WIND_REF)
    p_start.b = PMP.optimal_torque_slope(p_start)  # b = I ω
    lower = [first(L_BOUNDS), first(I_BOUNDS), first(B_BOUNDS), first(DPHI_BOUNDS)]
    upper = [last(L_BOUNDS), last(I_BOUNDS), last(B_BOUNDS), last(DPHI_BOUNDS)]

    @info "opt: reference at $(WIND_REF) m/s over Θ = ($(join(SYMS, ", ")))"
    opt_ref = opt_step(p_start, SYMS, lower, upper)
    p_ref = opt_ref.p  # the Θ optimum at the reference wind speed

    ## Numerical: wind sweep over Θ, from the 9 m/s optimum
    winds = range(5, 15, length=11)
    i_ref = findfirst(x -> x >= WIND_REF, winds)
    @info "opt: wind sweep, $(length(winds)) speeds from $(first(winds)) to $(last(winds)) m/s"
    sweep = wind_sweep(SYMS, p_ref, winds, lower, upper, i_ref)

    ## Numerical: Δθ sensitivity at the optimum, with the nominal angle included
    # 2° to 20° in steps of 2°, which already contains the nominal (6°)
    dtheta_angles = collect(range(DTHETA_RANGE..., length=10))
    @info "opt: Δθ sweep with Θ fixed, $(length(dtheta_angles)) angles"
    dtheta = dtheta_study(p_ref, dtheta_angles)

    ## Numerical: Δθ sweep with Θ re-optimized at each angle
    @info "opt: Δθ sweep with Θ re-optimized, $(length(dtheta_angles)) angles"
    dtheta_opt = dtheta_opt_sweep(p_ref, SYMS, dtheta_angles, lower, upper;
        i0=argmin(abs.(dtheta_angles .- rad2deg(p_ref.Δθ))))

    ## Numerical: 10×10 shape map at the optimal ℓ, I, b
    @info "opt: shape map $(GRID_N)×$(GRID_N) at the optimal ℓ, I, b"
    shape = shape_grid(p_ref,
        range(1 - SHAPE_WINDOW, 1 + SHAPE_WINDOW, length=GRID_N),
        range(rad2deg.(DPHI_BOUNDS)..., length=GRID_N))  # shape_grid takes degrees

    ## Numerical: objective landscape around the Θ optimum
    @info "opt: objective landscape around the Θ optimum"
    LANDSCAPE = landscape(p_ref, SYMS; lower, upper)

    ## Numerical: time-series cycles, initial guess against optimum
    names = (:init, :opt)
    sols = map((p_start, p_ref)) do p
        LC.compute_limit_cycle(PMP.build_vbpara(p); sense=SENSE, tol=TOL, save_everystep=true)
    end

    ## Console report
    print_final_table(p_start, p_ref, winds, sweep, dtheta, dtheta_opt, shape, LANDSCAPE)
    print_latex_tables(p_start, p_ref, lower, upper)

    # Heatmap axes, precomputed for `paper.jl` (see `plot_shape_grid`)
    x_map = 100 .* (shape.factors .- 1)  # Δθ deviation (%) on the shape heatmaps
    return (; p_start, lower, upper, opt_ref, p_ref, winds, i_ref, sweep,
        dtheta_angles, dtheta, dtheta_opt, shape,
        P_map=shape.P, Tmin_map=shape.Tmin, Tmax_map=shape.Tmax,
        amax_map=shape.amax, x_map, shape_dphi=shape.dphi_deg,
        LANDSCAPE, sols, names)
end

"""
Print the console tables of the optimization section: final parameters and
limit-cycle statistics, the wind sweep, the Δθ sweep, the shape map and the
objective landscape."""
function print_final_table(p_start, p_ref, winds, sweep, dtheta, dtheta_opt, shape, LANDSCAPE)
    stats = (init=cycle_stats(LC.compute_limit_cycle(PMP.build_vbpara(p_start);
        sense=SENSE, tol=TOL, save_everystep=true)),
        opt=cycle_stats(LC.compute_limit_cycle(PMP.build_vbpara(p_ref);
            sense=SENSE, tol=TOL, save_everystep=true)))
    cols = ((p_start, stats.init), (p_ref, stats.opt))
    fmt(x) = string(round(x; sigdigits=6))
    row(name, f) = println(rpad(name, 30), join([lpad(f(c), 14) for c in cols]))

    println("\nfinal parameters")
    println(rpad("", 30), join(lpad.(["init", "Θ"], 14)))
    for sym in SYMS
        row("$(LABELS[sym]) ($(UNITS[sym]))", c -> fmt(TEX[sym][3](c[1][sym])))
    end
    row("Δθ (deg, fixed)", c -> fmt(rad2deg(c[1].Δθ)))
    println()
    row("period (s)", c -> fmt(c[2].period))
    row("objective P (kW)", c -> fmt(c[2].P / 1000))
    row("min tension (kN)", c -> fmt(c[2].Tmin / 1000))
    row("max tension (kN)", c -> fmt(c[2].Tmax / 1000))
    row("max acceleration (m/s²)", c -> fmt(c[2].amax))
    row("max kite speed (m/s)", c -> fmt(c[2].vmax))
    println()

    println("wind sweep over Θ (reference at $(round(Int, WIND_REF)) m/s)")
    println(rpad("v (m/s)", 8), rpad("fixed P", 10), rpad("opt P", 10), rpad("gain", 8),
        rpad("Tmax opt", 10), rpad("ℓ (m)", 9), rpad("I", 9), rpad("b", 9), "Δφ (deg)")
    for k in eachindex(winds)
        P0, P = sweep.no_opt[k].P / 1000, sweep.opt[k].P / 1000
        θ = sweep.theta[k]
        println(rpad(round(Int, winds[k]), 8), rpad(round(P0; digits=3), 10),
            rpad(round(P; digits=3), 10), rpad(f"{(P / P0 - 1) * 100:+.1f} %", 8),
            rpad(round(sweep.opt[k].Tmax / 1000; digits=2), 10),
            rpad(round(θ[1]; digits=2), 9), rpad(round(θ[2]; digits=1), 9),
            rpad(round(TEX[:b][3](θ[3]); digits=1), 9), round(rad2deg(θ[4]); digits=2))
    end
    println()

    println("Δθ sweep: Θ fixed against re-optimized (nominal Δθ = $(round(rad2deg(p_ref.Δθ); digits=1)) deg)")
    println(rpad("Δθ (deg)", 10), rpad("P0 (kW)", 9), rpad("P (kW)", 9), rpad("gain", 8),
        rpad("Tmax0 (kN)", 12), rpad("Tmax (kN)", 11), rpad("amax0", 9), "amax")
    for (d0, d) in zip(dtheta, dtheta_opt)
        gain = (d.P / d0.P - 1) * 100
        println(rpad(f"{d0.Δθ:.1f}", 10), rpad(f"{d0.P / 1000:.2f}", 9),
            rpad(f"{d.P / 1000:.2f}", 9),
            rpad(abs(gain) < 0.5 ? "—" : f"{gain:+.0f} %", 8),
            rpad(f"{d0.Tmax / 1000:.2f}", 12), rpad(f"{d.Tmax / 1000:.2f}", 11),
            rpad(f"{d0.amax:.0f}", 9), f"{d.amax:.0f}")
    end
    println()

    println("shape map (no re-optimization, $(GRID_N)×$(GRID_N) limit cycles)")
    finite(Z) = filter(isfinite, Z)
    println("  valid: ", count(isfinite, shape.P), " of ", length(shape.P), " cells")
    for (name, Z, unit) in (("P", shape.P, "kW"), ("Tmax", shape.Tmax, "kN"),
        ("Tmin", shape.Tmin, "kN"))
        lo, hi = extrema(finite(Z))
        println("  ", rpad(name, 5), round(lo; digits=2), " … ", round(hi; digits=2), " ", unit)
    end
    k = argmax(ifelse.(isfinite.(shape.amax), shape.amax, -Inf))
    println("  max a: ", round(shape.amax[k]; digits=0), " m/s² at Δθ = ",
        round((shape.factors[k[2]] - 1) * 100; digits=0), " %, Δφ = ",
        round(shape.dphi_deg[k[1]]; digits=1), " deg")
    println()

    println("landscape at the Θ optimum (power relative to the optimum)")
    println(rpad("parameter", 12), rpad("at -10 %", 12), rpad("at +10 %", 12), "box (÷/×)")
    for (i, sym) in enumerate(SYMS)
        println(rpad(LABELS[sym], 12), rpad(round(LANDSCAPE.ys[i][1]; digits=4), 12),
            rpad(round(LANDSCAPE.ys[i][end]; digits=4), 12),
            round.(LANDSCAPE.box[i]; digits=3))
    end
end
