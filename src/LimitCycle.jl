module LimitCycle

#=
Shooting = [α0, dα0, dτ0, T]


What do we want ?
 1. Obtain a limit cycle for some parameters and a tolerance

Solution :
 1. Callback with a closure that remembers the previous activation state and computes the distance to the curent one. If tolerance is satisfied, stop.



What should go in there :

# TODO:
#  - one function that computes the next n passes through τ=0 mod 2π with saveeverystep=false and ContinuousCallback(..., save=true)
#  - one function (elsewhere) that compute the convergence rate of the output of the previous
#  - one function that computes a (limit) cycle between the n-th pass through τ=0 and the next (cf. https://github.com/ljad-cnrs/keep/blob/41bcac967955c0fd5b0be11c631bdba10b8b76c9/test/05.5_LimitCycleFromEquilibrium.jl)
#  - one function that computes a limit cycle by zeroing a shooting function
=#

using DiffEqCallbacks: ContinuousCallback
using OrdinaryDiffEqVerner: terminate!
using StaticArrays
using Clustering: dbscan
using Setfield
using Statistics: mean
using QuasiMonteCarlo: sample, HaltonSample
import ForwardDiff

import KEEP: DEFAULT_TOLERANCE, TAU0
import KEEP.PointMass4 as PM4

export build_poincare_callback,
    distance_on_section,
    build_shooting,
    endpoint_residuals,
    shooting_residuals,
    shoot,
    compute_limit_cycle,
    refine_shooting,
    average_power,
    unpack_shooting,
    all_limit_cycles

"""
Orthogonaly projects a state to the Poincaré section."""
function project_to_section(u)
    return SA[rem2pi(u[1], RoundNearest), u[3], u[4]]
end

"""
The Poincaré section is defined as {x | poincare_section(x) = 0}."""
function poincare_section(τ)
    sin((τ - TAU0) / 2)
end

"""
Compute the distance between two states after projecting them to the Poincaré section."""
function distance_on_section(u1, u2)
    return maximum(abs, project_to_section(u1) - project_to_section(u2))
end

"""
Ideally, tol should be greater than the tolerance used for integrating, like 10 times bigger.

tol = eps(1.): never force termination
tol = Inf: Stop after the first pass

Starting on the Poincarré section will not save the state trough the callback.

In practice with Vern9 it is not needed, but you shouldn't rely on this happy behavior.
Also, tol should be greater than 1e-12 because of numerical errors accumulating during one cycle.

Let tol_int be the tolerance used for integration. Then tol > max(10 * tol_int, 1e-12)."""
function build_affect(tol=10DEFAULT_TOLERANCE)
    function affect!(integrator)
        u = integrator.sol.u
        length(u) < 2 && return
        distance_on_section(u[end-1], u[end]) > tol && return
        terminate!(integrator)
    end
end

"""
Return a callback that saves the states when passing through section = 0. If two subsequent states are within `tol` of each other, terminate the integration."""
function build_poincare_callback(tol=DEFAULT_TOLERANCE; section=poincare_section)
    # `abstol` is not a root-find tolerance: it is the re-trigger guard of
    # DiffEqBase.nudge_tprev. Below ~eps(t) the guard stops working, the same crossing
    # is saved twice, and the shooting then has T = 0. It must stay finite: `tol = Inf`
    # means "stop at the first crossing" and must not become an infinite tolerance.
    ContinuousCallback((τ, t, integrator) -> section(τ), build_affect(tol); save_positions=(true, false), idxs=2, abstol=clamp(tol / 10, 1e-12, 1e-3))
end


"""
shooting = [α0, dα0, dτ0, T] -> u0, T"""
function unpack_shooting(shooting)
    α0, dα0, dτ0, T = shooting
    return SA[rem2pi(α0, RoundNearest), TAU0, dα0, dτ0, 0], T
end

"""
u0, t0, t1 -> shooting = [α0, dα0, dτ0, T]"""
function build_shooting(u0, t0, t1)
    α0, dα0, dτ0 = project_to_section(u0)
    T = t1 - t0
    return SA[α0, dα0, dτ0, T]
end

"""
sol contains only states that are on the Poincaré section, eg. with `build_poincare_callback()`.

sol -> u0, t0, t1 -> shooting = [α0, dα0, dτ0, T]

Return `nothing` when the integration saved less than two distinct crossings (diverging seed, or a trajectory that does not cross the section)."""
function build_shooting(sol)
    length(sol.t) < 2 && return nothing
    u0 = sol.u[end]
    t0, t1 = sol.t[[end - 1, end]]
    t1 > t0 || return nothing
    return build_shooting(u0, t0, t1)
end

"""Guess whether τ should be increasing or decreasing."""
default_sense(shooting) = ifelse(shooting[3] > 0, +, -)
default_sense(u0, u1) = ifelse(u1[2] > u0[2], +, -)


"""
Sense is either `+` or `-`, the expected sign of dτ."""
function endpoint_residuals(u0, u1; sense=default_sense(u0, u1))
    Δα, Δdα, Δdτ = project_to_section(u1 - u0)
    Δτ = u1[2] - u0[2] - sense(2π)
    return SA[Δα, Δτ, Δdα, Δdτ]
end
endpoint_residuals(sol; sense=default_sense(build_shooting(sol))) = endpoint_residuals(sol.u[1], sol.u[end]; sense)

"""
Compute the residuals of a limit cycle from a shooting by integrating. Specify

Sense is either `+` or `-`, the expected sign of dτ."""
function shooting_residuals(shooting, vbp; sense=default_sense(shooting), tol=DEFAULT_TOLERANCE)
    u0, T = unpack_shooting(shooting)
    res = endpoint_residuals(PM4.integrate(u0, T, vbp; tol=tol / 10); sense=sense)
    if maximum(abs, res) < tol
        @warn "LimitCycle.jl:shooting_residuals: Residuals are too low compared to integration tolerance `tol`, it is recommended to reduce `tol` to increase accuracy."
    end
    return res
end

function shoot(shooting, vbp; tol=DEFAULT_TOLERANCE, kwargs...)
    u0, T = unpack_shooting(shooting)
    return PM4.integrate(u0, T, vbp; tol=tol, kwargs...)
end

"""
    refine_shooting(shooting, vbp; sense, tol, maxiter=8, target=tol / 100)

Refine `shooting = [α0, dα0, dτ0, T]` by Newton's method on the closure residual
`endpoint_residuals(shoot(shooting))` (4 equations, 4 unknowns, Jacobian by
ForwardDiff through the integrator), with a backtracking line search on the max
norm of the residual. Stop below `target`, after `maxiter` iterations, or when no
step decreases the residual.

Settling by integration converges only linearly to the cycle, and stops on the
distance between two successive crossings, not on the distance to the cycle: the
shooting it gives closes up to an error that depends on the contraction rate and on
floating-point details. Newton converges quadratically from it, so the refined cycle
closes at the level of the integration tolerance (issue #79).

Return `(shooting, residual)`. The refined shooting is only kept if it decreases the
residual, keeps the sign of `dτ0` and has `T > 0`; otherwise the input is returned,
so that refining never makes a cycle worse."""
function refine_shooting(shooting, vbp; sense, tol=DEFAULT_TOLERANCE, maxiter=8, target=tol / 100)
    F(s) = endpoint_residuals(shoot(s, vbp; tol=tol); sense=sense)
    valid(s) = sign(s[3]) == sign(shooting[3]) && s[4] > 0
    s, r = shooting, F(shooting)
    res = maximum(abs, r)
    for _ in 1:maxiter
        res <= target && break
        δ = try
            ForwardDiff.jacobian(F, s) \ r
        catch
            break  # singular Jacobian (non-hyperbolic cycle) or failed integration
        end
        improved = false
        t = 1.0
        while t >= 1 / 64
            s_new = s - t * δ
            if valid(s_new)
                r_new = F(s_new)
                res_new = maximum(abs, r_new)
                if res_new < res
                    s, r, res, improved = s_new, r_new, res_new, true
                    break
                end
            end
            t /= 2
        end
        improved || break
    end
    return s, res
end

"""
Compute a limite cycle by integrating from u0.

The shooting obtained by settling is refined by [`refine_shooting`](@ref) unless
`refine = false`.

Return its shooting parameters [α0, dα0, dτ0, T]."""
function compute_limit_cycle(u0, vbp; tf=1e9, tol=DEFAULT_TOLERANCE, refine=true, kwargs...)
    cb = build_poincare_callback(10tol)
    sol = PM4.integrate(u0, tf, vbp; callback=cb, save_end=false, tol=tol)
    shooting = build_shooting(sol)
    shooting === nothing && return nothing
    if refine
        shooting, _ = refine_shooting(shooting, vbp; sense=shooting[3] > 0 ? (+) : (-), tol=tol)
    end
    lc = shoot(shooting, vbp; tol=tol, kwargs...)
    # `build_shooting` keeps the last two section crossings, which are unrelated when the
    # trajectory has not settled: check that what comes back closes on itself, in the
    # direction it actually travels, before handing it out as a limit cycle.
    sense = shooting[3] > 0 ? (+) : (-)
    res = maximum(abs, endpoint_residuals(lc; sense=sense))
    if res > max(100tol, 1e-6)
        @error "compute_limit_cycle: the trajectory is not a limit cycle" res T = lc.t[end] Δτ = lc.u[end][2] - lc.u[1][2] dτ0 = shooting[3] tol
        return nothing
    end
    return lc
end

"""
Compute the attracting limit cycle with `sense(dτ) > 0`, by integrating from the
Poincaré-section state `u0 = (α, τ, dα, dτ) = (0, TAU0, 0, sense(1))`.

Both pumping cycles are attracting, so the sign of `dτ` selects one and the
magnitude of the seed only needs to be a reasonable guess. Use `all_limit_cycles`
for an inventory of every cycle in a box; the two calls are independent."""
function compute_limit_cycle(vbp; sense, tol=DEFAULT_TOLERANCE, kwargs...)
    u0 = PM4.build_u(; dτ=sense(1.0))
    lc = compute_limit_cycle(u0, vbp; tol=tol, kwargs...)
    lc === nothing &&
        error("compute_limit_cycle: no closed limit cycle with dτ = $(sense(1.0)) (never settled, or not periodic)")
    sense(lc.u[1][4]) > 0 || error("compute_limit_cycle: converged to a cycle with dτ of the wrong sign")
    return lc
end


"""
shooting -> average power"""
function average_power(shooting, vbp; tol=DEFAULT_TOLERANCE)
    lc = shoot(shooting, vbp; tol=tol)
    return lc.u[end][end] / lc.t[end]
end

"""
Comptue all limit cycles by sampling `N` initial conditions within a region of the Poincarré section.
 - A rough integration with low precision (sqrt(tol))
 - Clusterize and take centroids
 - Integrate each centroid with high precision (tol)"""
function all_limit_cycles(vbp; αmin=0., αmax=1π, vmax=10, N=100, tol=DEFAULT_TOLERANCE, kwargs...)
    dαmin, dαmax = dτmin, dτmax = -vmax, vmax

    samples = sample(N, [αmin, dαmin, dτmin], [αmax, dαmax, dτmax], HaltonSample())
    u0s = [SA[α, TAU0, dα, dτ, 0] for (α, dα, dτ) in eachcol(samples)]
    # A seed that diverges, or whose trajectory does not close, returns `nothing` and
    # the clustering cannot use it.
    rough = [compute_limit_cycle(u0, vbp; tol=sqrt(tol), refine=false) for u0 in u0s]
    rough_shootings = filter(!isnothing, [lc === nothing ? nothing : build_shooting(lc) for lc in rough])
    isempty(rough_shootings) && error("no limit cycle seed converged")

    clusters = dbscan(reduce(hcat, rough_shootings), 100sqrt(tol))
    rough_centroids = [mean(rough_shootings[c.core_indices]) for c in clusters.clusters]

    return filter(!isnothing,
        [compute_limit_cycle(unpack_shooting(s)[1], vbp; tol=tol, kwargs...) for s in rough_centroids])
end

end  # module