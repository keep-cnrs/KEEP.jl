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
using OrdinaryDiffEqTsit5: terminate!
using StrFormat
using StaticArrays
using Clustering: dbscan
using Setfield
using Statistics: mean
using QuasiMonteCarlo: sample, HaltonSample

import KEEP: DEFAULT_TOLERANCE, TAU0
import KEEP.PointMass4 as PM4

export build_poincare_callback,
    distance_on_section,
    build_shooting,
    endpoint_residuals,
    shooting_residuals,
    shoot,
    compute_limit_cycle,
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

In practice with Tsit5 it is not needed, but you shouldn't rely on this happy behavior.
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
    ContinuousCallback((τ, t, integrator) -> section(τ), build_affect(tol); save_positions=(true, false), idxs=2, abstol=tol / 10)
end


"""
shooting -> u0, T"""
function unpack_shooting(shooting)
    α0, dα0, dτ0, T = shooting
    return SA[rem2pi(α0, RoundNearest), TAU0, dα0, dτ0, 0], T
end

"""
u0, t0, t1 -> shooting"""
function build_shooting(u0, t0, t1)
    α0, dα0, dτ0 = project_to_section(u0)
    T = t1 - t0
    return SA[α0, dα0, dτ0, T]
end

"""
sol contains only states that are on the Poincaré section, eg. with `build_poincare_callback()`.

sol -> u0, t0, t1 -> shooting"""
function build_shooting(sol)
    u0 = sol.u[end]
    t0, t1 = sol.t[[end - 1, end]]
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
Compute a limite cycle by integrating from u0.

Return its shooting parameters [α0, dα0, dτ0, T]."""
function compute_limit_cycle(u0, vbp; tf=1e9, tol=DEFAULT_TOLERANCE, kwargs...)
    cb = build_poincare_callback(10tol)
    sol = PM4.integrate(u0, tf, vbp; callback=cb, save_end=false, tol=tol)
    return shoot(build_shooting(sol), vbp; tol=tol, kwargs...)
end
function compute_limit_cycle(vbp; sense, tol=DEFAULT_TOLERANCE, kwargs...)
    all_lcs = all_limit_cycles(vbp; tol=tol, kwargs...)
    return first(filter(lc -> sense(lc.u[1][4]) > 0, all_lcs))
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
function all_limit_cycles(vbp; αmin=0., αmax=1π, vmax=5, N=100, tol=DEFAULT_TOLERANCE, kwargs...)
    dαmin, dαmax = dτmin, dτmax = -vmax, vmax

    samples = sample(N, [αmin, dαmin, dτmin], [αmax, dαmax, dτmax], HaltonSample())
    u0s = [SA[α, TAU0, dα, dτ, 0] for (α, dα, dτ) in eachcol(samples)]
    rough_shootings = [build_shooting(compute_limit_cycle(u0, vbp; tol=sqrt(tol))) for u0 in u0s]

    clusters = dbscan(reduce(hcat, rough_shootings), 100sqrt(tol))
    rough_centroids = [mean(rough_shootings[c.core_indices]) for c in clusters.clusters]

    limit_cycles = [compute_limit_cycle(unpack_shooting(s)[1], vbp; tol=tol, kwargs...) for s in rough_centroids]
    return limit_cycles
end

end  # module