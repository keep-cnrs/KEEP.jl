#= cf. 05.4_functions.jl
What should go in there :
ddq(q, dq, p)
ddq(q, p) = ddq(q, 0, p)
aligned_α(τ, p)

steadiest state from (τ, τ -> α(τ)) 
steady state from (τ, orientation)
steady state from (α, τ)

aligned_and_opposite_steady_states
=#

module SteadyState

using NonlinearSolve
using StaticArrays
using SplitApplyCombine
using QuasiMonteCarlo: sample, HaltonSample
using Statistics: mean
using Clustering: dbscan

import KEEP.PointMass4 as PM4
using KEEP: DEFAULT_TOLERANCE

export ddq, ddq_partial, steady_state, aligned_and_opposite_steady_states
export all_steady_states, all_steady_states_halton

ddq(q, dq, p) = PM4.dynamics(SA[q[1], q[2], dq[1], dq[2], 0], p, 0)[[3, 4]]
ddq_partial(q, p) = ddq(q, (0, 0), p)
aligned_α(τ, p) = PM4.τ_to_θφ(τ, p)[2]


"""
Return (α, τ) that locally minimizes the norm of the acceleration such that α = α_func(τ)"""
function steadiest_τ(α_func, τ_init, p; tol=DEFAULT_TOLERANCE)
    state_from_τ = τ -> SA[α_func(τ), τ]
    opt_func = (τ_arr, p) -> ddq_partial(state_from_τ(only(τ_arr)), p)
    τ_opt = only(solve(NonlinearLeastSquaresProblem(opt_func, [τ_init], p); abstol=tol, reltol=tol).u)
    return rem2pi.(state_from_τ(τ_opt), RoundNearest)
end


"""
Return (α, τ) that locally minimizes the norm of the acceleration such that the arm and the lines are aligned in Oxy. `orientation` controls whether they point in the same or opposite direction, taking values `Val(:same)` or `Val(:opposite)`"""
function steadiest_aligned_state(τ_init, p, orientation::Val{:same}=Val(:same); tol=DEFAULT_TOLERANCE)
    return steadiest_τ(τ -> aligned_α(τ, p), τ_init, p; tol=tol)
end

function steadiest_aligned_state(τ_init, p, orientation::Val{:opposite}; tol=DEFAULT_TOLERANCE)
    return steadiest_τ(τ -> π + aligned_α(τ, p), τ_init, p; tol=tol)
end


"""
Find a steady state starting from q = [α_init, τ_init]"""
function steady_state(α_init::T1, τ_init::T2, p; tol=DEFAULT_TOLERANCE) where {T1, T2}
    T = promote_type(T1, T2)  # couldn't convert IrrationalConstant.fourπ to some reversediff thing
    root = solve(NonlinearProblem(ddq_partial, T[α_init, τ_init], p); abstol=tol, reltol=tol).u
    return rem2pi.(root, RoundNearest)
end


"""
Find a steady state starting from τ

`orientation` is `Val(:same)` or `Val(:opposite)`"""
function steady_state(τ_init, p; orientation=Val(:same), tol=DEFAULT_TOLERANCE)
    α_aligned_, τ_aligned = steadiest_aligned_state(τ_init, p, orientation; tol=tol)
    steady_state(α_aligned_, τ_aligned, p, tol=tol)
end


"""
Find two steady states near τ by:
 1. finding the aligned configuration [α(τ), τ] that minimizes the norm of the acceleration such that the arm and the lines are aligned in Oxy
 2. finding a steady state starting from this aligned configuration"""
function aligned_and_opposite_steady_states(τ_init, p; tol=DEFAULT_TOLERANCE)
    return Tuple(
        steady_state(steadiest_aligned_state(τ_init, p, Val(orientation))..., p; tol=tol)
        for orientation in (:same, :opposite)
    )
end

"""
From a list of steady states, find the unique steady states using clustering."""
function unique_steady_states(steady_states; tol=DEFAULT_TOLERANCE)
    ss_yx = [flatten(sincos.(ss)) for ss in steady_states]
    clustering = dbscan(reduce(hcat, ss_yx), 10tol)
    
    by((y1, x1, y2, x2)) = x1 + tol * y1
    centroids_yx = sort([mean(ss_yx[c.core_indices]) for c in clustering.clusters], by=by)
    centroids = [[atan(c[1], c[2]), atan(c[3], c[4])] for c in centroids_yx]
    return centroids
end

"""
Take unique steady states, augment them using the symmetry of the system, and then take unique steady states again"""
function process_steady_states(steady_states, tol=DEFAULT_TOLERANCE)
    steady_states = unique_steady_states(steady_states; tol=tol)
    steady_states = [[-α, τ - π] for (α, τ) in steady_states]
    steady_states = unique_steady_states(steady_states; tol=tol)
    return steady_states
end

"""
Find all steady states using `aligned_and_opposite_steady_states` where α is guessed from τ"""
function all_steady_states(p; τmin=0, τmax=π, N=100, tol=DEFAULT_TOLERANCE)
    τs = range(τmin, τmax, length=N+2)[2:end-1]
    ss = flatten([aligned_and_opposite_steady_states(τ, p; tol=tol) for τ in τs])
    return process_steady_states(ss)
end

"""
Find all steady states by sampling α and τ from a Halton sequence."""
function all_steady_states_halton(p; αmin=-π, αmax=π, τmin=0, τmax=π, N=100, tol=DEFAULT_TOLERANCE)
    ατs = eachcol(sample(N, Float64[αmin, τmin], Float64[αmax, τmax], HaltonSample()))
    ss = [steady_state(α, τ, p; tol=tol) for (α, τ) in ατs]
    return process_steady_states(ss)
end

end # module