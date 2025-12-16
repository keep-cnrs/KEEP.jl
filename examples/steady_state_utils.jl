# using NonlinearSolve
using SplitApplyCombine
using ADTypes: AutoForwardDiff
using NonlinearSolve: NewtonRaphson
# using StaticArrays
# using Statistics
# using Distances: pairwise, Euclidean
# using Plots

import KEEP.PointMass4 as PM4
import KEEP: DEFAULT_TOLERANCE

function ddq(q, p)
    u = zeros(eltype(q), 5)
    u[[1, 2]] = q
    return PM4.dynamics(u, p, 0)[[3, 4]]
end

function aligned_α(τ, p)
    return PM4.τ_to_θφ(τ, p)[2]
end

# TODO: annuler uniquement ddα ou ddτ avec uniquement α ou τ ?

"""
Return (α, τ) that locally minimizes the norm of the acceleration such that α = α_func(τ)"""
function steadiest_τ(α_func, τ_init, p; tol=DEFAULT_TOLERANCE)
    state_from_τ = τ -> SA[α_func(τ), τ]
    opt_func = (τ_arr, p) -> ddq(state_from_τ(only(τ_arr)), p)
    τ_opt = only(solve(NonlinearLeastSquaresProblem(opt_func, [τ_init], p); abstol=tol, reltol=tol).u)
    return rem2pi.(state_from_τ(τ_opt), RoundNearest)
end

"""
Return (α, τ) that locally minimizes the norm of the acceleration such that the arm and the lines are aligned in Oxy. `orientation` controls whether they point in the same or opposite direction, taking values `Val(:same)` or `Val(:opposite)`"""
function steadiest_aligned_state(τ_init, p, orientation::Val{:same}=Val(:same); tol=DEFAULT_TOLERANCE)
    return steadiest_τ(τ -> aligned_α(τ, p), τ_init, p; tol=tol)
end

function steadiest_aligned_state(τ_init, p, orientation::Val{:opposite}; tol=DEFAULT_TOLERANCE)
    return steadiest_τ(τ -> π+aligned_α(τ, p), τ_init, p; tol=tol)
end

"""
Find a steady state starting from q = [α_init, τ_init]"""
function steady_state(α_init, τ_init, p; tol=DEFAULT_TOLERANCE)
    root = solve(NonlinearProblem(ddq, [α_init, τ_init], p); abstol=tol, reltol=tol).u
    return rem2pi.(root, RoundNearest)
end

"""
Find a steady state starting from τ

orientation is `Val(:same)` or `Val(:opposite)`"""
function steady_state(τ_init, p; orientation=Val(:same), tol=DEFAULT_TOLERANCE)
    α_aligned_, τ_aligned = steadiest_aligned_state(τ_init, p, orientation; tol=tol)
    steady_state(α_aligned_, τ_aligned, p, tol=tol)
end

"""
Find two steady states near τ by:
 1. finding the aligned configuration [α(τ), τ] that minimizes the norm of the acceleration such that the arm and the lines are aligned in Oxy
 2. finding a steady state starting from this aligned configuration"""
function aligned_and_opposite_steady_states(τ_init, p; tol=DEFAULT_TOLERANCE)
    # Compact version
    return Tuple(steady_state(steadiest_aligned_state(τ_init, p, Val(orientation))..., p; tol=tol) for orientation in (:same, :opposite))

    # # Non-compact version
    # α_same_init, τ_same_init = steadiest_aligned_state(τ_init, p, Val(:same), tol=tol)
    # ss_same = steady_state(α_same_init, τ_same_init, p, tol=tol)

    # α_opposite_init, τ_opposite_init = steadiest_aligned_state(τ_init, p, Val(:opposite), tol=tol)
    # ss_opposite = steady_state(α_opposite_init, τ_opposite_init, p, tol=tol)

    # return (ss_same, ss_opposite)
end

# ## Limit Cycle
# """
# given the parameters p, find an equilibrium state q0, add a small negative speed in dτ and integrate 8 cycles to find the limit cycle"""
# function get_limit_cycle(p, dτ_sign::Union{typeof(+), typeof(-)})
#     q0 = steady_state(1, p)
#     q0[4] = dτ_sign(1)
#     cb1 = ContinuousCallback((u, t, integrator) -> u[2] - dτ_sign(8π), identity, save)
#     cb2 = DiscreteCallback((u, t, integrator) -> u[2] - dτ_sign(10π), terminate!)


    ## Plots

function init_plot2d(p)
    # Eight
    θ0 = p.l+p.r
    Δφ = θ0 * p.Δφ
    Δθ = p.Δθ * Δφ / p.Δφ
    p_plot = build_para(θ0 = θ0, Δθ = Δθ, φ0 = 0, Δφ = Δφ)
    τ = range(-π, π, length=100)
    θ, φ = invert(PM4.τ_to_θφ.(τ, Ref(p_plot)))
    plot(φ, θ, aspect_ratio=:equal, color=:black, grid=false, axis=false, label="")

    # Circle
    α = range(-π, π, length=100)
    x, y = p.l .* invert(sincos.(α))
    plot!(y, x; color=:black, label="")

    l = p_plot.l
    offset = l / 4
    
    # ordre: y, x, φ, θ
    x_axes = [-l, -l - offset, p_plot.φ0 - p_plot.Δφ, p_plot.φ0 - p_plot.Δφ - offset]
    y_axes = [-l - offset, -l, p_plot.θ0 - p_plot.Δθ - offset, p_plot.θ0 - p_plot.Δθ]
    
    u_axes = [2l, 0, 2l, 0]
    v_axes = [0, 2l, 0, 2l]
    xy_axes = hcat(x_axes, y_axes)
    uv_axes = hcat(u_axes, v_axes)
    quiver!(eachcol(xy_axes)..., quiver=Tuple(eachcol(uv_axes)), c=:black)
    
    x_lbl = @. x_axes + u_axes / 2 - v_axes / 4
    y_lbl = @. y_axes + v_axes / 2 - u_axes / 4
    labels = text.(["y", "x", "φ", "θ"], :black, :center, 10)
    annotate!(x_lbl, y_lbl, labels)
    plot!([-1, -2p.l], alpha=0, label="")  # Else "y" is not shown
    
    return p_plot
end

function add_point_plot2d!(q, p_plot, plot_args...; plot_kwargs...)
    α, τ = q
    θ, φ = PM4.τ_to_θφ(τ, p_plot)
    x, y, _ = p_plot.l * PM4.compute_αhat(α)
    return plot!([0, y, φ], [0, x, θ]; plot_args..., plot_kwargs...)
end

## Graphs

"""
Compute the distance between vectors of angles by first mapping them to the unit circle."""
function compute_adjacency(angular_samples, ε)
    samples = [flatten([sincos(α) for α in angles]) for angles in eachcol(angular_samples)]
    pairwise(Euclidean(), samples) .< ε
end

"""
Depth first search to find connected components from a binary adjacency matrix."""
function connected_components(adjacency)
    n = size(adjacency, 1)
    visited = falses(n)
    component_index = zeros(Int, n)
    next_component = 1

    function dfs!(v, current_component)
        visited[v] = true
        component_index[v] = current_component
        for neighbor in 1:n
            if adjacency[v, neighbor] && !visited[neighbor]
                dfs!(neighbor, current_component)
            end
        end
    end

    for v in 1:n
        if !visited[v]
            dfs!(v, next_component)
            next_component += 1
        end
    end
    return component_index, next_component - 1
end

"""
Add complement steady states and remove duplicates.

The complement of (α, τ) is (-α, -τ).
Remove duplicates by keeping only the centers of each connected component, two states being considered connected when their euclidean distance is less than ε.

1. Reduce: Clusterize the steady states
2. Expand: Add the complement of each cluster center
3. Reduce: Remove potential duplicates (complements that were already in the centers)"""
function reduce_expand_reduce(steady_states, ε)
    # Clusterize
    steady_states .= rem2pi.(steady_states, RoundNearest)
    clusters, nb_clusters = connected_components(compute_adjacency(steady_states, ε))
    
    # Compute the center of each cluster
    cluster_centers = [mean(steady_states[:, clusters .== i_cluster], dims=2) for i_cluster in 1:nb_clusters]
    cluster_centers = combinedims(vec.(cluster_centers))
    
    # Add the complement of each center
    complement_α = -cluster_centers[1:1, :]
    complement_τ = π .+ cluster_centers[2:2, :]
    reduced_steady_states = [cluster_centers [complement_α; complement_τ]]
    
    # Clusterize the complemented centers
    reduced_steady_states .= rem2pi.(reduced_steady_states, RoundNearest)
    reduced_clusters, nb_reduced_clusters = connected_components(compute_adjacency(reduced_steady_states, ε))
    
    # Compute the new centers
    reduced_cluster_centers = [mean(reduced_steady_states[:, reduced_clusters .== i_cluster], dims=2) for i_cluster in 1:nb_reduced_clusters]
    reduced_cluster_centers = combinedims(vec.(reduced_cluster_centers))
    
    return reduced_cluster_centers
end