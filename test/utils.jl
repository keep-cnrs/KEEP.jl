using LinearAlgebra: norm
using Plots
using Logging
using StrFormat
using SplitApplyCombine

using KEEP.TorqueFunction
import KEEP.PointMass4 as PM4
import KEEP.PointMass10 as PM10

"""Pretty string of an object"""
function pretty_string(obj)
    io = IOBuffer()
    show(io, "text/plain", obj)
    String(take!(io))
end

"""arr has samples in columns and variables in rows. Returns the relative negative log norm difference between each pair of samples"""
function rel_neg_log_norm_diff(arr)
    norms = map(norm, eachslice(arr, dims=2))
    arr1 = reshape(arr, size(arr)..., 1)
    arr2 = permutedims(arr1, (1, 3, 2))
    diff_norms = map(norm, eachslice(arr1 .- arr2, dims=(2, 3)))
    @. -log10(diff_norms / max(norms, norms'))
end

function should_verbose()
    return Logging.min_enabled_level(current_logger()) <= Info
end

"""
Quadratic decay outside [-Ωmax, Ωmax] instead of instantaneous
"""
function advanced_torque_func(dα, p)
    (; Ωmin, Ωmax, Ωlim, Cmax) = p
    abs_out = begin
        if Ωmin < abs(dα) < Ωmax
            (Cmax / (Ωmax - Ωmin)) * (abs(dα) - Ωmin)
        elseif Ωmax < abs(dα) < Ωlim
            Cmax * ((abs(dα) - Ωlim) / (Ωmax - Ωlim))^10
        else
            zero(dα)
        end
    end
    return sign(dα) * abs_out
end

function vecvec2mat(vecvec)
    return reduce(hcat, vecvec)
end

# function plot_position10(sol, p; t=sol.t, figure=plot())
#     q = sol.(t; idxs=1:5)
#     A = PM10.compute_posA.(q, (p,))
#     K1 = PM10.compute_pos1.(q, (p,))
#     K2 = PM10.compute_pos2.(q, (p,))

#     # indices : x/y/z, O/A/K2, indice du segment
#     segments = stack(vecvec2mat.((zero(A), A, K2)); dims=2)

#     plot!(figure)
#     plot!(xlabel=raw"$x$", ylabel=raw"$y$", zlabel=raw"$z$", legend=:topleft, size=(600, 600))
#     plot!(eachslice(segments; dims=1)..., c=:black, alpha=20 / length(t), label="")
#     plot!(Tuple.(A), lw=1, label=raw"$A(α)$", c=1)
#     plot!(Tuple.(K1), lw=1, label=raw"$K_1(R, τ)$", c=2)
#     plot!(Tuple.(K2), lw=1, label=raw"$K_2(α, θ2, φ2)$", c=3)
#     return plot!()
# end


# # TODO: Have one function that initializes the plot, and one function that adds a position to it, then loop over all positions
# # Reuse these for the animation
# function plot_position4(sol, p; t=sol.t, figure=plot())
#     q = sol.(t; idxs=1:5)
#     α = first.(q)
#     A = p.l .* PM4.compute_αhat.(α)
#     K = PM4.compute_OK.(PM4.compute_Rτ.(q, (p,)), (p,))

#     # indices : x/y/z, O/A/K2, indice du segment
#     segments = stack(vecvec2mat.((zero(A), A, K)); dims=2)
#     plot!(figure)
#     plot!(xlabel=raw"$x$", ylabel=raw"$y$", zlabel=raw"$z$", legend=:topleft, size=(600, 600))
#     plot!(eachslice(segments; dims=1)..., c=:black, alpha=20 / length(t), label="")
#     plot!(Tuple.(A), label=raw"$A(α)$", c=1)
#     plot!(Tuple.(K), label=raw"$K(α, τ)$", c=2)
#     return plot!()
# end


# function plot_avg_power10(sol, p; t=sol.t, figure=plot())
#     dα_ind, P_ind = 8, 11
#     t0, tf = extrema(t)

#     dα = sol.(t, idxs=dα_ind)
#     power = dα .* torque_function.(dα, (p,))
#     avg_power = sol(tf, idxs=P_ind) / (tf - t0)
#     plot!(figure, t, power, c=:black, lw=2, label="")
#     plot!([0, tf], avg_power .* [1, 1], c=:red, lw=2, label=f"Average power = \%.0f(avg_power) W")
#     plot!(legend=:outerbottom, xlabel="t (s)", ylabel="Power (W)")
# end


# function plot_avg_power4(sol, p; t=sol.t, figure=plot(), kwargs...)
#     dα_ind, P_ind = 3, 5
#     L, M, T = lmt(p)
#     t0, tf = extrema(t)

#     dα = sol.(t, idxs=dα_ind)
#     power = dα .* torque_function.(dα .* T, (p,)) .* (M * L^2 * T^-2)
#     avg_power = sol(tf, idxs=P_ind) / (tf - t0)
#     plot!(figure, t, power, c=:black, lw=2, label="")
#     plot!([0, tf], avg_power .* [1, 1], c=:red, lw=2, label=f"Average power = \%.0f(avg_power) W")
#     default_kwargs = (legend=:outerbottom, xlabel="t (s)", ylabel="Power (W)")
#     return plot!(; merge(default_kwargs, kwargs)...)
# end


# function animate_position4(sol, p; tspan=extrema(sol.t), fps=30, trail_length=1, trail_step=1 / 5fps, start_camera=(30, 30), end_camera=(40, 30), figsize=(800, 800))
#     # Compute box
#     K = PM4.compute_OK.(PM4.compute_Rτ.(sol.u, (p,)), (p,))
#     xmax, ymax, zmax = 1.1 * maximum(K)

#     t0, tf = tspan
#     anim = @animate for t_ in t0:(1/fps):tf
#         t_trail = [max(0, t_ - trail_length):trail_step:t_; t_]
#         q = sol.(t_trail; idxs=1:5)
#         α = first(q[end])
#         A = p.l * PM4.compute_αhat(α)
#         trail_K = PM4.compute_OK.(PM4.compute_Rτ.(q, (p,)), (p,))
#         segment = [zero(A), A, trail_K[end]] |> vecvec2mat

#         width = range(0.0, 1.0, length=length(t_trail))
#         τ = (t_ - t0) / (tf - t0)
#         camera = start_camera .* (1 - τ) .+ end_camera .* τ
#         plot(eachrow(segment)..., lw=1, c=:black, label="")
#         plot!(Tuple.(trail_K), msw=0, lw=4 * width, alpha=width, label="")
#         scatter!(Tuple.([A]), msw=0, ms=4, label="")
#         plot!(xlims=(0, xmax), ylims=(-ymax, ymax), zlims=(0, zmax),
#             xlabel="x", ylabel="y", zlabel="z",
#             camera=camera,
#             titlefont="monospace", title=f"t = \%.2f(t_) s / \%.2f(tf) s",
#             size=figsize
#         )
#     end fps = fps
#     return anim
# end


function pairwise_dist(X)
    G = X * X'
    sq_norms = sum(abs2, X, dims=2)
    D_sq = @. sq_norms + sq_norms' - 2 * G
    @. D_sq = sqrt(max(D_sq, 0.0))

    return D_sq
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
    cluster_centers = [mean(steady_states[:, clusters.==i_cluster], dims=2) for i_cluster in 1:nb_clusters]
    cluster_centers = combinedims(vec.(cluster_centers))

    # Add the complement of each center
    complement_α = -cluster_centers[1:1, :]
    complement_τ = π .+ cluster_centers[2:2, :]
    reduced_steady_states = [cluster_centers [complement_α; complement_τ]]

    # Clusterize the complemented centers
    reduced_steady_states .= rem2pi.(reduced_steady_states, RoundNearest)
    reduced_clusters, nb_reduced_clusters = connected_components(compute_adjacency(reduced_steady_states, ε))

    # Compute the new centers
    reduced_cluster_centers = [mean(reduced_steady_states[:, reduced_clusters.==i_cluster], dims=2) for i_cluster in 1:nb_reduced_clusters]
    reduced_cluster_centers = combinedims(vec.(reduced_cluster_centers))

    return reduced_cluster_centers
end
