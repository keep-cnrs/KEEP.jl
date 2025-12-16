using NonlinearSolve
using SplitApplyCombine
using StaticArrays
using Statistics
using Plots
using Distances: pairwise, Euclidean
using LinearAlgebra: norm, I
using Test
using QuasiMonteCarlo: sample, SobolSample

import KEEP.PointMass4 as PM4

using KEEP.PointMassPara

include("steady_state_utils.jl")

default(label="")  # Plots default
tol = 1e-12
vbp = build_vbpara()

## One aligned and one opposite steady state from one τ value
τ_init =  .5
ss = aligned_and_opposite_steady_states(τ_init, vbp; tol=tol)
@test all(norm.(ddq.(ss, Ref(vbp))) .< tol)

## Aligned and opposite steady states from multiple τ values
N = 100
τmin, τmax = 0, π
τs = range(τmin, τmax, length=N)[2:end-1]

complete_steady_states = combinedims(flatten(aligned_and_opposite_steady_states.(τs, Ref(vbp); tol=tol)))

steady_states = reduce_expand_reduce(complete_steady_states, 10tol)

# Verify that they are all equilibriums
@test all(norm.(ddq.(eachcol(steady_states), Ref(vbp))) .< tol)

# Verify that these equilibriums are all different
@test all(compute_adjacency(steady_states, tol) - I .== 0)

# sort steady states by φ
sort_inds_φ = sortperm([aligned_α(τ, vbp) for (_, τ) in eachcol(steady_states)])
sort_inds_α = sortperm([-cos(α) for (α, _) in eachcol(steady_states)])
steady_states .= steady_states[:, sort_inds_α]

## Steady states from a grid of α and τ
αmin, αmax = -π, π
αs, τs = eachrow(sample(N, [αmin, τmin], Float64[αmax, τmax], SobolSample()))
complete_steady_states_QM = combinedims(steady_state.(αs, τs, Ref(vbp); tol=tol))
steady_states_QM = reduce_expand_reduce(complete_steady_states_QM, 10tol)
steady_states_QM .= steady_states_QM[:, sortperm([-cos(α) for (α, _) in eachcol(steady_states_QM)])]

## Do we get the same steady states from τ only and from both α and τ ? YES !
@test length(steady_states) == length(steady_states_QM)
@test norm(minimum(pairwise(Euclidean(), steady_states, steady_states_QM; dims=2), dims=2)) < sqrt(tol)  # All equal up to order

# true to plot one at a time, by hand
i = 0
p_plot = init_plot2d(vbp)
if false
    add_point_plot2d!(steady_states[:, @show i+=1], p_plot)
end

# Plot all at once
p_plot = init_plot2d(vbp)
for ss in reverse(eachcol(steady_states))
    local col, ls, lw = ifelse(cos(ss[1]) > 0,
            (:green, :solid, 3),
            (:red, :dash, 2))
    add_point_plot2d!(ss, p_plot, color=col, linestyle=ls, linewidth=lw, marker=:blue, markersize=4)
end
plot!(title="8 equilibriums")
display(plot!())
# savefig("docs/media/equilibriums.pdf")

split_ss = [filter(ss -> cos(ss[1]) > 0, eachcol(steady_states)), filter(ss -> cos(ss[1]) < 0, eachcol(steady_states))]

# Animation positive x/negative x
if false
    anim = @animate for steady_states in split_ss
        @show steady_states
        p_plot = init_plot2d(vbp)
        col, ls, lw = ifelse(cos(steady_states[1][1]) > 0,
                (:green, :solid, 3),
                (:red, :solid, 2))
        add_point_plot2d!.(steady_states, Ref(p_plot), color=col, linestyle=ls, linewidth=lw, marker=:blue, markersize=4)
        plot!(title="8 equilibriums")
    end
    display(gif(anim, fps=1/2))
    # gif(anim, "docs/media/equilibriums-clustered.gif"; fps=1/2)

    # Animation one by one
    anim = @animate for (i, ss) in enumerate(eachcol(steady_states))
        p_plot = init_plot2d(vbp)
        add_point_plot2d!(ss, p_plot, color=ifelse(cos(ss[1]) > 0, :green, :red), marker=:blue, lw=2)
        plot!(title="Equilibrium $i")
    end
    display(gif(anim, fps=1))
    # gif(anim, "docs/media/equilibriums.gif"; fps=1)
end
