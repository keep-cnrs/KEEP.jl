# L'initialisation de l'état dans l'ancien code (2*10d) supposait que α et dα étaient nuls.
# Pour tester la dynamique pour α et dα non-nuls, on intègre la dynamique puis on fait les tests à un temps quelconque.

using Logging
using Test
using Clustering
using Plots
using MultivariateStats
using LinearAlgebra: norm
using StaticArrays
using StatsBase

import KEEP.PointMass10 as PM10
import KEEP.PointMass4 as PM4
using KEEP.PointMassPara

include("utils.jl");

## Create the initial state
p = build_para()
vbp = build_vbpara(p)
vbp_normed = normalize_vbpara(vbp)
p_normed = build_para(vbp_normed)

L, M, T = lmt(p)

# indices of 4d state in the 10d state
inds_10d = SA[3, 2, 8, 7, 11]
state_dim_10d = [L, 1, 1, 1, 1, L / T, 1 / T, 1 / T, 1 / T, 1 / T, M * L^2 * T^-2]
state_dim_4d = state_dim_10d[inds_10d]

α = 0.
τ = 1.
dα = 0.
dτ = 3.

u0_4d = SA[α, τ, dα, dτ, 0.]
u0_10d = PM10.init_u(τ, dτ, p)
u0_4d_normed = u0_4d ./ state_dim_4d
u0_10d_normed = u0_10d ./ state_dim_10d
# u0_10d_normed = PM10.init_u(u0_4d_normed[[2, 3, 4]]..., p_normed)

# Check that the two ways of creating u0_10d_normed are consistent
@test u0_10d_normed ≈ PM10.init_u(τ, dτ * T, p_normed)

## Compare at t=0
# 10d & 10d_normed
# 4d & 4d_normed
# 4dVB & 4dVB_normed
# dy[time of comparison]_[code used][VBPara? less parameters][normalised? no dimensions]
du0_10d = PM10.dynamics!(similar(u0_10d), u0_10d, p, 0.)[inds_10d]
du0_10d_normed = PM10.dynamics!(similar(u0_10d_normed), u0_10d_normed, p_normed, 0.)[inds_10d] .* state_dim_4d / T  # PM10 using Para with normalisation
du0_4dVB = PM4.dynamics(u0_4d, vbp, 0.)  # PM4norm using VBPara without normalisation
du0_4dVB_normed = PM4.dynamics(u0_4d_normed, vbp_normed, 0.) .* state_dim_4d / T  # PM4norm using VBPara with normalisation

# equivalent to `@test all(eachcol(du0s) .≈ (du0s[:, 1],))`
@test du0_10d ≈ du0_10d_normed
@test du0_10d ≈ du0_4dVB
@test du0_10d ≈ du0_4dVB_normed

# Logarithm of the relative difference
names = ["10d" "10d_normed" "4dVB" "4dVB_normed"]
du0s = [du0_10d du0_10d_normed du0_4dVB du0_4dVB_normed]
str = join(split(pretty_string(rel_neg_log_norm_diff(du0s)), '\n')[2:end], '\n')
@info """[Dynamique à l'état initial]
Nombre de chiffres significatifs communs entre les dynamiques (1 pas)
[10d, 10d normalisé, 4d, 4d normalisé] \n\n""" * str


## Compare at t=1
# dα(t=1) is sufficient to reach energy production

# Generate the state at t=1. Use 10d because there is no code
# for 4d -> 10d, and 10d -> 4d is a simple projection
tf = 20.
cb = PM10.build_manifold_projection(u0_10d)
sol_cb = PM10.integrate(u0_10d, tf, p; callback=cb)
y1_10d_cb = sol_cb(tf)

y1_4d = SA[y1_10d_cb[inds_10d]...]
y1_10d_normed = y1_10d_cb ./ state_dim_10d
y1_4d_normed = y1_4d ./ state_dim_4d


# Compare 10d and 4d at t=20
dy1_10d = PM10.dynamics!(similar(y1_10d_cb), y1_10d_cb, p, tf)[inds_10d]
dy1_10d_normed = PM10.dynamics!(similar(y1_10d_normed), y1_10d_normed, p_normed, tf)[inds_10d] .* state_dim_4d / T
dy1_4dVB = PM4.dynamics(y1_4d, vbp, tf)
dy1_4dVB_normed = PM4.dynamics(y1_4d_normed, vbp_normed, tf) .* state_dim_4d / T

@test dy1_10d ≈ dy1_10d_normed
@test dy1_10d ≈ dy1_4dVB
@test dy1_10d ≈ dy1_4dVB_normed

dy1s = [dy1_10d dy1_10d_normed dy1_4dVB dy1_4dVB_normed]
str = join(split(pretty_string(rel_neg_log_norm_diff(dy1s)), '\n')[2:end], '\n')
relres_cb = norm(PM10.manifold_residuals!(similar(u0_10d, 6), y1_10d_cb, p)) / norm(y1_10d_cb)
@info """[Dynamique à un état "dense" le long d'une trajectoire obtenue par intégration]
Nombre de chiffres significatifs communs entre les dynamiques à t=20 (1 pas depuis la solution du modèle 10d)
[10d, 10d normalisé, 4d, 4d normalisé] \n\n""" * str * "\n\n[Validation du point de calcul, obtenu par intégration avec callback]
Norme relative des résidus pour le point utilisé : " * string(relres_cb)


## Test de l'intégration de la dynamique
# 10d & 10d_normed
# 4d & 4d_normed
# 4dVB & 4dVB_normed

tf_normed = tf / T

sol10d = PM10.integrate(u0_10d, tf, p)
y1_10d = sol10d(tf)[inds_10d]

sol10d_normed = PM10.integrate(u0_10d_normed, tf_normed, p_normed)
y1_10d_normed = sol10d_normed(tf / T)[inds_10d] .* state_dim_4d

sol4dVB = PM4.integrate(u0_4d, tf, vbp)
y1_4dVB = sol4dVB(tf)

sol4dVB_normed = PM4.integrate(u0_4d_normed, tf_normed, vbp_normed)
y1_4dVB_normed = sol4dVB_normed(tf / T) .* state_dim_4d

@test y1_10d ≈ y1_10d_normed
@test y1_10d ≈ y1_4dVB
@test y1_10d ≈ y1_4dVB_normed

y1s = collect([y1_10d y1_10d_normed y1_4dVB y1_4dVB_normed])
str = join(split(pretty_string(rel_neg_log_norm_diff(y1s)), '\n')[2:end], '\n')
@info """[États après intégration]
Nombre de chiffres significatifs communs entre les états après une intégration de 20 secondes :
[10d, 10d normalisé, 4d, 4d normalisé] \n\n""" * str

relres_no_cb = norm(PM10.manifold_residuals!(similar(u0_10d, 6), sol10d(tf), p)) / norm(sol10d(tf))
@info "[Résidus callback vs no callback]
Résidus à t=20 avec callback : " * string(relres_cb) * " (point utilisé plus haut)
Résidus à t=20 sans callback : " * string(relres_no_cb)
@test relres_cb < relres_no_cb

## Visualisation : quels méthodes donnent des résultats similaires entre eux ?


"""
        generate_logsim_and_order(samples)

    Calculate the log-similarity of the provided sample matrix and order the samples by similarity.

    # Arguments
    - `samples::Matrix`: A matrix where each column represents a different sample.

    # Returns
    - `logsim::Matrix`: A matrix representing the log-similarity values between samples.
    - `order::Vector{Int}`: An ordering of the samples that groups similar samples close to one another.

    # Methodology
    The function computes the pairwise Euclidean distance between columns of the `samples` matrix.
    It then converts these distances into log-similarity values and normalizes them to the range [0, 1].
    Finally, it uses hierarchical clustering to order the samples by similarity.
"""
function generate_logsim_and_order(samples)
    # Calculate the pairwise distance between columns of the samples matrix
    dist = [norm(samp1 - samp2) for samp1 in eachcol(samples), samp2 in eachcol(samples)]

    # Transform distances to log-similarity values
    logsim = -log10.(dist .+ minimum(dist[dist.>0.]))
    min_logsim, max_logsim = extrema(logsim)
    logsim = (logsim .- min_logsim) ./ (max_logsim - min_logsim)

    # Order the samples by similarity using hierarchical clustering
    order = hclust(1 .- logsim; branchorder=:r).order

    return logsim, order
end

"""
        logsim_heatmap(logsim, order, names)

    Create a heatmap of the log-similarity matrix ordered by similarity.

    # Arguments
    - `logsim`: Matrix of log-similarity values.
    - `order`: Ordering vector to reorder the rows and columns of `logsim`.
    - `names`: Names corresponding to the samples to label the axes.
    - `title`: Title of the heatmap.

    # Examples
    ```julia
    logsim_heatmap(logsim, order, names)
    ```
"""
function logsim_heatmap(logsim, order, names, title="")
    # Reorder log-similarity matrix and names according to the ordering
    ordered_logsim = logsim[order, order]
    ordered_names = names[order]

    heatmap(ordered_logsim,
        xticks=(1:length(ordered_names), ordered_names),
        yticks=(1:length(ordered_names), ordered_names),
        xmirror=true,
        xrotation=45,
        yflip=true,
        color=:Purples,
        title=title,
        aspect_ratio=1,
        size=(600, 600)
    )
end

## Pour voir les intensités en passant la souris sur le heatmap
# plotlyjs()
# gr()

## Heatmap pour la dynamique
logsim_du0, order_du0 = generate_logsim_and_order(du0s)
logsim_dy1, order_dy1 = generate_logsim_and_order(dy1s)

display(logsim_heatmap(logsim_du0, order_dy1, names, "Dynamique à t=0 — log-similarité"))
display(logsim_heatmap(logsim_dy1, order_dy1, names, "Dynamique à t=20 — log-similarité"))

## Heatmap sur l'état final de l'intégration
logsim_y1, order_y1 = generate_logsim_and_order(y1s)
display(logsim_heatmap(logsim_y1, order_dy1, names, "État à t=20 — log-similarité"))

## PCA sur l'état final de l'intégration

# rescale les données avec StatsBase
rescaler = fit(ZScoreTransform, y1s, dims=2)
y1s_rescaled = StatsBase.transform(rescaler, y1s)
pca = fit(PCA, y1s_rescaled, mean=0, pratio=1)
pc1, pc2 = eachrow(predict(pca, y1s_rescaled))

min_pc1, max_pc1 = extrema(pc1)
min_pc2, max_pc2 = extrema(pc2)
# Normalized PCA projection
# 10d_normed is omitted because it is an outlier
x_starts = 0.8 .* range(min_pc1, max_pc1, length(names))
y_starts = zeros(length(names))
pca_plot = plot([x_starts pc1]', [y_starts pc2]', labels=names, c=(1:length(names))', legend=:outerright)
scatter!(pca_plot, collect(pc1)', collect(pc2)', c=(1:length(names))', msw=0, labels="")
plot!(pca_plot, title="PCA de l'état à t=20", xlabel="PC1", ylabel="PC2")
display(pca_plot)
