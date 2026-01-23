#=
Define a shooting by :
    - α, dα, dτ the starting position
    - T the period

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

using StaticArrays: SA
using Setfield
using Test
using Logging
using StrFormat

import KEEP.PointMass4 as PM4
import KEEP.Visualisation as VIS
using KEEP.PointMassPara: build_para, build_vbpara, VBPara
using KEEP.LimitCycle
using KEEP: DEFAULT_TOLERANCE, TAU0

# using Plots
# default(formatter=:plain, label="", lw=3)

# xmax = 50
# plot()
# hline([0], c=:grey, lw=1)
# vline!([TAU0] .+ (0:2π:xmax), c=:grey, lw=1, label = "τ = TAU0 mod 2π")
# plot!(section, 0, xmax, label="Poincare section")
# # Visually check that section is 0 when τ ≡ TAU0 mod 2π


tol = 1e-9
u0 = SA[0.0, TAU0, 0.0, 0.0, 0.0]
vbp = build_vbpara()
tf = 100
cb = build_poincare_callback(eps(1.))
sol = PM4.integrate(u0, tf, vbp; callback=cb, save_end=false, tol=tol)


## Plot convergence toward a fixed point in the Poincaré section
u = sol.u[1:end-1]
u_ref = sol.u[end]
errors = [distance_on_section(u, u_ref) for u in u]

# Visually, decreases linearly (in log-plot) between i=1 to 6
i1, i2 = 2, 5
convergence_rate = -(log10(errors[i2]) - log10(errors[i1])) / (i2 - i1)

@info "Convergence rate = $(convergence_rate)"

using Plots
default(label="", lw=3)
plot(title="Convergence of the Poincaré map", xlabel="iteration", ylabel="distance to last iterate", yscale=:log10, yticks=exp10.(-16:2:0))
scatter!(errors, m=:x)
plot!([i1, i2], [errors[i1], errors[i2]], c=:red, lw=3, label=f"Convergence rate = \%.2f(convergence_rate)")
display(plot!())

plot!(xlim=(0, 10))
plot!(title="Convergence of the Poincaré map — zoom")
display(plot!())



#=
If needed:
 - a function that given u0, uf and tf, computes the residuals (no integration)
 - a function that given a shooting, computes the residuals (with integration), integrate then use the previous function
=#

# sol is a limit cycle shooting with precision tol
# When integrating, it should verify α0 = αf, ...
# This should be correct up to 2 tol: tol for u0 + tol for uf (approximation?)
shooting = build_shooting(sol)
@test maximum(abs, shooting_residuals(shooting, vbp; sense=-, tol=tol / 2)) < tol  # It is a limit cycle with negative dτ


# Visualise the path
lc = shoot(shooting, vbp; tol=tol, save_everystep=true)
VIS.plot_trajectory_4D(lc, points_per_second=100)


## Finding all limit cycles
using Clustering: dbscan
using StrFormat
using Statistics: mean
using QuasiMonteCarlo: sample, HaltonSample

αmin, αmax = 0., 1π
vmax = 5.
dαmin, dαmax = dτmin, dτmax = -vmax, vmax

samples = sample(100, [αmin, dαmin, dτmin], [αmax, dαmax, dτmax], HaltonSample())
u0s = [SA[α, TAU0, dα, dτ, 0] for (α, dα, dτ) in eachcol(samples)]
shootings = [build_shooting(compute_limit_cycle(u0, vbp; tol=sqrt(tol))) for u0 in u0s]

# Number of clusters wrt. clustering parameter
data = reduce(hcat, shootings)
εs = exp10.(-10:0.1:0)
nb_clusters = [length(dbscan(data, ε).counts) for ε in εs]

with_logger(NullLogger()) do
    plot(scale=:log10, xlabel="ε", ylabel="Number of clusters", title="Number of clusters vs dbscan parameter ε", xticks=exp10.(-10:2:0))
    hline!([length(shootings) 2], c=[:black :red], alpha=0.5, label=["n = #samples" "n = 2"])
    plot!(εs, nb_clusters, c=:green3, label="Number of clusters")
    display(plot!())
end

# Clustering result
clusters = dbscan(data, 10sqrt(tol))
centroids = [mean(shootings[c.core_indices]) for c in clusters.clusters]
labels = ["dτ " * ifelse(c[3] > 0, ">", "<") * " 0" for c in centroids]
i = 1:length(centroids)
prop = clusters.counts / sum(clusters.counts)
Pavg = [average_power(s, vbp; tol=tol) for s in centroids]
plot(title="Clustering phase space wrt. limit cycles", xlabel="cluster index", xticks=(1:length(centroids), labels))
plot!(i, prop, seriestype=:bar, seriescolor=1:length(clusters.counts))
annotate!(i, prop / 2, [f"\%.0f(100 * p)%" for p in prop], :bottom)
annotate!(i, prop / 2 .+ 0.1, [f"P = \%.0f(Pavg) W" for Pavg in Pavg])
display(plot!())

# Recomputing the limit cycles with high precision
refined_shootings = [build_shooting(compute_limit_cycle(unpack_shooting(s)[1], vbp; tol=tol)) for s in centroids]

err_before = maximum(maximum.(abs, shooting_residuals.(centroids, Ref(vbp); tol=sqrt(tol))))
err_after = maximum(maximum.(abs, shooting_residuals.(refined_shootings, Ref(vbp); tol=tol / 10)))

@test err_after < err_before

@info f"
Error of the centroids using sqrt(tol) = \%.1e(sqrt(tol)) to find limit cycles (fast computation of many limit cycles)
\t\%.1e(err_before)
Error after recomputing the limit cycles using tol = \%.1e(tol) (accurate computation of a few limit cycles after clustering)
\t\%.1e(err_after)"


# In one call
limit_cycles = all_limit_cycles(vbp; αmin=αmin, αmax=αmax, vmax=vmax, N=100, tol=tol)

@test maximum(maximum.(abs, endpoint_residuals.(limit_cycles))) < 2tol

## Limit cycles by root finding

"""
Sometimes produces infinities and NaNs in the ODE solver, and throws an error, not usable."""
function _limit_cycle_rootfind(shooting, vbp; sense=default_sense(shooting), tol=tol)
    u0, _ = unpack_shooting(shooting)
    shooting1 = compute_limit_cycle(u0, vbp; tol=0.1)
    prob = NonlinearProblem((s, p) -> shooting_residuals(s, p; sense=sense, tol=tol / 100), shooting1, vbp)
    return solve(prob)
end
