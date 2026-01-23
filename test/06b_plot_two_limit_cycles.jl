using Logging
using Plots
using StaticArrays
using StrFormat

import KEEP.PointMass4 as PM4
using KEEP.PointMassPara
using KEEP.LimitCycle
import KEEP: TAU0
using KEEP.TorqueFunction
using KEEP.Visualisation: plot_trajectory_4D

p = build_vbpara()

tol = 1e-13
lc_neg = compute_limit_cycle(SA[0.0, TAU0, 0.0, -10, 0.0], p; tol=tol, save_everystep=true)
lc_pos = compute_limit_cycle(SA[0.0, TAU0, 0.0, 10, 0.0], p; tol=tol, save_everystep=true)

@info f"Average power for the negative cycle (good one): \%.2f(lc_neg.u[end][end] / lc_neg.t[end])\n" * f"Average power for the positive cycle (bad one) : \%.2f(lc_pos.u[end][end] / lc_pos.t[end])"

t_pos = range(lc_pos.t[1], lc_pos.t[end], length=101)
t_neg = range(lc_neg.t[1], lc_neg.t[end], length=101)

α_pos = lc_pos.(t_pos; idxs=1)
τ_pos = rem.(lc_pos.(t_pos; idxs=2), 2π, RoundDown)
i_pos = sortperm(τ_pos)
α_neg = lc_neg.(t_neg; idxs=1)
τ_neg = rem.(lc_neg.(t_neg; idxs=2), 2π, RoundDown)
i_neg = sortperm(τ_neg)

fig = plot(title="α(τ)", xlabel="τ", ylabel="α", legend=:topright)
plot!(2π .- τ_neg[i_neg], -α_neg[i_neg], label="cycle OK")
plot!(τ_pos[i_pos], α_pos[i_pos], label="cycle KO")
display(fig)

dα_pos = lc_pos.(t_pos; idxs=3)
dα_neg = lc_neg.(t_neg; idxs=3)

fig = plot(title="dα(τ)", xlabel="τ", ylabel="dα", legend=:topright)
plot!(2π .- τ_neg[i_neg], -dα_neg[i_neg], label="cycle OK")
plot!(τ_pos[i_pos], dα_pos[i_pos], label="cycle KO")
display(fig)


plot_trajectory_4D(lc_pos, points_per_second=300)
plot!(title="Bad limit cycle")
display(plot!())

plot_trajectory_4D(lc_neg, points_per_second=300)
plot!(title="Good limit cycle")
display(plot!())


## For presentation
fig = plot(title=raw"2 limit cycles", xlabel=raw"$τ$", ylabel=raw"$α(τ)$", legend=:topright, size=(300, 200))
plot!(2π .- τ_neg[i_neg], -α_neg[i_neg], lw=3, label=raw"cycle OK $←$")
plot!(τ_pos[i_pos], α_pos[i_pos], lw=3, label=raw"cycle KO $→$")
display(fig)
# savefig(fig, "docs/media/two_limit_cycles_alpha.pdf")
