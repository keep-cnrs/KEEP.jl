using QuasiMonteCarlo: sample, HaltonSample, SobolSample
using Plots
using SplitApplyCombine
using Setfield

import KEEP.SteadyState as SS
import KEEP.LimitCycle as LC
using KEEP.PointMassPara

N = 1001
xmax = 5
Xs = eachcol(sample(N, [-xmax, -xmax], [xmax, xmax], HaltonSample()))
Xs = eachcol(sample(N, [-xmax, -xmax], [xmax, xmax], SobolSample()))
Xs = [2xmax .* rand(2) .- xmax for i in 1:N]
ys = [round(Int, sqrt(sum(x.^2))) for x in Xs]
scatter(invert(Xs)..., group=ys, c=ys, ratio=:equal, markerstrokewidth=0)


vbp = build_vbpara()
syms = [:r, :Δφ]
μ = 3
lb = vbp[syms] ./ μ |> collect
ub = vbp[syms] .* μ |> collect
lb = [10, .1]
ub = [100, 1.5]

# lb = [15, .4]
# ub = [35, .55]


Ps = sample(1000, lb, ub, HaltonSample()) |> eachcol
@time ys = [length(SS.all_steady_states(@set vbp[syms] = p)) for p in Ps]

inds = copy(ys)
unique_vals = sort(unique(ys))
for (i, y) in enumerate(unique_vals) inds[ys .== y] .= i end
cmap = palette(:tab10)

scatter(invert(Ps)..., color=cmap[inds],
xlabel=syms[1], ylabel=syms[2], title="Number of steady states", markerstrokewidth=0, label="")
scatter!(length(unique_vals), markerstrokewidth=0, color=cmap[1:length(unique_vals)]', label=unique_vals')

syms = [:r, :I_eq, :torque_slope, :Δθ, :Δφ]
μ = 3
lb = vbp[syms] ./ μ
ub = vbp[syms] .* μ

Ps = sample(100, lb, ub, HaltonSample()) |> eachcol
@time ys = [
    try
        length(LC.all_limit_cycles(@set vbp[syms] = p))
    catch e
        -1
    end
    for p in Ps]
unique(ys)