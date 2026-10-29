using Pkg
Pkg.activate(@__DIR__)

using KEEP: PointMass10 as PM10
using KEEP: PointMass4 as PM4
using KEEP: PointMassPara as PMP
using KEEP: Visualization as VIZ

using OrdinaryDiffEqVerner: Vern9
using OrdinaryDiffEqTsit5: Tsit5

using BenchmarkTools
using Plots
using PyFormattedStrings
using Logging

using LinearAlgebra: norm

include("jtcam_theme.jl")

## Setup
p = PMP.build_para()
τ0, dτ0 = 1e-10, 10
u0_10 = PM10.init_u(τ0, dτ0, p)
u0_4 = PM4.build_u(α=0, dα=0, dτ=dτ0; τ=τ0)
vbp = PMP.build_vbpara(p)

# Projection callback: its tolerance is independent of the integrator tolerance.
cb_tol = 1e-10
callback = PM10.build_manifold_projection(u0_10; save=true, tol=cb_tol)
callback_no_save = PM10.build_manifold_projection(u0_10; save=false, tol=cb_tol)

# Physical plots at loose tolerance, to exaggerate the drift of the raw 10D model.
tf_short, tf, tol = 20, 100, 1e-3
sol10 = PM10.integrate(u0_10, tf, p, Tsit5(); save_everystep=true, tol)
sol10_cb = PM10.integrate(u0_10, tf, p, Tsit5(); save_everystep=true, tol, callback=callback_no_save)

## 10D, Tsit5, tol=1e-3, tf = 20 vs 100 s vs callback
points_per_second = 30
fig_traj_short = VIZ.plot_trajectory_10D(sol10; tspan=(0, tf_short), points_per_second)
fig_traj = VIZ.plot_trajectory_10D(sol10; points_per_second)
fig_traj_cb = VIZ.plot_trajectory_10D(sol10_cb; points_per_second)
params = (; xlabel="x (m)", ylabel="y (m)", zlabel="z (m)", legend=false)
plot!(fig_traj_short; title="tf = 20 s", params..., legend=:topleft)
plot!(fig_traj; title="tf = 100 s", params...)
plot!(fig_traj_cb; title="tf = 100 s + callback", params...)
tf_fig = plot(fig_traj_short, fig_traj, fig_traj_cb, layout=(1, 3), size=SIZE_FULL)
# plot!(leftmargin=2Plots.mm)
savefig(tf_fig, joinpath(@__DIR__, "figs", "trajectory_10D_short_long_cb.pdf"))
display(tf_fig)

## Interpolated residuals
tf_long = 1000
sol10 = PM10.integrate(u0_10, tf_long, p, Tsit5(); save_everystep=true, tol=1e-6)
sol10_cb_f = PM10.integrate(u0_10, tf_long, p, Tsit5(); save_everystep=false, tol=1e-6, callback=callback)
sol10_cb_t = PM10.integrate(u0_10, tf_long, p, Tsit5(); save_everystep=true, tol=1e-6, callback=callback_no_save)

ts = logrange(1, sol10.prob.tspan[2], 501)
res = [PM10.manifold_residuals!(similar(u0_10, 6), sol10(t), p) for t in ts]
res_cb = [PM10.manifold_residuals!(similar(u0_10, 6), sol10_cb_t(t), p) for t in ts]
ylim = exp10.((-9.5, 0.5))
xticks_ = exp10.(-18:1:10)
yticks_ = exp10.(-18:3:10)

plot(; xscale=:log10, xticks=xticks_, xformatter=x->"", ylim, yscale=:log10, yticks=yticks_, ylabel="Residual on position (m)", legend=:topleft)
plot!(ts, norm.(getindex.(res, Ref(1:3))), label="Without callback")
fig_pos = plot!(ts, norm.(getindex.(res_cb, Ref(1:3))), label="With callback")

plot(; xscale=:log10, xticks=xticks_, xlabel="Time (s)", ylim, yscale=:log10, yticks=yticks_, ylabel="Residual on speed (m/s)")
plot!(ts, norm.(getindex.(res, Ref(4:6))))
fig_speed = plot!(ts, norm.(getindex.(res_cb, Ref(4:6))))
plot(fig_pos, fig_speed, layout=(2, 1), size=SIZE_HALF .* (1, 1.5))
plot!(bottommargin=3Plots.mm, leftmargin=6Plots.mm)
display(current())
savefig(current(), joinpath(@__DIR__, "figs", "residuals_10D_interp.pdf"))


## --- Model equivalence at maximum accuracy (Alg = Vern9, abstol = reltol = 1e-10) ---
# Short window: the raw 10D model's constraint drift grows with time, keep it small.
tf_cmp = 100.0
tc = range(0, tf_cmp, length=1000)
pos10(ss) = [PM10.compute_pos1(q, p) for q in ss.(tc, idxs=1:5)]
pos4(ss) = [PM4.compute_OK(PM4.compute_Rτ(q, vbp), vbp) for q in ss.(tc, idxs=1:2)]
dist(a, b) = norm.(a .- b)

s10 = PM10.integrate(u0_10, tf_cmp, p, Vern9(); tol=1e-10, save_everystep=true)
s10cb = PM10.integrate(u0_10, tf_cmp, p, Vern9(); tol=1e-10, save_everystep=true, callback)
s4 = PM4.integrate(u0_4, tf_cmp, vbp, Vern9(); tol=1e-10, save_everystep=true)
s4f = PM4.integrate(u0_4, tf_cmp, vbp, Vern9(); tol=1e-12, save_everystep=true)  # numerical floor

p10_t, p10cb_t, p4_t, p4f_t = pos10(s10), pos10(s10cb), pos4(s4), pos4(s4f)
floor_d = dist(p4_t, p4f_t)

# Pairwise kite-position difference (m) against the numerical floor
v1 = plot(title="Pairwise kite-position difference (Vern9, 1e-10)",
    xlabel="Time (s)", ylabel="Position difference (m)", yscale=:log,
    ylims=(1e-13, 1e-9), left_margin=8Plots.mm, bottom_margin=6Plots.mm)
plot!(v1, tc, dist(p10_t, p4_t), label="10D vs 4D")
plot!(v1, tc, dist(p10cb_t, p4_t), label="10D + callback vs 4D")
plot!(v1, tc, dist(p10_t, p10cb_t), label="10D vs 10D + callback")
plot!(v1, tc, floor_d, label="Numerical floor", ls=:dash, c=:black)

# Common significant digits of the kite position, relative to the line length r
names = ["10D", "10D + callback", "4D"]
digits(a, b) = -log10(maximum(dist(a, b)) / p.r)
D = [a === b ? NaN : digits(a, b) for a in (p10_t, p10cb_t, p4_t), b in (p10_t, p10cb_t, p4_t)]
v3 = heatmap(D, xticks=(1:3, names), yticks=(1:3, names), xmirror=true, yflip=true,
    color=:Purples, aspect_ratio=1, title="Common significant digits",
    colorbar_title="matching digits")
for i in 1:3, j in 1:3
    isnan(D[i, j]) || annotate!(v3, i, j, text(round(D[i, j], digits=1), 10, :black))
end

cmp = plot(v1, v3, layout=(1, 2), size=SIZE_FULL)
savefig(cmp, joinpath(@__DIR__, "figs", "_comparison.pdf"))
display(cmp)

## Models converge to the same reference trajectory
pos10(sol) = [PM10.compute_pos1(u[1:5], sol.prob.p) for u in sol.u]
pos4(sol) = [PM4.compute_OK(PM4.compute_Rτ(u[1:2], sol.prob.p), sol.prob.p) for u in sol.u]

saveat = range(0, tf_short, 1001)  # sample times
tols = exp10.(-2:-1:-12)
ref = PM4.integrate(u0_4, tf_short, vbp, Vern9(); tol=1e-14, saveat)
pref = pos4(ref)

err10 = [maximum(dist(pos10(PM10.integrate(u0_10, tf_short, p, Vern9(); tol, saveat)), pref)) for tol in tols]
err10cb = [maximum(dist(pos10(PM10.integrate(u0_10, tf_short, p, Vern9(); tol, saveat, callback=callback_no_save)), pref)) for tol in tols]
err4 = [maximum(dist(pos4(PM4.integrate(u0_4, tf_short, vbp, Vern9(); tol, saveat)), pref)) for tol in tols]

plot(xlabel="Tolerance (relative and absolute)", ylabel="Norm of error (m)", legend=:bottomright)
plot!(tols, err10, marker=:circle, label="10D")
plot!(tols, err10cb, marker=:square, label="10D + callback")
plot!(tols, err4, marker=:diamond, label="4D")
plot!(xscale=:log10, xticks=exp10.(-15:3:15), yscale=:log10, yticks=exp10.(-15:3:15))
fig = plot!(xflip=true, ymirror=true, legend=:topright, size=SIZE_HALF)
plot!(bottommargin=4Plots.mm, rightmargin=6Plots.mm)
savefig(fig, joinpath(@__DIR__, "figs", "three_models_convergence.pdf"))
display(fig)

#=
Max(OK_i(t) - OK_ref(t)), t in [0, t_f]

We compare the error of all three implementations (10D, 10D + callback, and 4D) across a range of tolerances (10^-2 to 10^-12) against a reference run (4D, tolerance at 10^-14). As the tolerance gets tighter, all runs converge to the reference trajectory, showing they are all equivalent.

Note: The horizontal axis is reversed so that tighter integration tolerances progress to the right.
=#
## --- Timing at the production tolerance -----------------------------------------
btime4 = @belapsed PM4.integrate(u0_4, tf, vbp; save_everystep=true, tol) seconds = 1
btime = @belapsed PM10.integrate(u0_10, tf, p; save_everystep=true, tol) seconds = 1
btime_cb = @belapsed PM10.integrate(u0_10, tf, p; save_everystep=true, tol, callback) seconds = 1

ratio_cb = btime_cb / btime
speedup = btime / btime4
speedup_cb = btime_cb / btime4
@info f"Timing
reference (10D, no callback) = {1e3*btime:.2g} ms
reference + callback = {1e3*btime_cb:.2g} ms
proposed model (4D) = {1e3*btime4:.2g} ms

reference to reference+callback = {ratio_cb:.2f}x slowdown
reference to proposed model = {speedup:.2f}x speedup
reference+callback to proposed model = {speedup_cb:.2f}x speedup"

@info f"Conclusion
The overdetermined system, which requires carefull initialization (cf. algebraic constraints) experiences drifting when numerically integrated, which is not suitable for finer analysis.
By adding a manifold projection callback, we effectively nullify the drift -- at the cost of performance (about {ratio_cb:.1f} times slower), and we would still be working with an inherently inconsistent state (residuals are not exactly 0)
The proposed model (4D, formulated in α and τ only) carries no algebraic constraints, stays consistent by construction, and reproduces the 10D configuration while running about {speedup:.1f} times faster than the 10D reference and {speedup_cb:.1f} times faster than the same reference with the projection callback."
