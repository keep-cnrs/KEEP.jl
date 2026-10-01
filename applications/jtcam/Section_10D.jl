using KEEP: PointMass10 as PM10
using KEEP: PointMass4 as PM4
using KEEP: PointMassPara as PMP
using KEEP: Visualization as VIZ

using OrdinaryDiffEqVerner: Vern9

using BenchmarkTools
using Plots
using PyFormattedStrings
using Logging

using LinearAlgebra: norm

p = PMP.build_para()
τ0, dτ0 = 1e-10, 10
u0_10 = PM10.init_u(τ0, dτ0, p)
u0_4 = PM4.build_u(; α=0, dα=0, dτ=dτ0, τ=τ0)
vbp = PMP.build_vbpara(p)

# Projection callback: its tolerance is independent of the integrator tolerance.
callback = PM10.build_manifold_projection(u0_10; save=true, tol=1e-12)

# Physical plots at loose tolerance, to exaggerate the drift of the raw 10D model.
tf, tol = 80, 1e-3
sol10 = PM10.integrate(u0_10, tf, p; save_everystep=true, tol)
sol10_cb = PM10.integrate(u0_10, tf, p; save_everystep=false, tol, callback)

# 10D, Tsit5, tol=1e-3, tf = 20 vs 80 s
p20 = VIZ.plot_trajectory_10D(sol10; tspan=(0, 20))
p80 = VIZ.plot_trajectory_10D(sol10; tspan=(0, 80))
plot!(
    p20;
    title="10D, tf=20s, Tsit5, tol=1e-3",
    xlabel="x (m)",
    ylabel="y (m)",
    zlabel="z (m)",
    legend=:topleft,
)
plot!(
    p80;
    title="10D, tf=80s, Tsit5, tol=1e-3",
    xlabel="x (m)",
    ylabel="y (m)",
    zlabel="z (m)",
    legend=false,
)
tf_fig = plot(p20, p80; layout=(1, 2), size=(1200, 600))
savefig(tf_fig, joinpath(@__DIR__, "figs", "trajectory_10D_tf.png"))
display(tf_fig)

# 10D, Tsit5, tol=1e-3, tf = 80 s, no callback vs. callback
p10 = VIZ.plot_trajectory_10D(sol10)
p10_cb = VIZ.plot_trajectory_10D(sol10_cb)
plot!(
    p10;
    title="10D, tf=80s, Tsit5, tol=1e-3",
    xlabel="x (m)",
    ylabel="y (m)",
    zlabel="z (m)",
)
plot!(
    p10_cb;
    title="10D + callback, tf=80s, Tsit5, tol=1e-3",
    xlabel="x (m)",
    ylabel="y (m)",
    zlabel="z (m)",
)
plot!(p10; legend=:topleft)
plot!(p10_cb; legend=false)
fig = plot(p10, p10_cb; layout=(1, 2), size=(1200, 600))
savefig(fig, joinpath(@__DIR__, "figs", "trajectory_10D_callback.png"))
display(fig)

# Single-figure alternative: kite trajectory solid until 20s then dotted, arm/lines always solid.
# Both segments share t=20 exactly, so the transition is seamless and 0-20s is never plotted twice.
p_nc = VIZ.plot_trajectory_10D(sol10; tspan=(0, 20))
VIZ.plot_trajectory_10D(sol10; tspan=(20, 80), kite_ls=:dot, new=false)
plot!(
    p_nc;
    title="10D: solid 0-20s, dotted 20-80s",
    xlabel="x (m)",
    ylabel="y (m)",
    zlabel="z (m)",
    legend=:topleft,
)
p_cb = VIZ.plot_trajectory_10D(sol10_cb)
plot!(
    p_cb;
    title="10D + callback, tf=80s, Tsit5, tol=1e-3",
    xlabel="x (m)",
    ylabel="y (m)",
    zlabel="z (m)",
    legend=false,
)
split_fig = plot(p_nc, p_cb; layout=(1, 2), size=(1200, 600))
savefig(split_fig, joinpath(@__DIR__, "figs", "trajectory_10D_tf_callback.png"))
display(split_fig)

plot(; title="Configuration residuals without vs. with callback")
plot!(
    sol10.t,
    [norm(PM10.manifold_residuals!(similar(u0_10, 6), u, p)) for u in sol10.u];
    label="No callback",
)
plot!(
    sol10_cb.t,
    [norm(PM10.manifold_residuals!(similar(u0_10, 6), u, p)) for u in sol10_cb.u];
    label="Callback",
)
hline!([tol]; label="Tolerance", c=:black)
res_fig = plot!(;
    xlabel="Time (s)",
    ylabel="Residual norm (m, m/s)",
    yscale=:log,
    yticks=exp10.(-15:3:3),
    legend=:right,
)
savefig(res_fig, joinpath(@__DIR__, "figs", "residuals_10D.png"))
display(res_fig)

# --- Model equivalence at maximum accuracy (Alg = Vern9, abstol = reltol = 1e-10) ---
# Short window: the raw 10D model's constraint drift grows with time, keep it small.
tf_cmp = 3.0
tc = range(0, tf_cmp; length=1000)
pos10(ss) = [PM10.compute_pos1(q, p) for q in ss.(tc; idxs=1:5)]
pos4(ss) = [PM4.compute_OK(PM4.compute_Rτ(q, vbp), vbp) for q in ss.(tc; idxs=1:2)]
dist(a, b) = norm.(a .- b)

s10 = PM10.integrate(u0_10, tf_cmp, p, Vern9(); tol=1e-10, save_everystep=true)
s10cb = PM10.integrate(u0_10, tf_cmp, p, Vern9(); tol=1e-10, save_everystep=true, callback)
s4 = PM4.integrate(u0_4, tf_cmp, vbp, Vern9(); tol=1e-10, save_everystep=true)
s4f = PM4.integrate(u0_4, tf_cmp, vbp, Vern9(); tol=1e-12, save_everystep=true)  # numerical floor

p10_t, p10cb_t, p4_t, p4f_t = pos10(s10), pos10(s10cb), pos4(s4), pos4(s4f)
floor_d = dist(p4_t, p4f_t)

# Pairwise kite-position difference (m) against the numerical floor
v1 = plot(;
    title="Pairwise kite-position difference (Vern9, 1e-10)",
    xlabel="Time (s)",
    ylabel="Position difference (m)",
    yscale=:log,
    ylims=(1e-13, 1e-9),
    left_margin=8Plots.mm,
    bottom_margin=6Plots.mm,
)
plot!(v1, tc, dist(p10_t, p4_t); label="10D vs 4D")
plot!(v1, tc, dist(p10cb_t, p4_t); label="10D + callback vs 4D")
plot!(v1, tc, dist(p10_t, p10cb_t); label="10D vs 10D + callback")
plot!(v1, tc, floor_d; label="Numerical floor", ls=:dash, c=:black)

# Common significant digits of the kite position, relative to the line length r
names = ["10D", "10D + callback", "4D"]
digits(a, b) = -log10(maximum(dist(a, b)) / p.r)
D = [
    a === b ? NaN : digits(a, b) for
    a in (p10_t, p10cb_t, p4_t), b in (p10_t, p10cb_t, p4_t)
]
v3 = heatmap(
    D;
    xticks=(1:3, names),
    yticks=(1:3, names),
    xmirror=true,
    yflip=true,
    color=:Purples,
    aspect_ratio=1,
    title="Common significant digits",
    colorbar_title="matching digits",
)
for i in 1:3, j in 1:3
    isnan(D[i, j]) || annotate!(v3, i, j, text(round(D[i, j]; digits=1), 10, :black))
end

cmp = plot(v1, v3; layout=(1, 2), size=(1100, 500))
savefig(cmp, joinpath(@__DIR__, "figs", "_comparison.png"))
display(cmp)

# --- Error vs integration tolerance, against a very tight 4D reference ---
tols = [1e-3, 1e-6, 1e-9, 1e-12]
ref = PM4.integrate(u0_4, tf_cmp, vbp, Vern9(); tol=1e-13, save_everystep=true)
pref = pos4(ref)

err10 = [
    maximum(
        dist(
            pos10(PM10.integrate(u0_10, tf_cmp, p, Vern9(); tol=t, save_everystep=true)),
            pref,
        ),
    ) for t in tols
]
err10cb = [
    maximum(
        dist(
            pos10(
                PM10.integrate(
                    u0_10, tf_cmp, p, Vern9(); tol=t, save_everystep=true, callback
                ),
            ),
            pref,
        ),
    ) for t in tols
]
err4 = [
    maximum(
        dist(
            pos4(PM4.integrate(u0_4, tf_cmp, vbp, Vern9(); tol=t, save_everystep=true)),
            pref,
        ),
    ) for t in tols
]

conv = plot(;
    title="Kite-position error vs integration tolerance",
    xlabel="Tolerance (abstol = reltol)",
    ylabel="Max kite-position error (m)",
    xscale=:log,
    yscale=:log,
    legend=:bottomright,
)
plot!(conv, tols, err10; marker=:circle, label="10D")
plot!(conv, tols, err10cb; marker=:square, label="10D + callback")
plot!(conv, tols, err4; marker=:diamond, label="4D")
savefig(conv, joinpath(@__DIR__, "figs", "convergence.png"))
display(conv)

# --- Timing at the production tolerance -----------------------------------------
btime4 = @belapsed PM4.integrate(u0_4, tf, vbp; save_everystep=true, tol) seconds = 0.1
btime = @belapsed PM10.integrate(u0_10, tf, p; save_everystep=true, tol) seconds=0.1
btime_cb = @belapsed PM10.integrate(u0_10, tf, p; save_everystep=true, tol, callback) seconds=0.1

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
