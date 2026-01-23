using Setfield
using Plots
using StrFormat
using StaticArrays
using NonlinearSolve: NewtonRaphson, ReturnCode
using ADTypes: AutoForwardDiff
using LinearAlgebra: norm
using Test

import KEEP.PointMass10 as PM10
import KEEP.PointMass4 as PM4
using KEEP.PointMassPara
using KEEP.TorqueFunction
import KEEP.Visualisation as VIS

include("utils.jl")
default(lw=3, formatter=:plain, label="")

τ0, dτ0 = 1e-1, 0
tf = 60.0
tf_small = 10  # for short term plots
step = 0.02
tol = 5e-4

t = 0:step:tf
t_small = 0:step:tf_small

## plot position 10d (20 s et 80 s)
p = build_para()
u0 = PM10.init_u(τ0, dτ0, p)
sol10 = PM10.integrate(u0, tf, p; save_everystep=true, tol=tol)
cb = PM10.build_manifold_projection(u0)
sol10_cb = PM10.integrate(u0, tf, p; save_everystep=false, tol=tol, callback=cb)
@test sol10.retcode == ReturnCode.Success
@test sol10_cb.retcode == ReturnCode.Success

VIS.plot_trajectory_10D(sol10; tspan=tf_small)
fig1 = plot!(title=f"10D model, t=0..\%d(tf_small)s")

VIS.plot_trajectory_10D(sol10; tspan=tf)
fig2 = plot!(title=f"10D model, t=0..\%d(tf)s")

VIS.plot_trajectory_10D(sol10_cb; tspan=tf)
fig3 = plot!(title=f"10D model with callback, t=0..\%d(tf)s")

plots = [fig1 fig2 fig3]
display(plot(plots..., size=500 .* size(plots'), layout=size(plots)))

## Residuals
plot(title="Residuals with vs. without callback")
plot!(sol10.t, [norm(PM10.manifold_residuals!(similar(u0, 6), u, p)) for u in sol10.u], label="Without")
plot!(sol10_cb.t, [norm(PM10.manifold_residuals!(similar(u0, 6), u, p)) for u in sol10_cb.u], label="With")
display(plot!())

## Plot position 4d (dimensionné et normalisé)
vbp = build_vbpara(p)
vbp_normed = normalize_vbpara(vbp)
sol4 = PM4.integrate(SA[0.0, τ0, 0.0, dτ0, 0.0], tf, vbp; save_everystep=true)

_, _, T = lmt(p)
sol4_normed = PM4.integrate(SA[0.0, τ0, 0.0, dτ0*T, 0.0], tf / T, vbp_normed; save_everystep=true)

VIS.plot_trajectory_4D(sol4; tspan=tf)
fig1 = plot!(title=f"4D model, t=0..\%d(tf)s")

VIS.plot_trajectory_4D(sol4_normed; tspan=tf / T)
fig2 = plot!(title=f"Normalised 4D model, t=0..\%d(tf)s")

plots = [fig1 fig2]
display(plot(plots..., size=500 .* size(plots'), layout=size(plots)))

## energy generation plot for the three of them [10D, 4D, 4DVB normalised]
VIS.plot_avg_power_10D(sol10; tspan=tf_small)
fig1 = plot!(title="10D model")

VIS.plot_avg_power_4D(sol4; tspan=tf_small)
fig2 = plot!(title="4D model")

VIS.plot_avg_power_4D(sol4_normed; tspan=tf_small ./ T)  # Same shape but different numbers since we don't have the reference magnitudes anymore
fig3 = plot!(title="Normalised 4D model")

plots = [fig1 fig2 fig3]
display(plot(plots..., size=(500, 350) .* size(plots'), layout=size(plots)))

# animate the position of the kite (only 4D)
if false
    fps = 30
    anim = VIS.animate_trajectory_4D(sol4; tspan=tf_small, fps=fps)
    display(gif(anim; fps=fps))
    Plots.apng(anim, fps=fps)
    # display(mov(anim; fps=fps))
    # display(mp4(anim; fps=fps))
    # webm(anim; fps=fps)  # Only one that supports 60 fps
end

## For presentation
VIS.plot_avg_power_4D(sol4; tspan=tf_small)
plot!(title="Power generation over time", size=(600, 300))
# savefig("docs/media/power_generation.pdf")

VIS.plot_trajectory_4D(sol4; tspan=tf)
plot!(size=(400, 300))
# savefig("docs/media/trajectory_4D.pdf")

plot!(xlim=(-vbp.l, 3vbp.l), ylim=(-2vbp.l, 2vbp.l), camera=(30, 45))
# savefig("docs/media/trajectory_4D_arm_zoom.pdf")



if false
    fps = 30
    anim = VIS.animate_trajectory_4D(@set sol4.prob.p.r = 6; tspan=(0, tf_small), fps=fps, end_camera=(30, 30), figsize=(400, 400))
    display(gif(anim; fps=fps))
    # gif(anim, "docs/media/trajectory_animated.gif"; fps=fps)
end
## END for presentation


using SplitApplyCombine

τ0 = 0.1
dτ0 = 10
u0 = SA[0.0, τ0, 0.0, dτ0, 0.0]
sol = PM4.integrate(u0, tf, vbp; save_everystep=true)
t = 0:0.01:2
K = invert([PM4.compute_OK(PM4.compute_Rτ(sol(t), vbp), vbp) for t in t])
plot(t, K, label=["x" "y" "z"])