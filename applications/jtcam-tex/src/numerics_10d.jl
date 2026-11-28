# Appendix A figures (10D model): drift of the raw 10D model and the effect of the
# manifold projection callback.
#
# Outputs: figs/trajectory_10D_short_long_cb.{pdf,png}, figs/residuals_10D_interp.{pdf,png}.
using KEEP: PointMass10 as PM10
using KEEP: PointMass4 as PM4
using KEEP: PointMassPara as PMP
using KEEP: Visualization as VIZ

using OrdinaryDiffEqTsit5: Tsit5

using Plots
using PyFormattedStrings
using LinearAlgebra: norm

include(joinpath(@__DIR__, "helpers.jl"))

function run_10d(; fast::Bool=false)
    #region setup_10d
    ## Setup
    p = PMP.build_para()
    τ0, dτ0 = 1e-10, 10
    u0_10 = PM10.init_u(τ0, dτ0, p)

    # Projection callback: its tolerance is independent of the integrator tolerance.
    cb_tol = 1e-10
    callback = PM10.build_manifold_projection(u0_10; save=true, tol=cb_tol)
    callback_no_save = PM10.build_manifold_projection(u0_10; save=false, tol=cb_tol)
    #endregion

    #region traj_10d
    # Physical plots at loose tolerance, to exaggerate the drift of the raw 10D model.
    tf_short, tf, tol = 20, 100, 1e-3
    sol10 = PM10.integrate(u0_10, tf, p, Tsit5(); save_everystep=true, tol)
    sol10_cb = PM10.integrate(u0_10, tf, p, Tsit5(); save_everystep=true, tol, callback=callback_no_save)

    ## 10D, Tsit5, tol=1e-3, tf = 20 vs 100 s vs callback
    points_per_second = 30
    fig_traj_short = VIZ.plot_trajectory_10D(sol10; tspan=(0, tf_short), points_per_second)
    fig_traj = VIZ.plot_trajectory_10D(sol10; points_per_second)
    fig_traj_cb = VIZ.plot_trajectory_10D(sol10_cb; points_per_second)
    params = (; xlabel="x (m)", ylabel="y (m)", zlabel="z (m)", legend=false, widen=false)
    plot!(fig_traj_short; title="tf = 20 s", params..., legend=:topleft)
    plot!(fig_traj; title="tf = 100 s", params...)
    plot!(fig_traj_cb; title="tf = 100 s + callback", params...)
    tf_fig = plot(fig_traj_short, fig_traj, fig_traj_cb, layout=(1, 3), size=SIZE_FULL)
    # one tick at each end of each axis, plus 0 when it is interior
    for i in 1:3
        plot!(tf_fig, subplot=i,
            xticks=edge_ticks(tf_fig[i], :x; include_zero=true),
            yticks=edge_ticks(tf_fig[i], :y; include_zero=true),
            zticks=edge_ticks(tf_fig[i], :z; include_zero=true))
    end
    files = save_figure(tf_fig, "trajectory_10D_short_long_cb")
    #endregion

    #region residuals_10d
    ## Interpolated residuals
    tf_long = fast ? 100 : 1000
    sol10 = PM10.integrate(u0_10, tf_long, p, Tsit5(); save_everystep=true, tol=1e-6)
    sol10_cb_t = PM10.integrate(u0_10, tf_long, p, Tsit5(); save_everystep=true, tol=1e-6, callback=callback_no_save)

    ts = logrange(1, sol10.prob.tspan[2], 501)
    res = [PM10.manifold_residuals!(similar(u0_10, 6), sol10(t), p) for t in ts]
    res_cb = [PM10.manifold_residuals!(similar(u0_10, 6), sol10_cb_t(t), p) for t in ts]
    ylim = exp10.((-9.5, 0.5))
    xticks_ = exp10.(-18:1:10)
    yticks_ = exp10.(-18:3:10)

    plot(; xscale=:log10, xticks=xticks_, xformatter=x -> "", ylim, yscale=:log10, yticks=yticks_, ylabel="Residual on position (m)", legend=:topleft)
    plot!(ts, norm.(getindex.(res, Ref(1:3))), label="Without callback")
    fig_pos = plot!(ts, norm.(getindex.(res_cb, Ref(1:3))), label="With callback")

    plot(; xscale=:log10, xticks=xticks_, xlabel="Time (s)", ylim, yscale=:log10, yticks=yticks_, ylabel="Residual on speed (m/s)")
    plot!(ts, norm.(getindex.(res, Ref(4:6))))
    fig_speed = plot!(ts, norm.(getindex.(res_cb, Ref(4:6))))
    plot(fig_pos, fig_speed, layout=(2, 1), size=SIZE_HALF .* (1, 1.5))
    plot!(bottommargin=3Plots.mm, leftmargin=6Plots.mm)
    append!(files, save_figure(current(), "residuals_10D_interp"))
    #endregion
    return files
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_10d()
end
