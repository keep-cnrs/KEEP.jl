# Section_10D.jl — compute-only library
#
# Numerical results and plot-data for the 10D-vs-4D validation figures.
# All plotting, styling and notes live in paper.jl; this file only computes.

using KEEP: PointMass10 as PM10
using KEEP: PointMass4 as PM4
using KEEP: PointMassPara as PMP

using OrdinaryDiffEqVerner: Vern9
using OrdinaryDiffEqTsit5: Tsit5

using BenchmarkTools
using LinearAlgebra: norm

"""
    compute_10D()

Compute every numerical result and plot-data array needed by the 10D validation figures.

# Returned fields
- `sol10`: 10D Tsit5 solution at tf = 100 s, tol = 1e-3, no callback (trajectory panels).
- `sol10_cb`: same as `sol10`, with the manifold-projection callback (no saving).
- `tf_short`: short trajectory horizon, 20 s (first trajectory panel).
- `tf`: trajectory horizon, 100 s.
- `points_per_second`: sampling rate passed to the trajectory plotter.
- `ts`: log-spaced sample times for the interpolated-residual figure.
- `res`: position/speed residual vectors on `ts`, raw 10D model.
- `res_cb`: position/speed residual vectors on `ts`, with projection callback.
- `ylim`: y-axis limits shared by the two residual panels.
- `xticks_`: x tick locations for the residual panels.
- `yticks_`: y tick locations for the residual panels.
- `tols`: integration tolerances for the convergence figure.
- `err10`: max position error of the 10D model at each tolerance.
- `err10cb`: max position error of the 10D + callback model at each tolerance.
- `err4`: max position error of the 4D model at each tolerance.
- `btime4`: elapsed time of the 4D production run (s).
- `btime`: elapsed time of the 10D production run (s).
- `btime_cb`: elapsed time of the 10D + callback production run (s).
- `ratio_cb`: slowdown from adding the projection callback to 10D.
- `speedup`: speedup of the 4D model over raw 10D.
- `speedup_cb`: speedup of the 4D model over 10D with the projection callback.
"""
function compute_10D()
    p = PMP.build_para()
    τ0, dτ0 = 1e-10, 10
    u0_10 = PM10.init_u(τ0, dτ0, p)
    u0_4 = PM4.build_u(; α=0, dα=0, dτ=dτ0, τ=τ0)
    vbp = PMP.build_vbpara(p)

    # Projection callback: its tolerance is independent of the integrator tolerance.
    cb_tol = 1e-10
    callback = PM10.build_manifold_projection(u0_10; save=true, tol=cb_tol)
    callback_no_save = PM10.build_manifold_projection(u0_10; save=false, tol=cb_tol)

    # Physical plots at loose tolerance, to exaggerate the drift of the raw 10D model.
    tf_short, tf, tol = 20, 100, 1e-3
    sol10 = PM10.integrate(u0_10, tf, p, Tsit5(); save_everystep=true, tol)
    sol10_cb = PM10.integrate(
        u0_10, tf, p, Tsit5(); save_everystep=true, tol, callback=callback_no_save
    )
    points_per_second = 30

    # Interpolated residuals
    tf_long = 1000
    sol10_long = PM10.integrate(u0_10, tf_long, p, Tsit5(); save_everystep=true, tol=1e-6)
    sol10_cb_f = PM10.integrate(
        u0_10, tf_long, p, Tsit5(); save_everystep=false, tol=1e-6, callback=callback
    )
    sol10_cb_t = PM10.integrate(
        u0_10, tf_long, p, Tsit5(); save_everystep=true, tol=1e-6, callback=callback_no_save
    )

    ts = logrange(1, sol10_long.prob.tspan[2], 501)
    res = [PM10.manifold_residuals!(similar(u0_10, 6), sol10_long(t), p) for t in ts]
    res_cb = [PM10.manifold_residuals!(similar(u0_10, 6), sol10_cb_t(t), p) for t in ts]
    ylim = exp10.((-9.5, 0.5))
    xticks_ = exp10.(-18:1:10)
    yticks_ = exp10.(-18:3:10)

    # Models converge to the same reference trajectory
    pos10(sol) = [PM10.compute_pos1(u[1:5], sol.prob.p) for u in sol.u]
    pos4(sol) =
        [PM4.compute_OK(PM4.compute_Rτ(u[1:2], sol.prob.p), sol.prob.p) for u in sol.u]
    dist(a, b) = norm.(a .- b)

    saveat = range(0, tf_short, 1001)  # sample times
    tols = exp10.(-2:-1:-12)
    ref = PM4.integrate(u0_4, tf_short, vbp, Vern9(); tol=1e-14, saveat)
    pref = pos4(ref)

    err10 = [
        maximum(
            dist(pos10(PM10.integrate(u0_10, tf_short, p, Vern9(); tol, saveat)), pref)
        ) for tol in tols
    ]
    err10cb = [
        maximum(
            dist(
                pos10(
                    PM10.integrate(
                        u0_10, tf_short, p, Vern9(); tol, saveat, callback=callback_no_save
                    ),
                ),
                pref,
            ),
        ) for tol in tols
    ]
    err4 = [
        maximum(dist(pos4(PM4.integrate(u0_4, tf_short, vbp, Vern9(); tol, saveat)), pref))
        for tol in tols
    ]

    #=
    Max(OK_i(t) - OK_ref(t)), t in [0, t_f]

    We compare the error of all three implementations (10D, 10D + callback, and 4D) across a range of tolerances (10^-2 to 10^-12) against a reference run (4D, tolerance at 10^-14). As the tolerance gets tighter, all runs converge to the reference trajectory, showing they are all equivalent.

    Note: The horizontal axis is reversed so that tighter integration tolerances progress to the right.
    =#

    ## --- Timing at the production tolerance ---
    btime4 = @belapsed PM4.integrate($u0_4, $tf, $vbp; save_everystep=true, tol=$tol) seconds =
        1
    btime = @belapsed PM10.integrate($u0_10, $tf, $p; save_everystep=true, tol=$tol) seconds =
        1
    btime_cb = @belapsed PM10.integrate(
        $u0_10, $tf, $p; save_everystep=true, tol=$tol, callback=$callback
    ) seconds = 1

    ratio_cb = btime_cb / btime
    speedup = btime / btime4
    speedup_cb = btime_cb / btime4

    # Notes carried into paper.jl (values: 1e3*btime, 1e3*btime_cb, 1e3*btime4 ms):
    # reference to reference+callback = ratio_cb x slowdown
    # reference to proposed model = speedup x speedup
    # reference+callback to proposed model = speedup_cb x speedup
    # The overdetermined system, which requires carefull initialization (cf. algebraic constraints), experiences
    # drifting when numerically integrated, which is not suitable for finer analysis. By adding a manifold
    # projection callback, we effectively nullify the drift -- at the cost of performance (about ratio_cb times
    # slower), and we would still be working with an inherently inconsistent state (residuals are not exactly 0).
    # The proposed model (4D, formulated in α and τ only) carries no algebraic constraints, stays consistent by
    # construction, and reproduces the 10D configuration while running about speedup times faster than the 10D
    # reference and speedup_cb times faster than the same reference with the projection callback.

    return (;
        sol10,
        sol10_cb,
        tf_short,
        tf,
        points_per_second,
        ts,
        res,
        res_cb,
        ylim,
        xticks_,
        yticks_,
        tols,
        err10,
        err10cb,
        err4,
        btime4,
        btime,
        btime_cb,
        ratio_cb,
        speedup,
        speedup_cb,
    )
end
