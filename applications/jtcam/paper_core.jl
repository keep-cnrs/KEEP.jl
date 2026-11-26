# paper_core.jl — shared numerical helpers for the JTCAM paper.
#
# Single source of truth for the helpers that were previously copied between
# paper.jl and the Section_*.jl files. Included by paper.jl at top level; the
# Section files only *compute* results and never plot.

using KEEP: PointMass4 as PM4, PointMassPara as PMP, LimitCycle as LC
using KEEP: DEFAULT_TOLERANCE, TAU0
using OrdinaryDiffEqVerner: Vern9, ODEProblem
using BifurcationKit
const BK = BifurcationKit
using QuasiMonteCarlo: sample, HaltonSample
using LinearAlgebra: norm, eigvals, ×
using ForwardDiff: ForwardDiff

# Work done by the generator is measured from the arm velocity; the absolute
# velocity of the kite is needed for acceleration and curvature diagnostics.
function kite_velocity(u, p)
    pos(q_) = PM4.compute_OK(PM4.compute_Rτ(q_, p), p)
    return ForwardDiff.jacobian(pos, u[1:2]) * u[3:4]
end

# Metrics of one limit cycle. `nothing` means no admissible orbit was found:
# return NaN rather than crashing so parameter sweeps can keep marching.
function cycle_metrics(sol; n=1001)
    sol === nothing && return (;
        P=NaN,
        period=NaN,
        Tmin=NaN,
        Tmax=NaN,
        αmax=NaN,
        amax=NaN,
        vmax=NaN,
        κmax=NaN,
        Rmin=NaN,
    )
    ts = range(sol.t[1], sol.t[end], n)
    us = sol.(ts)
    p = sol.prob.p
    tension = PM4.compute_line_tension.(us, Ref(p))
    ok = isfinite.(tension)
    vel = kite_velocity.(us, Ref(p))
    acc = PM4.kite_acceleration.(us, Ref(p))
    κ = [norm(v × a) / max(norm(v)^3, 1e-12) for (v, a) in zip(vel, acc)]
    return (;
        P=sol.u[end][5] / sol.t[end],
        period=sol.t[end],
        Tmin=minimum(tension[ok]),
        Tmax=maximum(tension[ok]),
        αmax=maximum(u -> abs(rem2pi(u[1], RoundNearest)), us),
        amax=maximum(norm, acc),
        vmax=maximum(norm, vel),
        κmax=maximum(κ),
        Rmin=1 / maximum(κ),
    )
end

# Keep a warm-started shooting guess inside the solver's admissible box.
clamp_guess(::Nothing) = (;)
function clamp_guess(s)
    return (;
        α0=clamp(s.α0, -π + 1e-6, π - 1e-6),
        dα0=clamp(s.dα0, -100 + 1e-6, 100 - 1e-6),
        dτ0=clamp(s.dτ0, 1e-6, 100 - 1e-6),
        tf=clamp(s.tf, 1.1e-3, 99.0),
        P=clamp(s.P, 1.1e-3, 9.9e4),
    )
end

# Floquet analysis of the Poincaré map: the sub-unit multiplier gives the
# exponential convergence rate toward the limit cycle, in digits per revolution.
const CB_POINC = LC.build_poincare_callback(eps(1.0))
function Φ_map(x, T; vbp)
    return Array(PM4.integrate(PM4.build_u(x[1], x[2], x[3], x[4], 0.0), T, vbp).u[end])[1:4]
end
function floquet_rate(sol, vbp)
    u0, T = LC.unpack_shooting(LC.build_shooting(sol))
    μ = eigvals(
        Matrix(ForwardDiff.jacobian(x -> Φ_map(x, T; vbp), collect(Array(u0))[1:4]))
    )
    λ = maximum(abs, μ[abs.(μ) .< 0.5])
    return (; slope=-log10(λ), λ, T)
end

# BifurcationKit shooting over the fundamental domain: the dynamics are cut off
# outside a velocity box, and a "twist" enforces one extra turn of tau.
const BK_M, BK_NSEEDS = 10, 100
vec4(u) = [u[1], u[2], u[3], u[4]]
function f_search!(du, u, p, t)
    if abs(u[3]) > BK_M + 2 || abs(u[4]) > BK_M + 2
        du .= 0
    else
        du .= PM4.dynamics(PM4.build_u(u[1], u[2], u[3], u[4], 0.0), p, t)[1:4]
    end
end
struct TwistedShooting{S,U} <: BK.AbstractBifurcationProblem
    sh::S
    u0::U
    twist::Float64
end
BK.getparams(w::TwistedShooting) = BK.getparams(w.sh)
BK.getlens(w::TwistedShooting) = BK.getlens(w.sh)
BK.getparam(w::TwistedShooting) = BK.getparam(w.sh)
BK.getu0(w::TwistedShooting) = w.u0
BK.getdelta(w::TwistedShooting) = BK.getdelta(w.sh)
function BK.residual(w::TwistedShooting, x, p)
    return (out=BK.po_residual(w.sh, x, p); out[end - 3] -= w.twist; out)
end
BK.jacobian(w::TwistedShooting, x, p) = BK.po_jacobian(w.sh, x, p)
root_state(x) = (α=rem2pi(x[1], RoundNearest), dα=x[3], dτ=x[4], T=x[5])
const npar = BK.NewtonPar(; tol=1e-9, max_iterations=15, verbose=false)
function same_cycle(a, b)
    return abs(rem2pi(a.α - b.α, RoundNearest)) < 1e-4 &&
           abs(a.dα - b.dα) < 1e-4 &&
           abs(a.dτ - b.dτ) < 1e-4 &&
           abs(a.T - b.T) < 1e-4
end

# Optimization study settings, shared by every section that touches Θ.
const TOL_OPT = 1e-5
const WIND_REF = 9.0
const SENSE_OPT = (+)
const SYMS_OPT = [:l, :I, :b, :Δφ]
const LOWER_OPT = [15.0, 62.0, 1.0, deg2rad(10.0)]
const UPPER_OPT = [200.0, 1.0e6, 1.0e6, deg2rad(90.0)]
nothing
