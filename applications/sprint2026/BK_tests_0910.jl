# Refactor of `BK_parametric_fixed.jl` (see that file for the physical-time BVP
# rationale and the multiple-shooting determinism bug). Changes:
#
# 1. ONE-LINE METHOD SWITCH: `build_disc()` returns the adaptive-mesh
#    collocation or the multiple-shooting discretizer; everything downstream
#    (model, jacobian type, discretize kwargs, pre-solve tolerance, continuation
#    box, orbit plotting) is derived from `typeof(disc)` by multiple dispatch.
# 2. The `FastPara` struct is replaced by a flat `NamedTuple` (`dynamics` reads
#    its parameters only by getproperty). Measured: 0.49 μs per `dynamics` call
#    (NamedTuple) vs 0.61 μs (ComponentArray) vs 0.46 μs (typed struct), zero
#    allocations; cache-miss construction ~5 μs, dominated by `build_vbpara`.
# 3. `T0` is derived from the OPTIMIZED parameter set (`lmt(vbp)[3]`) —
#    numerically identical (`l, v_ref ∉ SYMS`) but robust to `SYMS` changes.
#    Unit convention: `TF` (normalized) × `T0` = physical period at the
#    reference wind; the period recorded along branches (5th state component)
#    is physical seconds.
# 4. Both discretizations continue with `bothside=true` (correct in this BK
#    build). With BK ≤ 0.8.2 the leg merge reversed the backward branch,
#    masking sheet-jumps and mixing per-direction budgets — if that reappears,
#    use the explicit two-direction pattern of `BK_parametric_fixed.jl`.
# 5. Constants are CAPITAL_CASE to disambiguate them from mutable globals.
# 6. LOADABLE WITHOUT RUNNING. The expensive pipeline lives in `run(method)`
#    (plus the phase helpers below); `make_setup()` is called once and its
#    artifacts are threaded through explicitly instead of via mutable globals.
#    Running the file as a script (`julia --project=. BK_tests_0910.jl`) still
#    executes `run(build_disc())`; `include`-ing it only defines. The benchmark
#    harness `benchmark/bench.jl` relies on this. Run with an explicit
#    `--project`: the project is chosen by the launcher, not hardcoded here.

import OrdinaryDiffEqTsit5 as ODE
import OrdinaryDiffEqVerner as ODEV
using BifurcationKit
using LinearAlgebra
using ComponentArrays: ComponentArray as CA
using StaticArrays
using Plots

import KEEP
using KEEP: TAU0
using KEEP.PointMass4: dynamics
using KEEP.PointMassPara: build_vbpara, build_para, lmt
using KEEP.Optimization: optimize
using KEEP.LimitCycle: shoot as lc_shoot

const BVP = BifurcationKit.BVP

## ===================================================================== ##
## Setup (method-independent): optimization + sampled orbit              ##
## ===================================================================== ##
const SYMS = (:r, :I_eq, :torque_slope)

"Sample the optimized cycle at physical time `t`: (α, τ, dα, dτ) in rad, rad/s."
function x_optimization_phys(solution_sim, T0, t)
    s = solution_sim(t / T0, idxs=1:4)
    return SA[s[1], s[2], T0*s[3], T0*s[4]]
end

"""
Solve the optimization problem (fixed kite at the default wind) and build the
shared artifacts. Expensive (`optimize` ≈ 60 s via IPOPT, plus `lc_shoot`): call
once and pass the returned `setup` around. `nt_p0` is the DEFAULT-parameter
`NamedTuple` (the optimization only supplies the warm-start orbit, not the BVP
parameters).

`opt` lets a caller skip IPOPT by reusing a previous run's
`(; shooting, params_opt)` (see `scratch/_opt_result_tmp.jl`); default `nothing`
re-runs the optimization.
"""
function make_setup(; factor=5, opt=nothing)
    vbp0 = build_vbpara()
    p0 = build_para(vbp0)
    if opt === nothing
        lb = p0[SYMS] ./ factor
        ub = p0[SYMS] .* factor
        solution, stats, _ = optimize(p0, SYMS, lb, ub)
        params_opt = solution.params
        shooting = solution[1:4]
    else
        params_opt = opt.params_opt
        shooting = opt.shooting
    end

    vbp = build_vbpara(CA(p0; params_opt...))
    solution_sim = lc_shoot(shooting, vbp, save_everystep=true)
    tf = solution_sim.t[end]     # normalized period [units of T0]
    T0 = lmt(vbp)[3]             # characteristic time [s] of the OPTIMIZED set
    tf_physical = tf * T0        # physical period [s] at the reference wind

    nt_p0 = NamedTuple(p0)::PARA_NT
    u0 = SA[x_optimization_phys(solution_sim, T0, 0.0)..., tf_physical]
    return (
        vbp=vbp,
        p0=p0,
        nt_p0=nt_p0,
        u0=u0,
        state_size=length(u0),
        tf=tf,
        T0=T0,
        tf_physical=tf_physical,
        solution_sim=solution_sim,
    )
end

## ===================================================================== ##
## BVP model — physical time                                             ##
## ===================================================================== ##
# State u = (α, τ, dα, dτ, tf), integration variable s ∈ [0, 1] spans one period:
# du/ds = tf * [dα, dτ, ddα, ddτ, 0] with dα, dτ in rad/s and ddα, ddτ in rad/s².
# Reference implementation: rebuilds a ComponentArray on EVERY call (~5 μs).
# Kept only for validation of `F_fast`.
function F_reference(u, params, t=0)
    α, τ, dα, dτ, tf = u
    T = params.l / params.v_ref  # characteristic time of the CURRENT parameters
    u_dyn = SA[α, τ, T*dα, T*dτ, 0]
    _, _, ddα, ddτ, _ = dynamics(u_dyn, build_vbpara(params))
    out = tf .* SA[dα, dτ, ddα/T^2, ddτ/T^2, 0]
    # OrdinaryDiffEq requires typeof(du) === typeof(u); the multiple-shooting
    # discretizer feeds plain Vector slices of the unknown vector, so match
    # the container to the input.
    return u isa SVector ? out : Vector(out)
end

# Flat NamedTuple instead of the former FastPara struct: `dynamics` reads its
# parameters only by getproperty, so the ComponentArray machinery drops out
# (see header for the measured comparison). ponytail: single-entry cache keyed
# on v_ref — continuation evaluates many RHS at one parameter value, so we
# rebuild only when v_ref actually changes.
#
# Concrete param types: `build_para`/`build_vbpara` declare an `Any` return
# (`const Para = VBPara = Any` in KEEP.PointMassPara), so without these
# assertions `NamedTuple(...)` infers as the ABSTRACT `NamedTuple`; the RHS then
# dynamic-dispatches and `F_fast` allocates ~1 KiB per evaluation. The
# assertions recover inference without touching the package (see RESULTS.md).
const PARA_NT = typeof(NamedTuple(build_para()))
const NTPARA = typeof(NamedTuple(build_vbpara()))
const PARAM_CACHE = Ref{Union{Nothing,Tuple{Float64,NTPARA}}}(nothing)

@inline function get_fast_para(params)
    vref = params.v_ref
    c = PARAM_CACHE[]
    if c !== nothing && c[1] == vref
        return c[2]
    end
    p = NamedTuple(build_vbpara(params))::NTPARA
    PARAM_CACHE[] = (vref, p)
    return p
end

function F_fast(u, params, t=0)
    α, τ, dα, dτ, tf = u
    T = params.l / params.v_ref
    u_dyn = SA[α, τ, T*dα, T*dτ, 0]
    _, _, ddα, ddτ, _ = dynamics(u_dyn, get_fast_para(params))
    out = tf .* SA[dα, dτ, ddα/T^2, ddτ/T^2, 0]
    # OrdinaryDiffEq requires typeof(du) === typeof(u); the multiple-shooting
    # discretizer feeds plain Vector slices of the unknown vector, so match
    # the container to the input.
    return u isa SVector ? out : Vector(out)
end

function g(u0, uT, p)
    return SA[
        uT[1] - u0[1]       # α loop
        uT[3] - u0[3]       # dα loop
        uT[4] - u0[4]       # dτ loop
        u0[2] - TAU0        # τ init
        uT[2] - TAU0 - 2π   # τ end
    ]
end

# Record the physical period (5th state component) along branches. Without an
# explicit record, `plot(br)` falls back to `norm(x)`, which is incomparable
# between the two discretizations (different unknown-vector layouts) — the
# "multiple-shooting values are much smaller" artifact. `raw_x` unwraps the
# saved-solution wrapper that the adaptive-mesh collocation stores per point.
raw_x(x) = x isa BifurcationKit.BVPSavedSolutionAndState ? BifurcationKit.saved_solution(x) : x
record_period(x, p; kwargs...) = (tf=raw_x(x)[5],)

## ===================================================================== ##
## ONE-LINE METHOD SWITCH: collocation (adaptive mesh) or shooting       ##
## ===================================================================== ##
"ODE algorithm for the shooting flow and the Floquet/monodromy integrators.
Selectable via the `BK_ODE_ALG` env var (Tsit5 | Vern7 | Vern9) for the solver
comparison; defaults to Vern9 — 2–3.5× fewer f-evals than Tsit5 on this RHS
(benchmark/RESULTS_solvers.md)."
const ODE_ALG = let name = get(ENV, "BK_ODE_ALG", "Vern9")
    name == "Tsit5" ? ODE.Tsit5() :
    name == "Vern7" ? ODEV.Vern7() :
    name == "Vern9" ? ODEV.Vern9() :
    error("BK_ODE_ALG must be Tsit5|Vern7|Vern9, got $(repr(name))")
end

"Selected discretizer. Uncomment the alternative line to switch method."
build_disc() = BVP.Collocation(Ntst=30, m=5, meshadapt=true)  # ← collocation
# build_disc() = BVP.Shooting(10, ODE_ALG, true)               # ← multiple shooting

"RHS-function model for collocation; Shooting requires an ODEProblem model."
make_model(::BVP.Collocation, setup) = BVP.BVPModel(F_fast, g; n=setup.state_size)
function make_model(::BVP.Shooting, setup)
    odeprob = ODE.ODEProblem(F_fast, setup.u0, (0, 1), setup.nt_p0)
    return BVP.BVPModel(odeprob, g; n=setup.state_size)
end

"Shooting needs tight integration tolerances; collocation takes no kwargs.
Use the passed `disc`, not the selected method, so switching stays one-line."
make_bvp(model, disc::BVP.Collocation) = BVP.discretize(model, disc)
make_bvp(model, disc::BVP.Shooting) = BVP.discretize(model, disc; abstol=1e-12, reltol=1e-10)

"Analytical Jacobian for collocation; shooting differentiates the pinned residual.
(`FullSparse` is not usable: BK's BVP collocation sparse path hits an
uninitialised block cache, `getindex(::Nothing, :, i)` — see RESULTS.md #3.)"
make_jac(::BVP.Collocation) = BifurcationKit.DenseAnalytical()
# Shooting: use the manual block Jacobian (see below). `ManualJacFwd` reproduces
# the central-difference `ManualJacFD` blocks to ~1e-7 and is ~3× faster at
# steady state; keep FD as the AD-free fallback.
make_jac(::BVP.Shooting) = ManualJacFwd()

"Pre-solve tolerance / iteration budget tuned per method."
make_presolve(::BVP.Collocation) = (tol=1e-9, max_iter=200)
make_presolve(::BVP.Shooting) = (tol=1e-6, max_iter=100)

"Continuation parameter box, tuned per method. `detect_bifurcation=0` skips the
per-step eigenvalue solve: Finding 1 shows BK's BVP spectra are unusable here and
the script never reads `br.eig` (ground truth is `floquet_true`), yet
`compute_eigenvalues` was ~half the continuation time (RESULTS.md #2)."
function make_contpar(::BVP.Collocation, optn)
    ContinuationPar(
        p_min=0.0,
        p_max=25.0,
        dsmax=0.1,
        ds=0.01,
        detect_bifurcation=0,
        newton_options=optn,
        max_steps=100,
    )
end
function make_contpar(::BVP.Shooting, optn)
    ContinuationPar(
        p_min=0.1,
        p_max=50.05,
        dsmax=0.1,
        ds=0.01,
        detect_bifurcation=0,
        newton_options=optn,
        max_steps=100,
    )
end

"""
Dense orbit `(t, u)` of the unknown vector `x`, per discretizer. Both go
through `BVP.get_solution_bvp`; Shooting's internal reshape chokes on the
trailing auxiliary unknown, so pass it the pure DOF block.
"""
function orbit(disc::BVP.Shooting, bvp, x, p)
    n = BVP.state_dimension(bvp)
    s = BVP.get_solution_bvp(bvp, @view(x[1:(n*disc.M)]), p)
    return s.t, s.u
end
function orbit(disc::BVP.Collocation, bvp, x, p)
    s = BVP.get_solution_bvp(bvp, x, p)
    return s.t, s.u
end

"Plot the cycle (dτ vs physical time) — single implementation for both methods."
function plot_solution(x, p; iter, state, k...)
    bvp = BVP.get_bvp(BifurcationKit.getprob(iter))
    t, u = orbit(BVP.get_discretizer(bvp), bvp, raw_x(x), BifurcationKit.getparams(iter, state))
    tf = raw_x(x)[5]
    plot!(t ./ t[end] .* tf, u[4, :]; ylabel="dτ", k...)
end

# Upstream bug fix (determinism), see `BK_parametric_fixed.jl` header:
# `out = similar(X)` in the Shooting residual leaves the trailing auxiliary
# unknown UNINITIALIZED, so Newton consumed random memory — hence
# "un objet est mal initialisé lors du solve, exécuter plusieurs fois donne
# différents résultats". We pin it to zero (residual row X[end] = 0), which
# also makes the dense Jacobian nonsingular (J[end,end] = 1). Only fires for
# the Shooting discretizer.
function BVP.bvp_residual(d_bvp::BVP.DiscretizedBVP{<:BVP.BVPModel,<:BVP.Shooting}, X, p)
    model_ = BVP.get_model(d_bvp)
    disc_ = BVP.get_discretizer(d_bvp)
    n = BVP.state_dimension(model_)
    t0, tf = BVP.get_time_interval(model_)
    M = BVP.mesh_size(disc_)

    Xm = reshape(@view(X[1:(n*M)]), n, M)
    T = tf - t0

    out = similar(X)
    outm = reshape(@view(out[1:(n*M)]), n, M)
    BVP.bvp_residual_bare!(d_bvp, outm, Xm, p, T)
    out[end] = X[end]  # pin the auxiliary unknown to zero (deterministic residual)
    return out
end

## --------------------------------------------------------------------- ##
## Shooting Jacobian — manual block assembly                              ##
## --------------------------------------------------------------------- ##
# BK's `AutoDiffDense` differentiates the whole 51-vector residual with
# ForwardDiff, i.e. THROUGH the adaptive ODE solve — and `dynamics` already
# calls ForwardDiff internally, so the duals nest and one Jacobian takes
# essentially forever. The shooting Jacobian is block-structured, so assemble
# it directly. For unknowns X = [u₁ … u_M ; aux] (n·M + 1):
#   matching row block i = 1…M-1 :  Φᵢ (diagonal) ,  -I (to uᵢ₊₁)
#   boundary row block  M        :  ∂g/∂u₀  ,  (∂g/∂u_T) Φ_M
#   pinned row                   :  [0 … 0 1]
# with Φᵢ = ∂φᵢ/∂uᵢ the segment transition matrix. `ManualJacFwd` gets the
# blocks from ForwardDiff (still AD through the solve — tried first, measured),
# `ManualJacFD` from central finite differences (no AD; the method the
# period-map helpers already use).
struct ManualJacFwd end
struct ManualJacFD end

"Central finite-difference Jacobian of `f` at `x` (buffer reused across columns)."
function fd_mat(f, x; h=1e-7)
    n = length(x)
    J = Matrix{Float64}(undef, n, n)
    e = zeros(n)
    for i in 1:n
        e[i] = h
        J[:, i] = (f(x .+ e) .- f(x .- e)) ./ (2h)
        e[i] = 0.0
    end
    return J
end

block_jac(::ManualJacFwd, f, x) = BifurcationKit.ForwardDiff.jacobian(f, x)
block_jac(::ManualJacFD, f, x) = fd_mat(f, x)

function BVP.bvp_jacobian(d_bvp::BVP.DiscretizedBVP{<:BVP.BVPModel,<:BVP.Shooting},
    jac::Union{ManualJacFwd,ManualJacFD}, X, p)
    model_ = BVP.get_model(d_bvp)
    disc_ = BVP.get_discretizer(d_bvp)
    sh = BVP.get_cache(d_bvp)
    n = BVP.state_dimension(model_)
    M = BVP.mesh_size(disc_)
    t0, tf = BVP.get_time_interval(model_)
    T = tf - t0
    U = reshape(@view(X[1:(n*M)]), n, M)
    N = n * M + 1
    J = zeros(eltype(X), N, N)
    In = Matrix{eltype(X)}(LinearAlgebra.I, n, n)
    # Segment end-state, matching whichever flow the residual uses: the parallel
    # `EnsembleProblem` path integrates the whole matrix at once, the serial path
    # one vector at a time.
    parallel = BVP.is_parallel(disc_)
    flow = if parallel
        (u, i) -> BifurcationKit.evolve(sh.flow, reshape(u, n, 1), p, [sh.ds[i] * T])[1].u
    else
        (u, i) -> BifurcationKit.evolve(sh.flow, u, p, sh.ds[i] * T).u
    end
    for i in 1:(M-1)
        ri = (i - 1) * n .+ (1:n)
        J[ri, ri] .= block_jac(jac, u -> flow(u, i), U[:, i])
        J[ri, i*n .+ (1:n)] .= -In
    end
    u1, uM = U[:, 1], U[:, M]
    uT = flow(uM, M)
    rM = (M - 1) * n .+ (1:n)
    J[rM, 1:n] .= block_jac(jac, u -> model_.g(u, uT, p), u1)
    J[rM, rM] .+= block_jac(jac, u -> model_.g(u1, u, p), uT) *
                  block_jac(jac, u -> flow(u, M), uM)
    J[end, end] = one(eltype(X))
    return J
end

"""
Armijo-backtracking Newton for a BVPBifProblem. BifurcationKit's plain `Newton`
only takes full steps, which overshoot on this stiff problem (from the raw guess,
one step moves ‖δ‖∞ ≈ 47 and diverges); backtracking enters the quadratic basin.
Trial points whose integration fails count as infinite residuals, and the step
is capped at ‖Δ‖∞ = Δmax to keep them dynamically feasible.
"""
function damped_newton(prob, x0, params; tol=1e-8, max_iter=200, Δmax=1.0)
    x = copy(x0)
    res_trial(x_) =
        try
            norminf(BifurcationKit.residual(prob, x_, params))
        catch
            Inf
        end
    for _ in 1:max_iter
        r = BifurcationKit.residual(prob, x, params)
        nr = norminf(r)
        nr < tol && return x, nr
        δ = try
            BifurcationKit.jacobian(prob, x, params) \ (-r)
        catch
            break
        end
        all(isfinite, δ) || break
        λ = min(1.0, Δmax / (norminf(δ) + eps()))
        while λ > 1e-6 && res_trial(x .+ λ .* δ) > nr
            λ /= 2
        end
        x .+= λ .* δ
    end
    return x, norminf(BifurcationKit.residual(prob, x, params))
end

"Recorded PALC points satisfy |Δp| ≤ dsmax; larger gaps flag corrector
sheet-jumps (only possible where sheets crowd, i.e. near folds). Uses the
branch param column — `br.sol` is a single saved solution, not the path."
function sheet_jumps(br; factor=1.5)
    dsmax = br.contparams.dsmax
    return [i for i in 2:length(br)
                  if abs(br.branch.param[i] - br.branch.param[i-1]) > factor * dsmax]
end

## Build the problem
# Initial guess. Collocation is robust to the optimized-vs-default parameter
# mismatch in the sampled optimization orbit and converges from it; for shooting
# that raw sample leaves the damped pre-solve stuck (residual O(1)). So, for
# shooting only, first obtain a converged cycle with collocation at the SAME
# (default) parameters and sample THAT into the shooting points — the warm start
# the original `BK_parametric_fixed.jl` used.
function collocation_cycle(setup)
    disc_c = BVP.Collocation(Ntst=30, m=5, meshadapt=true)
    bvp_c = make_bvp(make_model(disc_c, setup), disc_c)
    x0_c = BVP.generate_solution(bvp_c,
        t -> vcat(x_optimization_phys(setup.solution_sim, setup.T0, setup.tf_physical * t), setup.tf_physical))
    prob_c = BVP.BVPBifProblem(bvp_c, x0_c, setup.nt_p0, (@optic _.v_ref);
        jacobian=make_jac(disc_c), record_from_solution=record_period)
    x0_c, res_c = damped_newton(prob_c, x0_c, setup.nt_p0; make_presolve(disc_c)...)
    prob_c = BifurcationKit.re_make(prob_c; u0=x0_c)
    sol_c = BifurcationKit.solve(prob_c, Newton(), NewtonPar(tol=1e-10, linesearch=true))
    @assert BifurcationKit.converged(sol_c) "collocation warm start did not converge"
    println("collocation warm-start residual: ", res_c)
    return bvp_c, sol_c
end

"Sample a converged collocation orbit into the M shooting points (+ dummy)."
function shooting_warm_start(bvp_s, disc_s, setup)
    bvp_c, sol_c = collocation_cycle(setup)
    bs = BVP.get_solution_bvp(bvp_c, sol_c.u, setup.nt_p0)
    ts = collect(bs.t);
    ts ./= ts[end]
    n = BVP.state_dimension(bvp_s)
    M = BVP.mesh_size(disc_s)
    x = zeros(n * M + 1)
    for i in 1:M
        s = (i - 1) / M
        j = clamp(searchsortedfirst(ts, s) - 1, 1, length(ts) - 1)
        θ = (s - ts[j]) / (ts[j+1] - ts[j])
        x[(i-1)*n .+ (1:n)] .= (1 - θ) .* view(bs.u, :, j) .+ θ .* view(bs.u, :, j + 1)
    end
    x[end] = 0.0  # auxiliary unknown, pinned by the residual override
    return x
end

## --------------------------------------------------------------------- ##
## Floquet multipliers — ground truth (period-map finite differences)      ##
## --------------------------------------------------------------------- ##
# DO NOT read stability from `br.eig`. The collocation jacobian is the full
# discretized problem; empirically ~660 of its ~755 eigenvalues are SPURIOUS
# collocation modes (|mu| ≈ 11, Re ≈ −11), and BK keeps only the top-`nev`
# by REAL PART. The recorded set is a mix of spurious and physical modes and
# the `1+λ` / `Re` counting in the bvp module does not yield Floquet
# multipliers. (A first pass reading `br.eig` produced a bogus
# "Neimark–Sacker at v_ref ≈ 7" story; the ground truth below contradicts it.)
#
# The genuine multipliers come from the period map: integrate the cycle over
# s ∈ [0,1] and differentiate the endpoint w.r.t. the initial (α, τ, dα, dτ)
# by central differences. The result then contains the exact trivial
# multiplier 1; all others are the real Floquet multipliers.

"4 mechanical dofs in physical time (velocities are physical rad/s)."
function phys_rhs4(v, p)
    T = p.l / p.v_ref
    α, τ, dα, dτ = v
    _, _, ddα, ddτ, _ = dynamics(SA[α, τ, T*dα, T*dτ, 0.0], get_fast_para(p))
    return SA[dα, dτ, ddα/T^2, ddτ/T^2]
end

# `fd_mat` (above) is the single implementation; `fd_jac` is its name in the
# Floquet/monodromy helpers.
const fd_jac = fd_mat

"""
Floquet multipliers of the converged cycle `x` at parameters `p` (the BVP
period `tf` is held fixed). Returns four multipliers; the one closest to +1
is the trivial multiplier.
"""
function floquet_true(disc, bvp, x, p)
    u = orbit(disc, bvp, raw_x(x), p)[2]
    tf = raw_x(x)[5]
    v0 = u[1:4, 1]
    pm = v -> begin
        pr = ODE.ODEProblem((uu, pp, t) -> F_fast(uu, pp, t), vcat(v, tf), (0.0, 1.0), p)
        s = ODE.solve(pr, ODE_ALG, abstol=1e-13, reltol=1e-13)
        copy(Array(s.u[end])[1:4])
    end
    return LinearAlgebra.eigvals(fd_jac(pm, v0; h=1e-7))
end

"Max |mu| over the non-trivial multipliers (trivial = closest to +1)."
function floquet_max_true(disc, bvp, x, p)
    mu = floquet_true(disc, bvp, x, p)
    j = argmin(abs.(mu .- 1))
    return maximum(abs(mu[k]) for k in eachindex(mu) if k != j)
end

"""
Flow of the 4 mechanical state components over one segment `ds` of the
NORMALISED BVP variable (`F_fast` advances that variable; physical dt = tf·ds),
holding the period `tf` fixed. Building block of `monodromy`.
"""
function flowseg(v, tf, ds, p)
    pr = ODE.ODEProblem((uu, pp, t) -> F_fast(uu, pp, t), vcat(v, tf), (0.0, ds), p)
    return copy(Array(ODE.solve(pr, ODE_ALG, abstol=1e-13, reltol=1e-13).u[end])[1:4])
end

"""
Segment-wise monodromy of a saved orbit: the ordered product of the per-segment
flow Jacobians, `M = J_k ⋯ J_1`. Unlike the single period-map finite difference
(`floquet_true`), this stays well conditioned for the strongly unstable long
saddle cycle, over which one-period integration misses closure by O(1) and the
naive FD returns garbage (see Finding 9). `u` is the 5×N state matrix from
`orbit`, `t` its NORMALISED BVP time grid (0 → 1; a segment of normalised length
`ds` advances the state by physical `tf·ds`), and `tf` the fixed period.

NOTE: the per-saved-step partition is REQUIRED. Coarsening it (fewer, longer FD
blocks, which is algebraically equivalent by the chain rule) keeps the short
cycle's |mu| to ~1e-6 but corrupts the long cycle's (1.6e7 → 2.6e7…1.7e8):
central differences over a long segment of the strongly unstable saddle lose
accuracy. Tried and rejected — benchmark/RESULTS.md #4.
"""
function monodromy(u, t, tf, p)
    M = Matrix{Float64}(I, 4, 4)
    for k in 1:(length(t)-1)
        ds = t[k+1] - t[k]          # already normalised (orbit output)
        ds <= 0 && continue
        M = fd_jac(x -> flowseg(x, tf, ds, p), u[1:4, k]; h=1e-7) * M
        all(isfinite, M) || break
    end
    return M
end

"Max |mu| over the non-trivial multipliers of a saved orbit, via `monodromy`."
function monodromy_max(u, t, tf, p)
    mu = LinearAlgebra.eigvals(monodromy(u, t, tf, p))
    j = argmin(abs.(mu .- 1))
    return maximum(abs(mu[k]) for k in eachindex(mu) if k != j)
end

"""
Fit the low-wind terminus of `(p, tf)` as a fold, `p = p_f + γ (tf - tf_f)²`,
returning the best `(rms, tff, pf, γ)`. A fold is the square-root turning point
of a saddle-node; a homoclinic instead gives `tf = c - b ln(p - p_f)`
(`fit_homoclinic`). Comparing the two residuals decides which (Finding 3).
"""
function fit_fold(pp, tt)
    best = (rms=Inf, tff=NaN, pf=NaN, γ=NaN)
    for tff in range(minimum(tt) - 2, maximum(tt) + 2; length=2001)
        x = (tt .- tff) .^ 2
        X = hcat(ones(length(x)), x)
        coef = X \ pp
        rms = sqrt(sum(abs2, pp .- X * coef) / length(pp))
        rms < best.rms && (best = (rms=rms, tff=tff, pf=coef[1], γ=coef[2]))
    end
    return best
end

"""
Fit the low-wind terminus as a homoclinic, `p = p_f + A exp(-tf/b)`, returning
the best `(rms, pf, A, b)`. The `log` needs `p > p_f`, so `p_f` is searched
just below the smallest `p` in the window.
"""
function fit_homoclinic(pp, tt)
    best = (rms=Inf, pf=NaN, A=NaN, b=NaN)
    for pf in range(minimum(pp) - 0.05, minimum(pp) - 1e-6; length=2001)
        y = log.(pp .- pf)
        X = hcat(ones(length(y)), tt)
        coef = X \ y
        rms = sqrt(sum(abs2, y .- X * coef) / length(y))
        rms < best.rms && (best = (rms=rms, pf=pf, A=exp(coef[1]), b=-1 / coef[2]))
    end
    return best
end

## --------------------------------------------------------------------- ##
## Phase helpers                                                          ##
## --------------------------------------------------------------------- ##
"Low-wind terminus: fold (saddle-node of cycles) vs homoclinic. Only the
shooting branch reaches this window — collocation stalls before the fold."
function lower_terminus(disc, br; verbose=true)
    if !(disc isa BVP.Shooting)
        verbose && println("lower-terminus fit: skipped for ", nameof(typeof(disc)), " (does not pass the fold)")
        return nothing
    end
    # The family folds at low wind: p reaches a minimum p_f at finite tf and
    # turns back. Fit both candidate laws to the low-period window and compare
    # residuals; a fold (√ turning point) wins by orders of magnitude (Finding 3).
    pp, tt = collect(br.branch.param), collect(br.branch.tf)
    m = (tt .> 4.0) .& (tt .< 30.0) .& (pp .> 2.4580) .& (pp .< 2.55)
    if sum(m) <= 5
        verbose && println("lower-terminus fit: shooting branch has too few low-wind points (", sum(m), ")")
        return nothing
    end
    A = fit_fold(pp[m], tt[m])
    B = fit_homoclinic(pp[m], tt[m])
    i0 = argmin(pp)   # true discrete turning point (may sit outside the fit window)
    if verbose
        println("lower-terminus fit over ", sum(m), " points (tf ∈ [",
            round(minimum(tt[m]), digits=2), ", ", round(maximum(tt[m]), digits=2),
            "] s, p ∈ [", round(minimum(pp[m]), digits=4), ", ",
            round(maximum(pp[m]), digits=4), "]):")
        println("  fold       : p = p_f + γ(tf - tf_f)²   RMS(Δp)   = ",
            round(A.rms, sigdigits=3), "   vertex p_f = ", round(A.pf, digits=5),
            " (tf_f = ", round(A.tff, digits=3), " s, extrapolated)")
        println("  homoclinic : p = p_f + A exp(-tf/b)    RMS(Δlog) = ",
            round(B.rms, sigdigits=3), "   p_f = ", round(B.pf, digits=5))
        println("  discrete minimum: p = ", round(pp[i0], digits=5), " at tf = ",
            round(tt[i0], digits=2), " s")
    end
    return (fold=A, homoclinic=B)
end

"""
Continue the cycle in v_ref by damped Newton (warm start, one solution per
step) and return rows `(v_ref, tf, res, |mu|)` with `|mu|` sorted descending.
Runs until the corrector stops converging — the branch's numerical horizon.
"""
function stability_sweep(disc, bvp, x0, vrefs, setup; tol=1e-9)
    x = x0
    rows = NamedTuple[]
    for vr in vrefs
        p = merge(setup.nt_p0, (v_ref=vr,))
        prb = BVP.BVPBifProblem(bvp, x, p, (@optic _.v_ref);
            jacobian=make_jac(disc), record_from_solution=record_period, plot_solution=plot_solution)
        x, res = damped_newton(prb, x, p; tol=tol, max_iter=300)
        push!(rows, (v_ref=vr, tf=raw_x(x)[5], res=res, mu=sort(abs.(floquet_true(disc, bvp, x, p)), rev=true)))
    end
    return rows
end

"Landscape: branch geometry (top) and TRUE Floquet stability (bottom)."
function landscape(disc, bvp, x0, br, setup, plt_land; make_plots=true)
    pc, tfc = br.branch.param, br.branch.tf
    # True multipliers, warm-started Newton on both sides of v_ref = 9. Rows with
    # res ≳ 1e-6 did not converge (past the branch's numerical horizon) → dropped.
    rows_down = stability_sweep(disc, bvp, x0, collect(9.0:-0.25:2.5), setup)
    rows_up = stability_sweep(disc, bvp, x0, collect(9.5:0.5:17.0), setup)
    good(rows) = [r.res < 1e-6 ? r.mu[2] : NaN for r in rows]   # mu sorted desc, [1] = trivial

    if make_plots
        pa = plot(pc, tfc; label=string(nameof(typeof(disc))), xlabel="v_ref", ylabel="tf [s]",
            title="limit-cycle branch (default params, optimization warm start)")
        pb = plot([r.v_ref for r in rows_down], good(rows_down); label="from v_ref=9 ↓",
            color=1, yscale=:log10, xlabel="v_ref", ylabel="max nontrivial |mu|",
            title="Floquet stability (ground truth)", ylims=(1e-4, 10))
        plot!(pb, [r.v_ref for r in rows_up], good(rows_up); label="from v_ref=9 ↑", color=2)
        hline!(pb, [1.0]; label="|mu| = 1 (stability boundary)", color=:gray, ls=:dash)
        plt_land = plot(pa, pb; layout=(2, 1), size=(650, 750))
        savefig(plt_land, joinpath(@__DIR__, "BK_tests_0910_landscape.png"))
    end
    return rows_down, rows_up, plt_land
end

"Coexisting limit cycles (independent of continuation): Poincaré sampling.
Reveals arcs a single continuation cannot see. Reuses KEEP.LimitCycle."
function coexist_cycles(setup; verbose=true)
    rows = []
    for vr in (9.0, 5.0, 3.0)
        vbp = build_vbpara(merge(setup.nt_p0, (v_ref=vr,)))
        for lc in KEEP.LimitCycle.all_limit_cycles(vbp; αmin=(-π), αmax=π, vmax=12, N=40)
            s = KEEP.LimitCycle.build_shooting(lc)
            push!(rows, (v_ref=vr, α0=s[1], dα0=s[2], dτ0=s[3], T=s[4],
                power=lc.u[end][5] / lc.t[end]))
        end
    end
    if verbose
        println("coexisting cycles (v_ref, α0, dα0, dτ0, T, power):")
        foreach(println, rows)
    end
    return rows
end

"""
Two cycles at ONE v_ref: the stable "short" cycle and the "long" saddle cycle
past the fold. `br.sol` saves every continuation step (save_sol_every_step
defaults to 1), so BOTH crossings of v_ref — descending on the short branch,
ascending on the long branch — are already present, with no separate
arclength driver. Requires having passed the fold: shooting does, collocation
(this build) stalls just above it.
"""
function two_cycles(disc, bvp, br, setup; VRT=2.47, make_plots=true, verbose=true)
    discname = nameof(typeof(disc))
    # dsmax=0.1 steps may skip VRT, so widen the window, then refine candidate
    # cycles at exactly VRT with a fixed-parameter damped Newton. The two-cycle
    # decision is made AFTER refinement: the raw tf spread over the window is
    # merely the branch's tf(p) variation (points from both legs share one sheet
    # when the fold was not passed), so only genuinely separated refined periods
    # count as two cycles.
    cand = [(i, s.p, raw_x(s.x)[5]) for (i, s) in enumerate(br.sol) if abs(s.p - VRT) < 0.12]
    cand = [c for c in cand if isfinite(last(c))]
    pv = merge(setup.nt_p0, (v_ref=VRT,))
    refine(x) = (pr=BVP.BVPBifProblem(bvp, x, pv, (@optic _.v_ref);
            jacobian=make_jac(disc), record_from_solution=record_period);
        damped_newton(pr, x, pv; tol=1e-9)[1])
    if isempty(cand)
        verbose && println("two-cycles: no candidate near v_ref=", VRT, " (", discname, ")")
        return nothing
    end
    sort!(cand, by=last)
    xshort = refine(raw_x(br.sol[cand[1][1]].x))
    xlong = refine(raw_x(br.sol[cand[end][1]].x))
    tf1, tf2 = raw_x(xshort)[5], raw_x(xlong)[5]
    if max(tf1, tf2) / min(tf1, tf2) < 1.3
        verbose && println("two-cycles: only ONE branch near v_ref=", VRT, " (", discname,
            " did not pass the fold): tf ≈ ",
            string(round(minimum((tf1, tf2)), digits=2)), " s")
        return nothing
    end
    t1, u1 = orbit(disc, bvp, raw_x(xshort), pv)
    t2, u2 = orbit(disc, bvp, raw_x(xlong), pv)
    # Segment-wise monodromy for BOTH crossings so the numbers are comparable.
    # The short cycle could also be measured by the naive period-map FD, but the
    # long saddle cycle cannot (it is unstable enough that one-period integration
    # misses closure by O(1) — Finding 9).
    mu1 = monodromy_max(u1, t1, tf1, pv)
    mu2 = monodromy_max(u2, t2, tf2, pv)
    if verbose
        println("two-cycles at v_ref=", VRT, ": max non-trivial |mu| = ",
            round(mu1, sigdigits=3), " (short), ", round(mu2, sigdigits=3), " (long)")
    end
    if make_plots
        plt2 = plot(layout=(2, 1), size=(800, 620))
        plot!(plt2[1], t1 ./ t1[end] .* tf1, u1[1, :]; label="short tf=$(round(tf1, digits=2)) s")
        plot!(plt2[1], t2 ./ t2[end] .* tf2, u2[1, :]; label="long  tf=$(round(tf2, digits=2)) s", lw=2,
            xlabel="t [s]", ylabel="α [rad]", title="$discname: two cycles at v_ref=$VRT")
        plot!(plt2[2], u1[1, :], u1[4, :]; label="short")
        plot!(plt2[2], u2[1, :], u2[4, :]; label="long", lw=2,
            xlabel="α [rad]", ylabel="dτ [rad/s]", title="phase portrait")
        savefig(plt2, joinpath(@__DIR__, "BK_tests_0910_two_cycles.png"))
    end
    verbose && println("two-cycles at v_ref=", VRT, ": tf = ", round(tf1, digits=2), " s (short), ",
        round(tf2, digits=2), " s (long); saved BK_tests_0910_two_cycles.png")
    return (tf_short=tf1, tf_long=tf2, mu_short=mu1, mu_long=mu2)
end

## ===================================================================== ##
## Full pipeline                                                          ##
## ===================================================================== ##
"Time `f()`; print wall/alloc/GC when `verbose`. Returns `(value, @timed result)`."
function timed(label, f; verbose=true)
    r = @timed f()
    verbose && println(label, ": ", round(r.time, digits=2), " s, ",
        round(r.bytes / 2^20, digits=1), " MiB, gc ", round(100r.gctime / r.time, digits=1), "%")
    return r.value, r
end

"""
Run the full study for discretizer `method` (a `BVP.Collocation` or
`BVP.Shooting` instance). Returns the artifacts needed for inspection and
benchmarking. `setup` defaults to a fresh `make_setup()`.
"""
function run(method=build_disc(); setup=make_setup(), make_plots=true, verbose=true)
    discname = nameof(typeof(method))
    model = make_model(method, setup)
    bvp = make_bvp(model, method)
    x0 = method isa BVP.Shooting ? shooting_warm_start(bvp, method, setup) :
         BVP.generate_solution(bvp,
        t -> vcat(x_optimization_phys(setup.solution_sim, setup.T0, setup.tf_physical * t), setup.tf_physical))

    prob = BVP.BVPBifProblem(bvp, x0, setup.nt_p0, (@optic _.v_ref);
        jacobian=make_jac(method),
        record_from_solution=record_period,
        plot_solution=plot_solution,
    )

    optn = NewtonPar(tol=1e-10, verbose=true, linesearch=true)

    pre = @timed damped_newton(prob, x0, setup.nt_p0; make_presolve(method)...)
    x0, res_pre = pre.value
    t_pre = pre.time
    verbose && println(discname, " damped pre-solve residual: ", res_pre, "  (", round(t_pre, digits=2), " s)")
    prob = BifurcationKit.re_make(prob; u0=x0)

    sol, tim_solve = timed("solve", () -> BifurcationKit.solve(prob, Newton(), optn); verbose)
    @assert BifurcationKit.converged(sol) "$discname Newton did not converge"

    ## Sanity: converged cycle at the reference wind (exercises the orbit plumbing)
    t_orb, u_orb = orbit(method, bvp, raw_x(sol.u), setup.nt_p0)
    if make_plots
        plot(t_orb ./ t_orb[end] .* raw_x(sol.u)[5], u_orb[4, :];
            label=string(discname), xlabel="t [s]", ylabel="dτ [rad/s]",
            title="converged cycle at v_ref = $(setup.nt_p0.v_ref)")
        savefig(joinpath(@__DIR__, "BK_tests_0910_orbit.png"))
    end

    ## Continuation
    # bothside=true: correct in this BK build — with BK ≤ 0.8.2 the merged branch
    # reversed the backward leg and masked sheet-jumps; if that reappears, use the
    # explicit two-direction pattern of `BK_parametric_fixed.jl` (header, item 4).
    optc = make_contpar(method, optn)
    br, tim_cont = timed("continuation",
        () -> continuation(prob, PALC(), optc; plot=false, verbosity=1, normC=norminf, bothside=true); verbose)

    verbose && println(discname, " sheet-jumps at indices ", sheet_jumps(br))
    if make_plots
        plot(br; label=string(discname))
        savefig(joinpath(@__DIR__, "BK_tests_0910_branch.png"))
    end

    lower_terminus(method, br; verbose)

    land = @timed landscape(method, bvp, sol.u, br, setup, nothing; make_plots)
    rows_down, rows_up, _ = land.value
    t_land = land.time
    verbose && println("stability sweep (ground truth): ", round(t_land, digits=2), " s")

    co = @timed coexist_cycles(setup; verbose)
    coexist = co.value
    t_coexist = co.time
    verbose && println("coexist: ", round(t_coexist, digits=2), " s")

    tw = @timed two_cycles(method, bvp, br, setup; make_plots, verbose)
    two = tw.value
    t_two = tw.time
    verbose && println("two-cycles: ", round(t_two, digits=2), " s")

    timings = (
        presolve=t_pre, solve=tim_solve, continuation=tim_cont,
        stability=t_land, coexist=t_coexist, twocycles=t_two,
    )
    return (; setup, method, discname, bvp, prob, sol, br, rows_down, rows_up, coexist, two, timings)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run(build_disc())
end

## --------------------------------------------------------------------- ##
## Findings                                                              ##
## --------------------------------------------------------------------- ##
# 1. BK'S BVP SPECTRA ARE UNUSABLE FOR STABILITY HERE. The collocation
#    jacobian is the full ~755×755 discretized operator; ~660 of its
#    eigenvalues are spurious collocation modes (|mu| ≈ 11, Re ≈ −11). BK
#    keeps only the top-`nev` by REAL PART, so `br.eig` mixes spurious and
#    physical modes; the shooting spectrum is similarly polluted (our
#    residual pin injects an exact eigenvalue 1). Neither gives Floquet
#    multipliers — do not read stability/Neimark–Sacker from them.
#
# 2. GROUND-TRUTH STABILITY (period-map finite differences, `floquet_true`).
#    The cycle has a trivial mu = 1 and three nontrivial multipliers. Over
#    the whole reachable branch (v_ref ≈ 2.5 → 17) every nontrivial |mu| is
#    far below 1 → the cycle is STRONGLY STABLE everywhere. The largest
#    nontrivial |mu| grows from ≈0.002 near v_ref ≈ 8 to 0.25 at
#    v_ref ≈ 2.48 and 0.011 at v_ref = 17, but never approaches 1: no fold,
#    no Neimark–Sacker, no period doubling inside the branch.
#
# 3. THE LOWER TERMINUS IS A FOLD OF LIMIT CYCLES (saddle-node), NOT A
#    HOMOCLINIC. The shooting branch `brS` (p, tf profile) has exactly one
#    turning point: p decreases to a minimum p = 2.4585 at a FINITE period
#    tf ≈ 14.6 s, then increases again while tf keeps growing to ≥33.7 s. A
#    limit point in v_ref forces the shooting Jacobian singular, i.e. a
#    nontrivial Floquet multiplier → +1 (fold). A homoclinic would instead
#    have tf → ∞ with v_ref decreasing monotonically and no turnaround.
#    `fit_fold`/`fit_homoclinic` (called after the branch for the Shooting
#    discretization) confirm this quantitatively: the √ turning-point law fits
#    the low-wind window far better than the homoclinic log-law.
#    The long-period side is a SECOND arc (the saddle branch) reached only
#    after the fold; collocation stalls before it (v_ref ≈ 2.47, tf ≈ 11 s).
#
# 4. COEXISTING CYCLES ARE A MIRROR PAIR. `KEEP.LimitCycle.all_limit_cycles`
#    (Poincaré sampling + clustering) finds TWO cycles per v_ref, related by
#    the system symmetry (α0, dτ0) → (−α0, −dτ0): e.g. at v_ref = 9,
#    (α0, P) = (−0.777, 9055 W) and (0.674, 8267 W); at v_ref = 3,
#    (∓0.30, ≈150 W). A single continuation sees only one of them; enumerate
#    with the sampler.
#
# 5. EQUILIBRIA AND HOPF. `KEEP.SteadyState.all_steady_states_halton` gives
#    8 equilibria (4 mirror pairs) at every v_ref ∈ [0.5, 15] tested. ALL are
#    unstable everywhere: the operating cycle (α0 ≈ −0.78) wraps an unstable
#    focus at α ≈ −0.87 whose complex pair decays as the wind drops
#    (Re ≈ +1.7 at v_ref = 10 → +0.14 at v_ref = 1.25, grazing 0 near
#    v_ref ≈ 1) but a real unstable mode remains, so there is NO clean Hopf
#    generating a stable cycle in the reachable range. Equilibria explain the
#    dynamics (the kite is always on a cycle) but do not localise the arc's
#    lower terminus.
#
# 6. REUSE KEEP's TOOLS. src/SteadyState.jl (`steady_state`,
#    `all_steady_states_halton`, `ddq_partial`) and src/LimitCycle.jl
#    (`all_limit_cycles`, `shoot`, `build_shooting`, `unpack_shooting`) do
#    equilibrium and cycle discovery properly — prefer them to ad-hoc Newton.
#
# 7. DIFFERENTIATION BACKEND / MEMORY. Use FORWARD-mode AD only at the single
#    level KEEP already does (PointMass4.jacobian = ForwardDiff). Do NOT wrap
#    the period map in ForwardDiff: `dynamics` already calls ForwardDiff
#    internally, so AD-through-the-solver nests duals inside the adaptive
#    integrator and blows up memory (observed >1 GB, crashed the machine).
#    Central finite differences on the 4×4 period map (h = 1e-7) reproduce the
#    AD-quality multipliers to ~1e-4 — ample, since they are ≈1e-2. The
#    KEEP-only companion `scratch/_floquet_keep.jl` does this and stays ~300 MB.
#
# 8. FIXED IN THIS REVISION: `make_bvp` ignored its `disc` argument (the
#    shooting path silently used the collocation discretizer) and
#    `sheet_jumps` indexed `br.sol` (a single saved solution, not the path).
#
# 9. TWO CYCLES AT ONE v_ref (a consequence of the fold). The family is a
#    CURVE (v_ref, tf), not a graph tf(v_ref): it turns at
#    v_ref* = 2.45847 (tf_f = 14.276 s), so v_ref = 2.47 cuts it twice — a
#    strongly stable "short" cycle (tf = 11.37 s, |mu| = {1, 0.294, 9e-4, ~0})
#    and a "long" cycle (tf = 28.82 s) that is strongly UNSTABLE at v_ref=2.47
#    (dominant |mu| ≈ 3.5e16, Floquet exponent ≈ 1.3 /s). The two are the
#    saddle/short pair of the saddle-node: the long cycle's multiplier passes
#    through +1 only exactly AT the fold and moves far from it at t fixed
#    v_ref = 2.47, so "long ⇒ mu ≈ +1" is true only in the limit v → v_ref*.
#    No separate arclength driver is needed: ContinuationPar defaults to
#    save_sol_every_step = 1, so the single multiple-shooting continuation
#    already stores BOTH crossings in br.sol; select those with
#    |s.p - v_ref*| < tol, sort by raw_x(s.x)[5], take the candidates nearest
#    v_ref on each sheet, then Newton-refine each at v_ref*. Collocation does
#    NOT pass the fold (it stalls at v_ref ≈ 2.4868, tf ≈ 10.56 s), so it
#    yields only the short cycle — the "try both discretizations" check.
#    MULTIPLIER CAVEAT: a naive single-shooting period-map FD CANNOT measure
#    the long cycle — it is so unstable over 28.8 s that the integrated return
#    misses closure by ~1e1 (garbage eigenvalues ~1e8). Compute it as the
#    PRODUCT of the shooting sub-interval flow Jacobians (segment-wise
#    monodromy) — `monodromy`/`monodromy_max` above, applied to both crossings
#    in the two-cycles block. The historical extraction/plotting scripts are
#    archived under `scratch/` (`_two_cycles_bk.jl`, `_two_cycles_floquet.jl`).
#    The long-cycle tf/|mu| quoted above is the archived run's; which long-sheet
#    point the in-script selection refines depends on the continuation budget.
#
# 10. SHOOTING JACOBIAN. BK's `AutoDiffDense` differentiates the whole residual
#     through the adaptive solve, nesting ForwardDiff inside `dynamics`, and
#     effectively never finishes. The manual block Jacobian above is used
#     instead: matching blocks are segment transition matrices, the boundary
#     row is ∂g/∂u, and the trailing auxiliary row is 1. Blocks come from
#     ForwardDiff (`ManualJacFwd`) or central differences (`ManualJacFD`); the
#     two agree to ~1e-7 (validated against a full FD of the residual) and
#     `ManualJacFD` is the AD-free fallback if the nested-AD memory ever bites.
#     Shooting also needs a collocation-seeded warm start: the raw optimized
#     sample sits too far off the default-parameter problem for the damped
#     pre-solve, whereas a collocation cycle at the same parameters converges
#     and seeds the shooting points (full shooting run ≈ 2.8 min).
#
# The core of it is one saddle-node (fold) of limit cycles at low wind:
# Optimization: pick the operating cycle at v_ref = 9 m/s + 3 parameters (r, I_eq, torque_slope).
# Continuation then sweeps v_ref. The stable "short" cycle runs tf ≈ 0.63 s @ 9 → 0.12 s @ 20 m/s.
# Low wind: the family folds at v_ref* ≈ 2.458 m/s (tf ≈ 14.3 s), and the second, long-period branch beyond it is a strongly unstable saddle cycle. That fold is the "hard to pass" part numerically (collocation stalls just above it; multiple shooting gets through).
# But that's not all that's in here — three secondary things:
# 1. Fold tangle near 9.2–9.8 m/s: several folds crowd sheets, so the Newton corrector can jump sheets.
# 2. Coexisting mirror cycles (symmetry α,dτ → −α,−dτ), found by Poincaré sampling, not by a single continuation.
# 3. Equilibria: 8 (4 mirror pairs), all unstable everywhere; no clean Hopf generates the cycle in the reachable range.
# So: the headline phenomenon is the low-wind fold + unstable long branch; the tangle, the mirror pair, and the all-unstable equilibria are the side structure.
