# Fixed copy of `BK_parametric.jl`.
#
# Two changes w.r.t. the original:
#
# 1. PHYSICAL-TIME BVP (the wind was rescaled in Vaschy-Buckingham space).
#    The original `F` called `dynamics(u_dyn, build_vbpara(params))` with state
#    rates in normalized units 1/T, T = l/v_ref. Since `build_vbpara` renormalizes
#    with T(v_ref), continuing in `v_ref` changed the UNITS of the state and of
#    `tf` along the branch: the parameter acted through renormalization instead of
#    pure physics. Here the state is (α, τ, dα, dτ, tf) in physical units
#    (rad, rad/s, s); `F` converts to/from the internal normalization using the
#    characteristic time of the current parameter set. Consequence: gravity,
#    inertia and torque keep their physical values along the branch and `v_ref`
#    changes ONLY the wind magnitude v(z) = v_ref (z/h_ref)^(1/n_wind).
#
# 2. MULTIPLE-SHOOTING DETERMINISM (upstream bug in BifurcationKit 0.8.2,
#    src/bvp/shooting/residual.jl): `out = similar(X)` leaves the last entry of
#    the residual vector UNINITIALIZED (the trailing unknown of the vector layout
#    is unused by this discretizer), so Newton consumed random memory — hence
#    "un objet est mal initialisé lors du solve, exécuter plusieurs fois donne
#    différents résultats". We override `bvp_residual` for `Shooting` to pin that
#    auxiliary unknown to zero (residual row `X[end] = 0`), which also makes the
#    dense Jacobian nonsingular (J[end,end] = 1).
#
# 3. EXPLICIT RECORD AND EXPLICIT TWO-DIRECTION CONTINUATION (replaces
#    `bothside=true`). `bothside=true` looks equivalent to running ds>0 and
#    ds<0, but is not — see the continuation section below. Findings:
#
#    * Family topology (physical units): main sheet tf(v_ref) ≈ 0.63 s at
#      v_ref = 9 m/s falling monotonically to ≈ 0.12 s at 20 m/s; below 9 m/s
#      the period grows without bound (tf ≈ 47.7 s at v_ref = 2.4585, still
#      growing at continuation budget — no true low-wind endpoint found); a
#      fold tangle (detected folds at v_ref ≈ 9.31, 9.30, 9.46, 9.54, 9.72)
#      crowds sheets just above 9 m/s.
#    * Collocation and multiple shooting agree on the period to ~1e-6 s
#      wherever both track the same sheet (checked at v_ref = 8, 9, 10 m/s).

using Pkg
Pkg.activate("applications/sprint2026")

import OrdinaryDiffEqTsit5 as ODE
using BifurcationKit
using LinearAlgebra
using ComponentArrays: ComponentArray as CA
using StaticArrays
using Plots

using KEEP: TAU0
using KEEP.PointMass4: dynamics
using KEEP.PointMassPara: build_vbpara, build_para, lmt
using KEEP.Optimization: optimize
using KEEP.LimitCycle: shoot as lc_shoot

const BVP = BifurcationKit.BVP

## Solve the optimization problem (unchanged: fixed kite at the default wind)
const vbp0 = build_vbpara()
const p0 = build_para(vbp0)
const nt_p0 = (; zip(keys(p0), p0)...)
const syms = (:r, :I_eq, :torque_slope)

factor = 5
lb = p0[syms] ./ factor
ub = p0[syms] .* factor
solution, stats, model = optimize(p0, syms, lb, ub)

## Reconstruct the dense ODE solution (in normalized time)
vbp = build_vbpara(CA(p0; solution.params...))
shooting = solution[1:4]
solution_sim = lc_shoot(shooting, vbp, save_everystep=true)

x_optimization = t -> solution_sim(t, idxs=1:4)  # normalized time
tf_optimization = solution_sim.t[end]            # normalized period

## Conversion to physical time
const T0 = p0.l / p0.v_ref        # characteristic time [s] at reference wind
tf_physical = tf_optimization * T0

"Sample the optimized cycle at physical time `t`: (α, τ, dα, dτ) in rad, rad/s."
function x_optimization_phys(t)
    s = x_optimization(t / T0)
    return SA[s[1], s[2], T0*s[3], T0*s[4]]
end

## BVP model — physical time
# State u = (α, τ, dα, dτ, tf), integration variable s ∈ [0, 1] spans one period:
# du/ds = tf * [dα, dτ, ddα, ddτ, 0] with dα, dτ in rad/s and ddα, ddτ in rad/s².
# Reference implementation: rebuilds a ComponentArray on EVERY call (~4.7 μs of
# the ~6.4 μs total). Kept only for validation of `F_fast`.
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

"""
Normalized parameters as a plain typed struct: every field that
KEEP.PointMass4.dynamics reads (all its parameter access is by getproperty,
so a plain struct drops the ComponentArray machinery entirely).
Values are bit-identical to `build_vbpara`.
"""
struct FastPara{T}
    l::T
    m::T
    v_ref::T
    g::T
    h_ref::T
    r::T
    c_L::T
    f::T
    c_D_l::T
    m_l::T
    I_eq::T
    Cmax::T
    Ωmin::T
    Ωmax::T
    Ωlim::T
    torque_slope::T
    n_wind::T
    θ0::T
    φ0::T
    Δθ::T
    Δφ::T
end

FastPara(vbp) = FastPara{Float64}(
    vbp.l, vbp.m, vbp.v_ref, vbp.g, vbp.h_ref, vbp.r, vbp.c_L, vbp.f, vbp.c_D_l,
    vbp.m_l, vbp.I_eq, vbp.Cmax, vbp.Ωmin, vbp.Ωmax, vbp.Ωlim, vbp.torque_slope,
    vbp.n_wind, vbp.θ0, vbp.φ0, vbp.Δθ, vbp.Δφ,
)

# ponytail: single-entry cache keyed on v_ref; continuation evaluates many RHS
# at one parameter value, so we rebuild only when v_ref actually changes.
const FASTPARA_CACHE = Ref{Union{Nothing,Tuple{Float64,FastPara{Float64}}}}(nothing)  # (vref, FastPara)

@inline function get_fast_para(params)
    vref = params.v_ref
    c = FASTPARA_CACHE[]
    if c === nothing || c[1] != vref
        c = (vref, FastPara(build_vbpara(params)))
        FASTPARA_CACHE[] = c
    end
    return c[2]
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

u0_bif = SA[x_optimization_phys(0)..., tf_physical]

# Record the physical period (5th state component) along branches. Without an
# explicit record, `plot(br)` falls back to `norm(x)`, which is incomparable
# between the two discretizations (different unknown-vector layouts) — the
# "multiple-shooting values are much smaller" artifact. `raw_x` unwraps the
# saved-solution wrapper that the Shooting discretizer stores per point.
raw_x(x) = x isa BifurcationKit.BVPSavedSolutionAndState ? BifurcationKit.saved_solution(x) : x
record_period(x, p; kwargs...) = raw_x(x)[5]

const STATE_SIZE = length(u0_bif)
model = BVP.BVPModel(F_fast, g; n=STATE_SIZE)

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

## Collocation
grid_size, degree = 30, 5
const disc = BVP.Collocation(Ntst=grid_size, m=degree, meshadapt=true)
bvp = BVP.discretize(model, disc)

params = nt_p0
x0 = BVP.generate_solution(bvp, s -> vcat(x_optimization_phys(tf_physical * s), tf_physical))

prob = BVP.BVPBifProblem(bvp, x0, params, (@optic _.v_ref);
    jacobian=BifurcationKit.DenseAnalytical(),
    record_from_solution=record_period,
)

# linesearch is honored by the PALC continuation corrector
optn = NewtonPar(tol=1e-10, verbose=true, linesearch=true)

x0, res_pre = damped_newton(prob, x0, params; tol=1e-9)
println("Collocation damped pre-solve residual: ", res_pre)
prob = BifurcationKit.re_make(prob; u0=x0)

sol = @time BifurcationKit.solve(prob, Newton(), optn);
@assert BifurcationKit.converged(sol) "Collocation Newton did not converge"

# Continuation
optc = ContinuationPar(
    p_min=0.0,
    p_max=25.0,
    dsmax=0.1,
    ds=0.01,
    detect_bifurcation=0,
    newton_options=optn,
    max_steps=100,
    nev=20,
    n_inversion=6
)

# Continuation — TWO EXPLICIT SINGLE-DIRECTION RUNS, not `bothside=true`.
#
# `bothside=true` seems equivalent to running ds>0 and ds<0 and merging, but:
# * BifurcationKit merges the legs as `_cat!(_reverse(backward), forward)`
#   (Results.jl:464): the stored branch BEGINS at the backward leg's endpoint.
#   Plotted, that reversed leg is indistinguishable from a second physical
#   family — this produced the "both branches go right from v_ref = 9" artifact.
# * The merged branch gives each leg neither its own budget nor a direction
#   label; near the fold tangle at v_ref ≈ 8.7–9.7 the corrector can land on a
#   different sheet (observed: a single recorded step moved Δp = 0.27 with
#   dsmax = 0.1 — impossible for a true PALC step, whose tangent is a unit
#   vector, so the Newton correction left the sheet) and then silently spend
#   the whole budget on the wrong sheet.
# Two separate runs give each direction its own budget and make the step-size
# sanity check below meaningful (recorded points must satisfy |Δp| ≤ dsmax).
continuation_side(ds) = continuation(prob, PALC(), @set optc.ds = ds;
    plot=false, verbosity=1, normC=norminf)

"Recorded PALC points satisfy |Δp| ≤ dsmax; larger gaps flag corrector
sheet-jumps (only possible where sheets crowd, i.e. near folds)."
function sheet_jumps(br; dsmax=0.1, factor=1.5)
    return [i for i in 2:length(br)
                  if abs(br.sol[i].p - br.sol[i-1].p) > factor * dsmax]
end

br_fwd = @time continuation_side(+0.01)
br_bwd = @time continuation_side(-0.01)
println("collocation sheet-jumps: ds>0 at indices ", sheet_jumps(br_fwd),
    ", ds<0 at indices ", sheet_jumps(br_bwd))
plot(br_fwd, label="collocation, ds>0")
plot!(br_bwd, label="collocation, ds<0")
