# OC_controlled_kite_v3.jl — hybrid: Symbolics (Lagrangian derivatives) + plain numeric aero
#
# v2 (runtime ForwardDiff inside `eom`) failed OptimalControl's sparsity detection:
#   MethodError: *(GradientTracer, Dual{...}) is ambiguous
# because SparseConnectivityTracer cannot trace through ForwardDiff-inside-the-RHS.
#
# v3 design:
#  * The EOM is M(q)·q̈ = Q(q,q̇,u) + ∂L/∂q − (∂²L/∂q̇∂q)·q̇.
#  * M, ∂L/∂q, (∂²L/∂q̇∂q)·q̇ and the virtual-work map J_CG come ONLY from the
#    Lagrangian (kinematics/inertia: cheap trig), derived once with Symbolics.jl
#    and codegenned via build_function(force_SA=true) → straight-line StaticArrays.
#  * The generalized forces Q come from the aero model, written as a PLAIN numeric
#    function (smooth tanh stall, sqrt, atan — no comparisons, no `\`, no ForwardDiff).
#    Q is never differentiated, so it never needs to be symbolic.
#  * Result: `eom` is a straight-line function. SCT / ForwardDiff / ReverseDiff all
#    trace through it trivially. No runtime AD, no huge symbolic aero expressions.
#
# Minimal-first: the OCP section below uses grid_size=8 and loose Ipopt tolerances
# so this script completes quickly and validates the pipeline. Bump grid_size /
# tighten tolerances (or use grid continuation) once it converges.
#
# Run: julia --project=applications/controlled OC_controlled_kite_v3.jl

using Pkg
Pkg.activate("applications/controlled")

using StaticArrays
using LinearAlgebra
using Symbolics
using ForwardDiff
using OrdinaryDiffEqTsit5
using Plots

using OptimalControl
using NLPModelsIpopt
using JLD2

# =========================================================
# Parameters
# =========================================================

function kite_params(; larm=2.0, lcg=20.0, mass=5.0, inertia=(25.0, 0.5, 2.5),
    inertiaarm=8.0, g=9.81, kgenerator=0.0, cgenerator=1.0, klines=0.0,
    S=2.5, rho=1.225, delta=30.0 * pi / 180)
    panels = @SMatrix [0.1 0.1; 1.25 -1.25; 0.0 0.0]
    return (;
        larm=larm, lcg=lcg, mass=mass, inertia=SA[inertia[1], inertia[2], inertia[3]],
        inertiaarm=inertiaarm,
        g=g, kgenerator=kgenerator, cgenerator=cgenerator, klines=klines,
        aero=(S=S, rho=rho, panels=panels, delta=delta),
    )
end

const pars = kite_params()

# =========================================================
# Symbolic derivation — Lagrangian part only (cheap trig)
# =========================================================

@variables q[1:4] dq[1:4]

# CG position and velocity (explicit chain rule)
pos_CG(q) = [pars.lcg * sin(q[2]) * cos(q[3]) + pars.larm * cos(q[1]),
    pars.lcg * sin(q[2]) * sin(q[3]) + pars.larm * sin(q[1]),
    pars.lcg * cos(q[2])]

function vel_CG(q, dq)
    l, b = pars.lcg, pars.larm
    cq, sq = cos.(q), sin.(q)
    return [
        -b * sq[1] * dq[1] + l * cq[2] * cq[3] * dq[2] - l * sq[2] * sq[3] * dq[3],
        b * cq[1] * dq[1] + l * cq[2] * sq[3] * dq[2] + l * sq[2] * cq[3] * dq[3],
        -l * sq[2] * dq[2],
    ]
end

# body angular velocity ω = W(q)·dq, written out explicitly
function omega(q, dq)
    cq, sq = cos.(q), sin.(q)
    return [
        dq[4] + cq[2] * dq[3],
        sq[4] * sq[2] * dq[3] + cq[4] * dq[2],
        cq[4] * sq[2] * dq[3] - sq[4] * dq[2],
    ]
end

function lagrangian(q, dq)
    vCG = vel_CG(q, dq)
    ω = omega(q, dq)
    T = 0.5 * pars.mass * (vCG[1]^2 + vCG[2]^2 + vCG[3]^2) +
        0.5 * (pars.inertia[1] * ω[1]^2 + pars.inertia[2] * ω[2]^2 + pars.inertia[3] * ω[3]^2) +
        0.5 * pars.inertiaarm * dq[1]^2
    V = pars.mass * pars.g * pos_CG(q)[3]
    return T - V
end

L = lagrangian(q, dq)
momentum = Symbolics.gradient(L, dq)          # ∂L/∂q̇
M_sym = Symbolics.jacobian(momentum, dq)      # M(q) = ∂²L/∂q̇²
dLdq_sym = Symbolics.gradient(L, q)           # ∂L/∂q
cross_sym = Symbolics.jacobian(momentum, q) * dq   # (∂²L/∂q̇∂q)·q̇
J_CG_sym = Symbolics.jacobian(pos_CG(q), q)   # 3x4 virtual-work map

# simplify keeps the generated code small (big compile-time win, faster duals)
M_sym = Symbolics.simplify.(M_sym)
dLdq_sym = Symbolics.simplify.(dLdq_sym)
cross_sym = Symbolics.simplify.(cross_sym)
J_CG_sym = Symbolics.simplify.(J_CG_sym)

println("symbolic: M=", size(M_sym), " dLdq=", size(dLdq_sym),
    " cross=", size(cross_sym), " J_CG=", size(J_CG_sym))

# Codegen straight-line StaticArray functions
Mf = Symbolics.build_function(M_sym, q; force_SA=true, expression=Val{false})[1]
dLdqF = Symbolics.build_function(dLdq_sym, q, dq; force_SA=true, expression=Val{false})[1]
crossF = Symbolics.build_function(cross_sym, q, dq; force_SA=true, expression=Val{false})[1]
JCGF = Symbolics.build_function(J_CG_sym, q; force_SA=true, expression=Val{false})[1]

# =========================================================
# Aerodynamics — plain numeric, branch-free, never differentiated
# =========================================================

wind(z) = 20.0

# smooth tanh window: ~1 inside |α5|≤30°, ~0 outside
stall_factor(alpha) = 0.25 * (1.0 + tanh(30.0 * (alpha + 5.0 * pi / 180.0 + 30.0 * pi / 180.0))) *
                      (1.0 - tanh(30.0 * (alpha + 5.0 * pi / 180.0 - 30.0 * pi / 180.0)))

function Caero(alpha)
    α5 = alpha + 5.0 * pi / 180.0
    return SVector(1.0 * α5 * stall_factor(alpha), 0.2 + 0.1 * α5^2, 0.0)
end

cross3(a, b) = SA[a[2]*b[3]-a[3]*b[2], a[3]*b[1]-a[1]*b[3], a[1]*b[2]-a[2]*b[1]]

function body_axes(q)
    sq, cq = sin.(q), cos.(q)
    Ihat = SVector(sq[2] * cq[3], sq[2] * sq[3], cq[2])
    Jhat = SVector(-cq[4] * sq[3] - sq[4] * cq[2] * cq[3],
        cq[4] * cq[3] - sq[4] * cq[2] * sq[3],
        sq[4] * sq[2])
    Khat = SVector(sq[4] * sq[3] - cq[4] * cq[2] * cq[3],
        -sq[4] * cq[3] - cq[4] * cq[2] * sq[3],
        cq[4] * sq[2])
    return Ihat, Jhat, Khat
end

# returns (F, M) forces/moment in the CG frame
function aerodynamics(q, dq, u)
    a = pars.aero
    Ihat, Jhat, Khat = body_axes(q)
    vCG = vel_CG(q, dq)
    vCGrel = vCG - SVector(wind(pos_CG(q)[3]), 0.0, 0.0)
    vCGb = SVector(dot(Ihat, vCGrel), dot(Jhat, vCGrel), dot(Khat, vCGrel))
    ω = omega(q, dq)

    Rp = a.panels
    Rp1, Rp2 = Rp[:, 1], Rp[:, 2]
    Vair1 = -(vCGb + cross3(ω, Rp1))
    Vair2 = -(vCGb + cross3(ω, Rp2))

    cdelta, sdelta = cos(a.delta), sin(a.delta)
    w1 = SVector(-sdelta, cdelta, 0.0)
    w2 = SVector(-sdelta, -cdelta, 0.0)
    l1 = SVector(cdelta, sdelta, 0.0)
    l2 = SVector(cdelta, -sdelta, 0.0)

    Vplane1 = Vair1 - w1 * dot(w1, Vair1)
    Vplane2 = Vair2 - w2 * dot(w2, Vair2)
    Vp1n = sqrt(sum(abs2, Vplane1))
    Vp2n = sqrt(sum(abs2, Vplane2))
    Vs1 = dot(w1, Vair1)
    Vs2 = dot(w2, Vair2)

    dragDir1 = Vplane1 / Vp1n
    dragDir2 = Vplane2 / Vp2n
    liftDir1 = -cross3(w1, dragDir1)
    liftDir2 = cross3(w2, dragDir2)

    i1 = atan(dot(l1, Vplane1), -Vplane1[3]) + u[1]
    i2 = atan(dot(l2, Vplane2), -Vplane2[3]) + u[2]
    Ca1 = Caero(i1)
    Ca2 = Caero(i2)

    c = 0.5 * a.rho * a.S
    F1 = c * (Vp1n^2 * (Ca1[1] * liftDir1 + Ca1[2] * dragDir1) + Vs1^2 * Ca1[3] * w1)
    F2 = c * (Vp2n^2 * (Ca2[1] * liftDir2 + Ca2[2] * dragDir2) + Vs2^2 * Ca2[3] * w2)
    F = F1 + F2
    M = cross3(Rp1, F1) + cross3(Rp2, F2)
    return F, M
end

# =========================================================
# EOM: straight-line, no ForwardDiff inside
# =========================================================

function cholesky_solve(M, b)
    l11 = sqrt(M[1, 1])
    l21 = M[2, 1] / l11
    l31 = M[3, 1] / l11
    l41 = M[4, 1] / l11
    l22 = sqrt(M[2, 2] - l21^2)
    l32 = (M[3, 2] - l31 * l21) / l22
    l42 = (M[4, 2] - l41 * l21) / l22
    l33 = sqrt(M[3, 3] - l31^2 - l32^2)
    l43 = (M[4, 3] - l41 * l31 - l42 * l32) / l33
    l44 = sqrt(M[4, 4] - l41^2 - l42^2 - l43^2)
    y1 = b[1] / l11
    y2 = (b[2] - l21 * y1) / l22
    y3 = (b[3] - l31 * y1 - l32 * y2) / l33
    y4 = (b[4] - l41 * y1 - l42 * y2 - l43 * y3) / l44
    x4 = y4 / l44
    x3 = (y3 - l43 * x4) / l33
    x2 = (y2 - l32 * x3 - l42 * x4) / l22
    x1 = (y1 - l21 * x2 - l31 * x3 - l41 * x4) / l11
    return SVector(x1, x2, x3, x4)
end

function eom(X, u)
    q = SVector{4}(X[1:4])
    dq = SVector{4}(X[5:8])

    # Generalized forces (numeric, branch-free): aero + generator + lines
    F, Ma = aerodynamics(q, dq, u)
    cq, sq = cos.(q), sin.(q)
    Wnum = @SMatrix [
        0.0 0.0 cq[2] 1.0;
        0.0 cq[4] sq[4] * sq[2] 0.0;
        0.0 -sq[4] cq[4] * sq[2] 0.0
    ]
    Qaero = JCGF(q)' * F + Wnum' * Ma
    Qgen = SVector(pars.kgenerator * q[1] + pars.cgenerator * dq[1], 0.0, 0.0, 0.0)
    Qline = SVector(0.0, 0.0, 0.0, pars.klines * q[4])
    Q = Qaero - Qgen + Qline

    # M·q̈ = Q + ∂L/∂q − (∂²L/∂q̇∂q)·q̇
    b = Q + dLdqF(q, dq) - crossF(q, dq)
    M = Mf(q)
    ddq = cholesky_solve(SMatrix{4,4}(M + 1e-9 * I), SVector{4}(b))

    return SVector{8}(dq[1], dq[2], dq[3], dq[4], ddq[1], ddq[2], ddq[3], ddq[4])
end

f(X, u) = eom(X, u)

# =========================================================
# Minimal sanity checks
# =========================================================
let X = SA[0.3, 0.5, 0.2, 0.6, 0.1, 0.2, 0.3, 0.4], u = SA[0.1, 0.0]
    y = f(X, u)
    @assert y isa SVector{8,Float64}
    J = ForwardDiff.jacobian(x -> f(x, u), X)
    @assert size(J) == (8, 8)
    g = ForwardDiff.gradient(x -> sum(f(x, u)), X)
    @assert length(g) == 8
    println("sanity: eom OK, jacobian=", size(J), ", ddq = ", y[5:8])
end
# =========================================================
# Helper functions for grid_size continuation
# =========================================================

"""
    cache_filepath(cache_dir, prefix, grid_size)

Returns file path prefix (without extension) for a given cache entry.
"""
cache_filepath(cache_dir::String, prefix::String, grid_size::Int) = joinpath(cache_dir, "$(prefix)_grid_$(grid_size)")

"""
    solve_cached(ocp; grid_size, init, cache_dir, prefix, cache=:exact, solve_options...)

Handles solution retrieval or execution for a single grid size according to `cache`:
- `:exact`: Directly loads the cached solution without solving; throws an error if missing.
- `:no` / `:latest`: Solves the OCP and exports the solution to disk.
"""
function solve_cached(ocp;
    grid_size::Int,
    init,
    cache_dir::String,
    prefix::String,
    cache::Symbol=:exact,
    solve_options...)
    fpath = cache_filepath(cache_dir, prefix, grid_size)
    jld2_file = fpath * ".jld2"

    if cache == :exact
        isfile(jld2_file) || error("Cached solution not found for grid = $grid_size: $jld2_file")
        @info "Loading cached solution: grid = $grid_size ($jld2_file)"
        return import_ocp_solution(ocp; filename=fpath)
    end

    @info "Solving OCP: grid = $grid_size..."
    sol = solve(ocp; init=init, grid_size=grid_size, solve_options...)
    export_ocp_solution(sol; filename=fpath)
    return sol
end

"""
    resolve_cache_plan(ocp, grid_schedule; init, cache_dir, prefix, cache)

Determines the initial guess and schedule based on the cache mode:
- `:no`: Ignores disk cache; returns `(init, grid_schedule)`.
- `:exact`: Returns `(init, grid_schedule)` for direct loading across all target grids.
- `:latest`: Finds the highest-grid cached solution matching `prefix` to use as `init`.
"""
function resolve_cache_plan(ocp, grid_schedule; init, cache_dir::String, prefix::String, cache::Symbol)
    cache in (:no, :exact, :latest) || throw(ArgumentError("cache must be :no, :exact, or :latest (got :$cache)"))

    if cache == :no || !isdir(cache_dir)
        return init, collect(grid_schedule)
    end

    if cache == :latest
        # Match any cache file with the given prefix and extract the highest grid number
        pattern = Regex("^$(prefix)_grid_(\\d+)\\.jld2\$")
        matched_grids = Int[]
        for f in readdir(cache_dir)
            m = match(pattern, f)
            m !== nothing && push!(matched_grids, parse(Int, m.captures[1]))
        end

        if !isempty(matched_grids)
            latest_grid = maximum(matched_grids)
            fpath = cache_filepath(cache_dir, prefix, latest_grid)
            @info "Found latest cached solution: grid = $latest_grid ($fpath.jld2)"
            cached_sol = import_ocp_solution(ocp; filename=fpath)
            return cached_sol, collect(grid_schedule)
        end
        return init, collect(grid_schedule)
    end

    # :exact mode verifies all required files exist before execution
    return init, collect(grid_schedule)
end

"""
    run_grid_homotopy(ocp, grid_schedule; init, cache=:exact, cache_dir="applications/controlled/solutions", prefix="ocp_half", solve_options...)

Executes grid continuation over `grid_schedule` (e.g. `(8, 20, 50)`).

# Cache Modes (`cache::Symbol`)
- `:no`: Optimizes through `grid_schedule` starting from `init`, ignoring existing cache files.
- `:exact`: Performs no optimization; loads each cached grid solution directly from disk.
- `:latest`: Uses the highest available cached grid solution matching `prefix` as `init`, then optimizes through `grid_schedule`.
"""
function run_grid_homotopy(ocp, grid_schedule; init,
    cache::Symbol=:exact,
    cache_dir::String="applications/controlled/solutions",
    prefix::String="ocp_half",
    solve_options...)
    mkpath(cache_dir)

    start_sol, remaining_grids = resolve_cache_plan(
        ocp, grid_schedule;
        init=init, cache_dir=cache_dir, prefix=prefix, cache=cache
    )

    if isempty(remaining_grids)
        println("✓ Latest solution already cached at finest grid ($((grid_schedule)[end])) — Objective: ", round(objective(start_sol), digits=4))
        return start_sol
    end

    # Sequential continuation across remaining grids
    sol_final = foldl(remaining_grids; init=start_sol) do warm_start, N
        sol_current = solve_cached(ocp;
            grid_size=N,
            init=warm_start,
            cache_dir=cache_dir,
            prefix=prefix,
            cache=cache,
            solve_options...
        )
        println("✓ Grid $N converged — Objective: ", round(objective(sol_current), digits=4))
        sol_current
    end

    return sol_final
end

# =========================================================
# Minimal OCP (grid_size=8, loose tolerances) — validate the pipeline fast
# =========================================================
TF_MIN, TF_MAX = 0.5, 20.0

# generator torque = -(k·α + c·α̇)  ⇒  generated power = torque·α̇ = α̇·(k·α + c·α̇)
generated_power(x) = x[5] * (pars.kgenerator * x[1] + pars.cgenerator * x[5])

ocp = @def begin
    p = (tf, p_) ∈ R², variable
    t ∈ [0, tf], time
    x = (α, θ, φ, β, α̇, θ̇, φ̇, β̇) ∈ R⁸, state
    u = (uₗ, uᵣ) ∈ R², control

    p_ == 0

    tf >= TF_MIN
    tf <= TF_MAX

    -30 * pi / 180 ≤ uₗ(t) ≤ 30 * pi / 180     # physical aileron limits
    -30 * pi / 180 ≤ uᵣ(t) ≤ 30 * pi / 180

    α(0) == 0  # Phase fixing
    α̇(0) >= 0  # Arm direction at t=0 (left)
    φ̇(0) >= 0  # Kite horizontal direction at t=0 (left -> pulling arm)
    θ̇(0) <= 0  # Kite vertical direction at t=0 (climbing)
    θ(t) >= 0  # Kite is always above ground, can be more complicated to set a minimal height
    x(0) - x(tf) == zeros(8)

    ẋ(t) == f(x(t), u(t))

    ∫(-generated_power(x(t)) + 1e-3 * sum(abs2, u(t))) → min
end

## Generate an initial guess

using KEEP.PointMassPara: build_vbpara
using KEEP.PointMass4: τ_to_θφ
using KEEP.LimitCycle: compute_limit_cycle

vbp = build_vbpara()
const lc = compute_limit_cycle(vbp, sense=(+), save_everystep=true)

fα(t) = lc(t)[1]
fθ(t) = τ_to_θφ(lc(t)[2], vbp)[1]
fφ(t) = τ_to_θφ(lc(t)[2], vbp)[2]

fpos(t) = begin
    LU = vbp.l
    a, θl, φl = fα(t), fθ(t), fφ(t)
    LU .* [cos(a) + vbp.r * sin(θl) * cos(φl),
        sin(a) + vbp.r * sin(θl) * sin(φl),
        vbp.r * cos(θl)]
end

fβ(t) = begin
    v = ForwardDiff.derivative(fpos, t)
    θl, φl = fθ(t), fφ(t)
    J0 = [-sin(φl), cos(φl), 0.0] # Side axis at β=0
    K0 = [-cos(θl) * cos(φl), -cos(θl) * sin(φl), sin(θl)] # Up axis at β=0
    atan(-(v[1] * J0[1] + v[2] * J0[2] + v[3] * J0[3]),
        (v[1] * K0[1] + v[2] * K0[2] + v[3] * K0[3]))
end

init = @init ocp begin
    # One symbol only on lhs, `θ, φ = τ_to_θφ(x[2], vbp)` is not allowed for example
    tf := lc.t[end]

    α(t) := fα(t)
    θ(t) := fθ(t)
    φ(t) := fφ(t)
    β(t) := fβ(t)
    α̇(t) := ForwardDiff.derivative(fα, t)
    θ̇(t) := ForwardDiff.derivative(fθ, t)
    φ̇(t) := ForwardDiff.derivative(fφ, t)
    β̇(t) := ForwardDiff.derivative(fβ, t)
end

# tf_init = 3
# init = @init ocp begin
#     tf := tf_init
#     # genuinely periodic guess: α oscillates (drives the generator), others steady
#     A = 0.5
#     ω = 2π / tf_init
#     x(t) := SA[
#         A*sin(ω*t),
#         pi/6,
#         0.0,
#         0.0,
#         A*ω*cos(ω*t),
#         0.0,
#         0.0,
#         0.0
#     ]
#     u(t) := SVector(0.0, 0.0)
# end

solve_options = (
    backend=:default,
    print_level=5,
    linear_solver="mumps",
    tol=1e-4,
    # maxiter=50
    # nlp_scaling_method="gradient-based"

)

sol = run_grid_homotopy(
    ocp,
    (8, 20, 50);
    init=init,
    cache=:latest,                     # :no | :exact | :latest
    cache_dir="applications/controlled/solutions",
    prefix="kite_lemniscate",
    solve_options...
)


plot(sol)

#= TODO

 - [x] Initialize with PM2 trajectory
 - [x] Add phase-fixing
 - [x] Add state constraints : going right then left or the inverse
 - [x] simulate only half the trajectory to enforce symmetry

 - [x] GPU (ExaModel+MadNLP): Metal not supported
 - [x] AppleAccelerate: seems marginally slower

 - [] load from cache: if ocp is the same as the one in cache: return directly, otherwise use as init for solve
 - [] Give the function parameters
 - [] Translate from PM2 params to controlled params
 - [] Transform state from 2 <-> 4
 - [] Take parameters p2 and p4
 - [] visualize
 - [] add tests matlab-v1 et v1-v3
=#

# =========================================================
# TLDR
# =========================================================
# 1. Runtime ForwardDiff inside eom cannot work: SCT sparsity tracing can't nest
#    inside ForwardDiff (v2's MethodError). Instead: symbolic derivatives of the
#    LAGRANGIAN ONLY (M, ∂L/∂q, cross-term, J_CG) are codegenned once; the aero
#    forces enter through Q as plain branch-free numeric code (never differentiated).
# 2. eom is a straight-line function → traceable by SCT/ForwardDiff/ReverseDiff.
# 3. Aero polar uses a smooth tanh window (no hard branch / abs).
# 4. backend=:default avoids the ReverseDiff-over-ForwardDiff DualMismatchError.
# To go full-scale: raise grid_size (with warm-start continuation) and tighten
# Ipopt tolerances; if MUMPS struggles, try linear_solver="ma57".

##
# =========================================================
# 3D Visualization & Animation Functions
# =========================================================

"""
    plot_kite_3d(sol; N_pts=300, n_arrows=6, n_tethers=4, camera=(40, 25), size=(900, 750), filename=nothing)

Generates a 3D plot of the kite trajectory showing:
- 3D CG trajectory with directional arrows
- Projections on bounding floor (XY), back wall (XZ), and side wall (YZ)
- Arm trajectory and sample tethers
"""
function plot_kite_3d(sol;
    N_pts::Int=300,
    n_arrows::Int=6,
    n_tethers::Int=4,
    camera=(40, 25),
    size=(900, 750),
    filename::Union{String,Nothing}=nothing)

    tf_sol = time_grid(sol)[end]
    t_eval = range(0, tf_sol, length=N_pts)

    q_eval = [state(sol)(t)[1:4] for t in t_eval]
    pos_eval = [pos_CG(q) for q in q_eval]
    arm_eval = [SVector(pars.larm * cos(q[1]), pars.larm * sin(q[1]), 0.0) for q in q_eval]

    X_cg = [p[1] for p in pos_eval]
    Y_cg = [p[2] for p in pos_eval]
    Z_cg = [p[3] for p in pos_eval]

    X_arm = [p[1] for p in arm_eval]
    Y_arm = [p[2] for p in arm_eval]
    Z_arm = [p[3] for p in arm_eval]

    # Compute bounding limits
    pad = 2.0
    all_x = [X_cg; X_arm; 0.0]
    all_y = [Y_cg; Y_arm; 0.0]
    all_z = [Z_cg; 0.0]

    xmin, xmax = minimum(all_x) - pad, maximum(all_x) + pad
    ymin, ymax = minimum(all_y) - pad, maximum(all_y) + pad
    zmin, zmax = 0.0, maximum(all_z) + pad

    p3d = plot(
        title="Kite 3D Periodic Limit Cycle",
        xlabel="X [m]", ylabel="Y [m]", zlabel="Z [m]",
        xlims=(xmin, xmax), ylims=(ymin, ymax), zlims=(zmin, zmax),
        camera=camera,
        size=size,
        grid=true,
        legend=:topright
    )

    # 1. Arm pivot and circular path
    scatter!(p3d, [0.0], [0.0], [0.0], color=:black, markersize=5, label="Pivot (0,0,0)")
    plot!(p3d, X_arm, Y_arm, Z_arm, color=:gray30, lw=2, linestyle=:dash, label="Arm path")

    # 2. Wall Projections
    plot!(p3d, X_cg, Y_cg, fill(zmin, N_pts), color=:gray60, lw=1.5, linestyle=:dot, label="XY-floor proj")
    plot!(p3d, X_cg, fill(ymax, N_pts), Z_cg, color=:gray60, lw=1.5, linestyle=:dashdot, label="XZ-back proj")
    plot!(p3d, fill(xmin, N_pts), Y_cg, Z_cg, color=:gray60, lw=1.5, linestyle=:dash, label="YZ-side proj")

    # 3. Sample Tethers
    if n_tethers > 0
        tether_indices = round.(Int, range(1, N_pts, length=n_tethers + 1))[1:n_tethers]
        for (i, idx) in enumerate(tether_indices)
            plot!(p3d, [X_arm[idx], X_cg[idx]], [Y_arm[idx], Y_cg[idx]], [Z_arm[idx], Z_cg[idx]],
                color=:orange, lw=1.2, alpha=0.5, label=(i == 1 ? "Tether sample" : false))
        end
    end

    # 4. Main Trajectory
    plot!(p3d, X_cg, Y_cg, Z_cg, color=:crimson, lw=3.5, label="Kite CG trajectory")

    # 5. Directional Arrows
    if n_arrows > 0
        arrow_indices = round.(Int, range(10, N_pts - 10, length=n_arrows))
        for idx in arrow_indices
            pt = pos_eval[idx]
            pt_next = pos_eval[idx+2]
            plot!(p3d, [pt[1], pt_next[1]], [pt[2], pt_next[2]], [pt[3], pt_next[3]],
                arrow=arrow(:closed, :head, 0.4, 0.4), color=:darkblue, lw=2.0, label=false)
        end
    end

    if filename !== nothing
        savefig(p3d, filename)
        println("Plot saved to $filename")
    end

    return p3d
end

"""
    animate_kite_3d(sol; fps=30, n_frames=nothing, camera=(40, 25), size=(850, 700), filename="kite_trajectory_animation.gif")

Renders and saves a 3D animated GIF of the kite trajectory.

# Arguments
- `sol`: Solution object or `NamedTuple` containing `time_grid` and `state`.

# Keyword Arguments
- `fps::Int = 30`: Playback frame rate (frames per second). Defaults to `30`.
- `n_frames::Union{Int, Nothing} = nothing`: Total animation frames. If `nothing` (default), 
  calculates `round(Int, tf * fps)` to guarantee 1:1 real-time playback matching trajectory duration.
- `camera = (40, 25)`: Camera viewing angles (azimuth, elevation).
- `size = (850, 700)`: Output figure dimensions in pixels.
- `filename::String = "kite_trajectory_animation.gif"`: Filepath for the generated GIF.
"""
function animate_kite_3d(sol;
    fps::Int=30,
    n_frames::Union{Int,Nothing}=nothing,
    camera=(40, 25),
    size=(850, 700),
    filename::String="kite_trajectory_animation.gif")

    tf_sol = time_grid(sol)[end]

    # Guarantee 1:1 real-time playback speed by default
    total_frames = n_frames === nothing ? max(2, round(Int, tf_sol * fps)) : n_frames

    N_static = 300
    t_static = range(0, tf_sol, length=N_static)

    q_static = [state(sol)(t)[1:4] for t in t_static]
    pos_static = [pos_CG(q) for q in q_static]
    arm_static = [SVector(pars.larm * cos(q[1]), pars.larm * sin(q[1]), 0.0) for q in q_static]

    X_cg = [p[1] for p in pos_static]
    Y_cg = [p[2] for p in pos_static]
    Z_cg = [p[3] for p in pos_static]

    X_arm = [p[1] for p in arm_static]
    Y_arm = [p[2] for p in arm_static]
    Z_arm = [p[3] for p in arm_static]

    pad = 2.0
    xmin, xmax = minimum([X_cg; X_arm; 0.0]) - pad, maximum([X_cg; X_arm; 0.0]) + pad
    ymin, ymax = minimum([Y_cg; Y_arm; 0.0]) - pad, maximum([Y_cg; Y_arm; 0.0]) + pad
    zmin, zmax = 0.0, maximum([Z_cg; 0.0]) + pad

    t_anim = range(0, tf_sol, length=total_frames)

    anim = @animate for t in t_anim
        q_curr = state(sol)(t)[1:4]
        p_cg = pos_CG(q_curr)
        p_arm = SVector(pars.larm * cos(q_curr[1]), pars.larm * sin(q_curr[1]), 0.0)

        # Kite body axes
        Ihat, Jhat, _ = body_axes(q_curr)
        wing_left = p_cg - 1.25 * Jhat
        wing_right = p_cg + 1.25 * Jhat
        chord_tail = p_cg - 0.40 * Ihat
        chord_nose = p_cg + 0.40 * Ihat

        p_frame = plot(
            title="Kite Cycle (t = $(round(t, digits=2)) s / $(round(tf_sol, digits=2)) s)",
            xlabel="X [m]", ylabel="Y [m]", zlabel="Z [m]",
            xlims=(xmin, xmax), ylims=(ymin, ymax), zlims=(zmin, zmax),
            camera=camera,
            size=size,
            legend=false
        )

        # Static background curves & shadows
        plot!(p_frame, X_cg, Y_cg, Z_cg, color=:crimson, lw=1.5, alpha=0.4)
        plot!(p_frame, X_cg, Y_cg, fill(zmin, N_static), color=:gray75, lw=1.0)
        plot!(p_frame, X_cg, fill(ymax, N_static), Z_cg, color=:gray75, lw=1.0)
        plot!(p_frame, fill(xmin, N_static), Y_cg, Z_cg, color=:gray75, lw=1.0)
        plot!(p_frame, X_arm, Y_arm, Z_arm, color=:gray50, lw=1.2, linestyle=:dash)

        # Generator arm & pivot
        scatter!(p_frame, [0.0], [0.0], [0.0], color=:black, markersize=5)
        plot!(p_frame, [0.0, p_arm[1]], [0.0, p_arm[2]], [0.0, p_arm[3]], color=:black, lw=4)
        scatter!(p_frame, [p_arm[1]], [p_arm[2]], [p_arm[3]], color=:black, markersize=4)

        # Dynamic Tether
        plot!(p_frame, [p_arm[1], p_cg[1]], [p_arm[2], p_cg[2]], [p_arm[3], p_cg[3]], color=:orange, lw=2.5)

        # Dynamic Kite Body
        plot!(p_frame, [wing_left[1], wing_right[1]], [wing_left[2], wing_right[2]], [wing_left[3], wing_right[3]], color=:darkblue, lw=4.5)
        plot!(p_frame, [chord_tail[1], chord_nose[1]], [chord_tail[2], chord_nose[2]], [chord_tail[3], chord_nose[3]], color=:red, lw=3.0)
        scatter!(p_frame, [p_cg[1]], [p_cg[2]], [p_cg[3]], color=:darkred, markersize=4)

        # Dynamic shadow dots on bounding planes
        scatter!(p_frame, [p_cg[1]], [p_cg[2]], [zmin], color=:gray40, markersize=4)
        scatter!(p_frame, [p_cg[1]], [ymax], [p_cg[3]], color=:gray40, markersize=4)
        scatter!(p_frame, [xmin], [p_cg[2]], [p_cg[3]], color=:gray40, markersize=4)
    end

    gif(anim, filename, fps=fps)
    println("Animation saved to $filename (Duration: $(round(total_frames / fps, digits=2)) s at $fps FPS)")
    return anim
end

## Half period
# Half-period bounds (half of full period [0.5, 20.0])
TF_HALF_MIN, TF_HALF_MAX = 0.25, 10.0

ocp_half = @def begin
    p = (tf, p_) ∈ R², variable
    t ∈ [0, tf], time
    x = (α, θ, φ, β, α̇, θ̇, φ̇, β̇) ∈ R⁸, state
    u = (uₗ, uᵣ) ∈ R², control

    p_ == 0

    tf >= TF_HALF_MIN
    tf <= TF_HALF_MAX

    -30 * pi / 180 ≤ uₗ(t) ≤ 30 * pi / 180
    -30 * pi / 180 ≤ uᵣ(t) ≤ 30 * pi / 180

    # Phase pinning at t = 0 (center crossing going to the left lobe)
    α(0) == 0
    α̇(0) >= 0
    θ̇(0) <= 0
    φ̇(0) <= 0

    α(t) >= 0
    # θ(t) <= π/2  # Kite above ground
    φ(t) >= 0  # Only y > 0 lobe

    # --- Half-period symmetry boundary condition: x(tf) == S(x(0)) ---
    α(tf) + α(0) == 0
    θ(tf) - θ(0) == 0
    φ(tf) + φ(0) == 0
    β(tf) + β(0) == 0
    α̇(tf) + α̇(0) == 0
    θ̇(tf) - θ̇(0) == 0
    φ̇(tf) + φ̇(0) == 0
    β̇(tf) + β̇(0) == 0

    ẋ(t) == f(x(t), u(t))

    ∫(-generated_power(x(t)) + 1e-3 * sum(abs2, u(t))) → min
end

# Initial guess over [0, T/2] extracted from the limit cycle
init_half = @init ocp_half begin
    tf := lc.t[end] / 2.0

    α(t) := fα(t)
    θ(t) := fθ(t)
    φ(t) := fφ(t)
    β(t) := fβ(t)
    α̇(t) := ForwardDiff.derivative(fα, t)
    θ̇(t) := ForwardDiff.derivative(fθ, t)
    φ̇(t) := ForwardDiff.derivative(fφ, t)
    β̇(t) := ForwardDiff.derivative(fβ, t)
end
# Reflection operators for half-period symmetry
const S_STATE = SA[-1.0, 1.0, -1.0, -1.0, -1.0, 1.0, -1.0, -1.0]
sym_state(x) = S_STATE .* x
sym_control(u) = SVector(u[2], u[1])

function reconstruct_full_solution(sol_half)
    th = time_grid(sol_half)[end]
    T = 2.0 * th

    x_half = state(sol_half)
    u_half = control(sol_half)

    # Periodic state & control closures across [0, 2*th]
    x_full(t) =
        let tm = mod(t, T)
            tm <= th ? x_half(tm) : sym_state(x_half(tm - th))
        end

    u_full(t) =
        let tm = mod(t, T)
            tm <= th ? u_half(tm) : sym_control(u_half(tm - th))
        end

    t_half = time_grid(sol_half)
    t_full = [t_half; t_half[2:end] .+ th]

    return (;
        time_grid=t_full,
        state=x_full,
        control=u_full,
        tf=T,
        objective=2.0 * objective(sol_half)
    )
end

# Make time_grid, state, control, objective work uniformly on the NamedTuple
OptimalControl.time_grid(s::NamedTuple) = s.time_grid
OptimalControl.state(s::NamedTuple) = s.state
OptimalControl.control(s::NamedTuple) = s.control
OptimalControl.objective(s::NamedTuple) = s.objective

# =========================================================
# Run Half-Period OCP with Grid Continuation
# =========================================================
println("--- Solving Half-Period OCP ---")

# 1. Run grid homotopy (8 -> 20 -> 50) with automatic warm-starting and JLD2 caching
sol_init_half = solve(ocp_half, init=init, maxiter=0, grid_size=50)
@time sol_half = run_grid_homotopy(
    ocp_half,
    (8, 20, 50);
    init=init_half,
    cache=:no,                     # :no | :exact | :latest
    cache_dir="applications/controlled/solutions",
    prefix="kite_half_lemniscate",
    solve_options...
)

plot(sol_init_half)
plot(sol_half)

# 2. Reconstruct the full symmetric figure-8 orbit
sol_full = reconstruct_full_solution(sol_half)
tf_full = sol_full.tf

println("\n=== Solution Summary ===")
println("Full period  (T)    = $(round(tf_full, digits=3)) s")
println("Total Power Output  = $(round(-sol_full.objective, digits=2)) W")

# Verify boundary periodic continuity
x0 = state(sol_full)(0.0)
xT = state(sol_full)(tf_full)
println("Periodic loop closure error ||x(0) - x(T)|| = ", norm(x0 - xT))

plot_kite_3d(sol_full)
animate_kite_3d(sol_full)

ts = range(0, tf_full, length=300)

# Evaluate states and body axes across the full figure-8
q_eval = [state(sol_full)(t)[1:4] for t in ts]
pos_eval = [pos_CG(q) for q in q_eval]
