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

tf_init = 3
init = @init ocp begin
    tf := tf_init
    # genuinely periodic guess: α oscillates (drives the generator), others steady
    A = 0.5
    ω = 2π / tf_init
    x(t) := SA[
        A*sin(ω*t),
        pi/6,
        0.0,
        0.0,
        A*ω*cos(ω*t),
        0.0,
        0.0,
        0.0
    ]
    u(t) := SVector(0.0, 0.0)
end

solve_options = (
    backend=:default,
    print_level=5,
    linear_solver="mumps",
    tol=1e-4,
    maxiter=50
    # nlp_scaling_method="gradient-based"

)
base_filename="applications/controlled/solutions/solution_$(hash(ocp))_"


USE_JLD2 = true
filename = base_filename * "8"
sol_ocp = try
    @assert USE_JLD2
    import_ocp_solution(ocp; filename)
catch e
    s = solve(ocp; init=init, grid_size=8,
        solve_options...,
        tol=1e-3
    )
    export_ocp_solution(s; filename)
    s
end

println("OCP solve done (grid 8): objective = ", objective(sol_ocp))
plot(sol_ocp)

# --- grid continuation: warm-start the fine grid from the coarse solution ---
for n in (20, 50)
    try
        filename = base_filename * string(n)
        sol_next = try
            @assert USE_JLD2
            import_ocp_solution(ocp; filename)
        catch e
            s = solve(ocp; init=sol_ocp,
                grid_size=n,
                solve_options...)
            export_ocp_solution(s; filename)
            s
        end
        global sol_ocp = sol_next
        println("OCP solve done (grid $n): objective = ", objective(sol_ocp))
    catch e
        println("grid $n failed: ", sprint(showerror, e)[1:min(end, 200)])
        break
    end
end

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

# Initialize with PM2 trajectory
# Add phase-fixing: OK
# Add state constraints : going right then left or the inverse

# GPU (ExaModel+MadNLP): Metal not supported
# AppleAccelerate: seems marginally slower

# Transform state from 2 <-> 4
# Take parameters p2 and p4
#
plot(sol_ocp)
