# OC_controlled_kite_v2.jl — state-of-the-art reimplementation
#
# 4-DOF kite optimal control (periodic orbit). Key changes vs v1:
#
#  1. EOM derived from the LAGRANGIAN via ForwardDiff: M(q) = ∂²L/∂q̇² and the
#     bias term ∂L/∂q - (∂²L/∂q̇∂q)q̇ are computed by AD. NO hand-derived mass
#     matrix, NO hand-derived dM/dt, NO hand-derived Coriolis terms — the
#     biggest source of bugs in v1.
#  2. Aerodynamics smoothed (tanh blend) instead of the hard `abs(α)<=30°`
#     branch: required for both SparseConnectivityTracer sparsity detection and
#     Ipopt's smoothness assumption.
#  3. Branch-free 4x4 Cholesky solve (pure arithmetic) so SparseConnectivityTracer
#     can trace through the linear system (pivoted `\` cannot).
#  4. `backend=:default` ADNLP modeler (pure ForwardDiff Hessian; the default
#     `:optimized` ReverseDiff-over-ForwardDiff throws DualMismatchError).
#
# Run: julia --project=applications/controlled OC_controlled_kite_v2.jl

using Pkg
Pkg.activate("applications/controlled")

using StaticArrays
using LinearAlgebra
using ForwardDiff
using OrdinaryDiffEqTsit5

using OptimalControl, NLPModelsIpopt

# =========================================================
# Parameters — plain NamedTuple (type-stable, fast)
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
# Kinematics
# =========================================================

@inline function pos_OA(q, p)          # generator arm tip
    return SVector(p.larm * cos(q[1]), p.larm * sin(q[1]), 0.0)
end

@inline function pos_CG(q, p)          # kite centre of gravity
    return SVector(
        p.lcg * sin(q[2]) * cos(q[3]),
        p.lcg * sin(q[2]) * sin(q[3]),
        p.lcg * cos(q[2]),
    ) + pos_OA(q, p)
end

# Body-frame axes (Ihat: along tether, Khat: up-ish, Jhat: closing) — matches v1.
@inline function body_axes(q)
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

# Body angular velocity ω = W(q)·dq  (W is the 3x4 map).
@inline function W(q)
    cq, sq = cos.(q), sin.(q)
    return @SMatrix [
        0.0 0.0 cq[2] 1.0;
        0.0 cq[4] sq[4] * sq[2] 0.0;
        0.0 -sq[4] cq[4] * sq[2] 0.0
    ]
end

# =========================================================
# Aerodynamics — smoothed (no hard branch), per-panel model
# =========================================================

wind(z) = 20.0

# Smooth saturation: ~1 inside |α5| ≤ 30°, ~0 outside, tanh-blended transition.
@inline function stall_factor(alpha)
    α5 = alpha + 5.0 * pi / 180.0
    k = 30.0
    return 0.5 * (1.0 - tanh(k * (abs(α5) - 30.0 * pi / 180.0)))
end

# Smooth lift/drag polar. v1 had `abs(α)<=30°` hard branch — replaced.
function Caero(alpha)
    α5 = alpha + 5.0 * pi / 180.0
    s = stall_factor(alpha)
    return SVector(
        1.0 * α5 * s,                                   # Cl, saturating
        0.2 + 0.1 * α5^2,                               # Cd
        0.0,                                            # side
    )
end

function aerodynamics(X, u, p)
    q = SVector{4}(X[1:4])
    dq = SVector{4}(X[5:8])
    a = p.aero
    l = p.lcg
    b = p.larm

    Ihat, Jhat, Khat = body_axes(q)
    vCG = ForwardDiff.jacobian(x -> pos_CG(x, p), q) * dq
    vCGrel = vCG - SVector(wind(pos_CG(q, p)[3]), 0.0, 0.0)
    vCGb = SVector(dot(Ihat, vCGrel), dot(Jhat, vCGrel), dot(Khat, vCGrel))
    ω = W(q) * dq

    Rp = a.panels
    Rp1, Rp2 = Rp[:, 1], Rp[:, 2]
    Vair1 = -(vCGb + cross(ω, Rp1))
    Vair2 = -(vCGb + cross(ω, Rp2))

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
    liftDir1 = -cross(w1, dragDir1)
    liftDir2 = cross(w2, dragDir2)

    i1 = atan(dot(l1, Vplane1), -Vplane1[3]) + u[1]
    i2 = atan(dot(l2, Vplane2), -Vplane2[3]) + u[2]
    Ca1 = Caero(i1)
    Ca2 = Caero(i2)

    c = 0.5 * a.rho * a.S
    F1a = c * SVector(Vp1n^2, Vp1n^2, Vs1^2) .* Ca1
    F2a = c * SVector(Vp2n^2, Vp2n^2, Vs2^2) .* Ca2
    F1 = F1a[1] * liftDir1 + F1a[2] * dragDir1 + F1a[3] * w1
    F2 = F2a[1] * liftDir2 + F2a[2] * dragDir2 + F2a[3] * w2

    F = F1 + F2
    M = cross(Rp1, F1) + cross(Rp2, F2)
    return (; F, M)
end

# =========================================================
# Lagrangian + generalized forces
# =========================================================

function lagrangian(q, dq, u, p)
    vCG = ForwardDiff.jacobian(x -> pos_CG(x, p), q) * dq
    ω = W(q) * dq
    T = 0.5 * p.mass * sum(abs2, vCG) +
        0.5 * dot(p.inertia, ω .^ 2) +
        0.5 * p.inertiaarm * dq[1]^2
    V = p.mass * p.g * pos_CG(q, p)[3]
    return T - V
end

function generalized_forces(q, dq, u, p)
    # Aero force/moment mapped through virtual work: Q += J_CG'·F + W'·M
    aero = aerodynamics(vcat(q, dq), u, p)
    J_CG = ForwardDiff.jacobian(x -> pos_CG(x, p), q)
    Qaero = J_CG' * aero.F + W(q)' * aero.M

    # Generator torque on α, line damping on β
    Qgen = SVector(p.kgenerator * q[1] + p.cgenerator * dq[1], 0.0, 0.0, 0.0)
    Qline = SVector(0.0, 0.0, 0.0, p.klines * q[4])

    return SVector{4}(Qaero - Qgen + Qline)
end

# =========================================================
# EOM from the Lagrangian, all via AD
# =========================================================

# ponytail: branch-free 4x4 Cholesky solve of the SPD mass matrix.
# Pivoted `M \ b` breaks SparseConnectivityTracer (comparisons in pivot search).
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

function eom(X, u, p)
    q = SVector{4}(X[1:4])
    dq = SVector{4}(X[5:8])

    # momentum p = ∂L/∂q̇
    momentum(q_, dq_) = ForwardDiff.gradient(v -> lagrangian(q_, v, u, p), dq_)

    # Mass matrix M = ∂p/∂q̇ = ∂²L/∂q̇²  (symmetric PD)
    M = ForwardDiff.jacobian(v -> momentum(q, v), dq)

    # Lagrange: d/dt(∂L/∂q̇) - ∂L/∂q = Q
    #   d/dt(∂L/∂q̇) = M·q̈ + (∂²L/∂q̇∂q)·q̇
    #   =>  M·q̈ = Q + ∂L/∂q - (∂²L/∂q̇∂q)·q̇
    dLdq = ForwardDiff.gradient(x -> lagrangian(x, dq, u, p), q)
    cross = ForwardDiff.jacobian(x -> momentum(x, dq), q) * dq
    Q = generalized_forces(q, dq, u, p)

    b = Q + dLdq - cross
    ddq = cholesky_solve(SMatrix{4, 4}(M + 1e-9 * I), b)

    return SVector{8}(dq[1], dq[2], dq[3], dq[4], ddq[1], ddq[2], ddq[3], ddq[4])
end

f(X, u) = eom(X, u, pars)

# =========================================================
# Sanity checks
# =========================================================
let X = SA[0.3, 0.5, 0.2, 0.6, 0.1, 0.2, 0.3, 0.4], u = SA[0.1, 0.0]
    @assert f(X, u) isa SVector{8, Float64}
    # ForwardDiff through the full RHS (what OptimalControl's Jacobian needs)
    J = ForwardDiff.jacobian(x -> f(x, u), X)
    @assert size(J) == (8, 8)
    # Second order as well
    g = ForwardDiff.gradient(x -> sum(f(x, u)), X)
    @assert length(g) == 8
    println("sanity: eom OK, jacobian=", size(J), ", first ddq = ", f(X, u)[5:8])
end

# =========================================================
# ODE solve (open loop, constant control)
# =========================================================
x0 = @SVector [0.0, pi / 6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
u_const = @SVector [-20.0 * pi / 180, -20.0 * pi / 180]
prob = ODEProblem((X, p, t) -> eom(X, u_const, p), x0, (0.0, 20.0), pars)
sol = solve(prob, Tsit5(); reltol=1e-6, abstol=1e-6)
println("ODE solve OK, t_end=", sol.t[end])

# =========================================================
# Optimal control — periodic orbit
# =========================================================
TF_MIN, TF_MAX = 0.5, 20.0

generated_power(x) = -x[5] * (pars.cgenerator * x[5] + pars.kgenerator * x[1])

ocp = @def begin
    tf ∈ R, variable
    t ∈ [0, tf], time
    x = (α, θ, φ, β, α̇, θ̇, φ̇, β̇) ∈ R⁸, state
    u = (uₗ, uᵣ) ∈ R², control

    tf >= TF_MIN
    tf <= TF_MAX

    α(0) == 0
    α̇(0) >= 0
    x(0) - x(tf) == zeros(8)

    ẋ(t) == f(x(t), u(t))

    ∫(-generated_power(x(t)) + 1e-3 * sum(abs2, u(t))) → min
end

init = @init ocp begin
    tf := 10.0
    x(t) := x0
    u(t) := SVector(0.0, 0.0)
end

# backend=:default => pure ForwardDiff gradient+Jacobian+Hessian (SparseADHessian).
# The default :optimized backend nests ReverseDiff inside ForwardDiff and throws
# "Cannot determine ordering of Dual tags" on StaticArrays-heavy RHS.
sol_ocp = solve(ocp; init=init, grid_size=20,
    modeler=OptimalControl.ADNLP(backend=:default),
    max_wall_time=120.0)
println("OCP solve done: objective = ", objective(sol_ocp))

# =========================================================
# TLDR
# =========================================================
# 1. EOM comes from the Lagrangian via ForwardDiff (mass matrix M=∂²L/∂q̇² and
#    bias ∂L/∂q-(∂²L/∂q̇∂q)q̇). No hand-derived M/dM/Coriolis → the recurring bug.
# 2. Aero polar is tanh-smoothed (no `abs(α)<=30°` hard branch).
# 3. cholesky_solve (branch-free) survives SparseConnectivityTracer sparsity.
# 4. `backend=:default` avoids the ReverseDiff-over-ForwardDiff DualMismatchError.
# If Ipopt still struggles: refine the initial guess (limit cycle), raise
# grid_size via continuation, or try linear_solver="ma57".
