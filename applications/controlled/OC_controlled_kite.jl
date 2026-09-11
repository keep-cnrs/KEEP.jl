using Pkg
Pkg.activate("applications/controlled")

using ComponentArrays: ComponentArray as CA
using StaticArrays
using LinearAlgebra
using OrdinaryDiffEqTsit5
using Plots

using OptimalControl, NLPModelsIpopt
using ForwardDiff

#=
The Coriolis term computed the mass matrix temporal derivative dM using positions q instead of velocities dq. This has been corrected to use dq.

The rotational kinetic energy term for the arm in the lagrangian function forgot to square the angular velocity dq(1). This has been fixed to dq[1]^2.
=#

# =========================================================
# Dynamics and Kinematics Functions
# =========================================================

function my_norm(x)
    return sqrt(sum(abs2, x))
end

@inline q_dq(X) = (X[SA[1, 2, 3, 4]], X[SA[5, 6, 7, 8]])

function aerodynamics(X, u, pars)
    q, dq = q_dq(X)

    b = pars.larm              # Length generator arm [m]
    l = pars.lcg               # distance bottom generator arm to CG ~ Length lines [m]
    S = pars.aero.S / 2.0      # Total wing surface [m^2]
    rho = pars.aero.rho          # Atmospheric density [kg / m^3]
    Rp = SMatrix{3,2}(pars.aero.panels) # Relative position
    delta = pars.aero.delta        # Dihedral angle [rad]

    # Rotation matrix components and angular velocity
    sq = sin.(q)
    cq = cos.(q)

    Ihat = SVector(sq[2] * cq[3], sq[2] * sq[3], cq[2])
    Jhat = SVector(
        -cq[4] * sq[3] - sq[4] * cq[2] * cq[3],
        cq[4] * cq[3] - sq[4] * cq[2] * sq[3],
        sq[4] * sq[2],
    )
    Khat = SVector(
        sq[4] * sq[3] - cq[4] * cq[2] * cq[3],
        -sq[4] * cq[3] - cq[4] * cq[2] * sq[3],
        cq[4] * sq[2],
    )

    W = SVector(
        dq[4] + cq[2] * dq[3],
        sq[4] * sq[2] * dq[3] + cq[4] * dq[2],
        cq[4] * sq[2] * dq[3] - sq[4] * dq[2],
    )

    # Position and velocity of the panels
    zCG = l * cq[2]
    vCG = SVector(
        -b * sq[1] * dq[1] + l * cq[2] * cq[3] * dq[2] - l * sq[2] * sq[3] * dq[3],
        b * cq[1] * dq[1] + l * cq[2] * sq[3] * dq[2] + l * sq[2] * cq[3] * dq[3],
        -l * sq[2] * dq[2],
    )
    vCGrel = vCG - SVector(wind(zCG), 0.0, 0.0)

    vCGb = SVector(dot(Ihat, vCGrel), dot(Jhat, vCGrel), dot(Khat, vCGrel))

    Rp1 = Rp[:, 1]
    Rp2 = Rp[:, 2]

    Vair1 = -(vCGb + cross(W, Rp1))
    Vair2 = -(vCGb + cross(W, Rp2))

    # Removing the side component of the velocity
    cdelta = cos(delta)
    sdelta = sin(delta)

    w1 = SVector(-sdelta, cdelta, 0.0)
    w2 = SVector(-sdelta, -cdelta, 0.0)
    l1 = SVector(cdelta, sdelta, 0.0)
    l2 = SVector(cdelta, -sdelta, 0.0)

    Vside1n = dot(w1, Vair1)
    Vside2n = dot(w2, Vair2)

    Vplane1 = Vair1 - w1 * Vside1n
    Vplane2 = Vair2 - w2 * Vside2n

    Vplane1n = my_norm(Vplane1)
    Vplane2n = my_norm(Vplane1)

    dragDir1 = Vplane1 / Vplane1n
    dragDir2 = Vplane2 / Vplane2n

    liftDir1 = -cross(w1, dragDir1)
    liftDir2 = cross(w2, dragDir2)

    # Computing aerodynamics coefficients as a function of the incidence of each panel
    i1 = atan(dot(l1, Vplane1), -Vplane1[3]) + u[1]
    i2 = atan(dot(l2, Vplane2), -Vplane2[3]) + u[2]

    Ca1 = Caero(i1)
    Ca2 = Caero(i2)

    F1a = 0.5 * rho * S * SVector(Vplane1n^2, Vplane1n^2, Vside1n^2) .* Ca1
    F2a = 0.5 * rho * S * SVector(Vplane2n^2, Vplane2n^2, Vside2n^2) .* Ca2

    F1 = F1a[1] * liftDir1 + F1a[2] * dragDir1 + F1a[3] * w1
    F2 = F2a[1] * liftDir2 + F2a[2] * dragDir2 + F2a[3] * w2

    # Total force and momentum
    F = F1 + F2
    M = cross(Rp1, F1) + cross(Rp2, F2)

    return (; F, M, i1, i2)
end

function generator(alpha, dalpha, pars)
    k = pars.kgenerator
    c = pars.cgenerator
    return -k * alpha - c * dalpha
end

function wind(z)
    TU = 1.0
    return 20.0 * TU
end

function Caero(alpha)
    return SVector(
        1.0 * (alpha + 5.0 * pi / 180.0) * (abs(alpha) <= 30.0 * pi / 180.0),
        0.2 + 0.1 * (alpha + 5.0 * pi / 180.0)^2,
        0.0,
    )
end

function lagrangian(X, pars)
    q, dq = q_dq(X)

    b = pars.larm
    l = pars.lcg
    m = pars.mass
    A, B, C = pars.inertia[1], pars.inertia[2], pars.inertia[3]
    Iarm = pars.inertiaarm
    g = pars.g

    cq = cos.(q)
    sq = sin.(q)

    rCG = SVector(b * cq[1] + l * sq[2] * cq[3], b * sq[1] + l * sq[2] * sq[3], l * cq[2])
    vCG = SVector(
        -b * sq[1] * dq[1] + l * cq[2] * cq[3] * dq[2] - l * sq[2] * sq[3] * dq[3],
        b * cq[1] * dq[1] + l * cq[2] * sq[3] * dq[2] + l * sq[2] * cq[3] * dq[3],
        -l * sq[2] * dq[2],
    )

    Ihat = SVector(sq[2] * cq[3], sq[2] * sq[3], cq[2])
    Jhat = SVector(
        -cq[4] * sq[3] - sq[4] * cq[2] * cq[3],
        cq[4] * cq[3] - sq[4] * cq[2] * sq[3],
        sq[4] * sq[2],
    )
    Khat = SVector(
        sq[4] * sq[3] - cq[4] * cq[2] * cq[3],
        -sq[4] * cq[3] - cq[4] * cq[2] * sq[3],
        cq[4] * sq[2],
    )

    w = SVector(
        dq[4] + cq[2] * dq[3],
        sq[4] * sq[2] * dq[3] + cq[4] * dq[2],
        cq[4] * sq[2] * dq[3] - sq[4] * dq[2],
    )

    # NOTE: The original MATLAB code missed squaring the angular velocity (dq[1]) 
    # of the arm in the rotational kinetic energy. We've corrected this to dq[1]^2.
    T_kin =
        0.5 * m * sum(abs2, vCG) +
        0.5 * dot(SVector(A, B, C), w .^ 2) +
        0.5 * Iarm * dq[1]^2
    V_pot = m * g * rCG[3]

    L = T_kin - V_pot

    return (; L, T=T_kin, V=V_pot, rCG, vCG, Ihat, Jhat, Khat, w)
end

# ponytail: branch-free 4x4 Cholesky solve of the SPD mass-matrix system.
# The pivoted `M \ b` breaks OptimalControl's sparsity detection (SparseConnectivityTracer
# can't evaluate the `if absi > amax` pivot comparison). Pure arithmetic here keeps the tracer happy.
# Alternative for benchmarking: ddq = -inv(M) * (dM * dq - dLdq - Q)  (StaticArrays 4x4 inv is also branch-free)
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
    return SA[x1, x2, x3, x4]
end

function eom(X, u, pars, t=0)
    q, dq = q_dq(X)

    b = pars.larm
    l = pars.lcg
    m = pars.mass
    A, B, C = pars.inertia[1], pars.inertia[2], pars.inertia[3]
    Iarm = pars.inertiaarm
    g = pars.g
    klines = pars.klines

    cq = cos.(q)
    sq = sin.(q)
    sq31 = sin(q[3] - q[1])
    cq31 = cos(q[3] - q[1])

    Mt = @SMatrix [
        (m*b^2+Iarm)/2 0 0 0;
        b*l*cq[2]*sq31 (m*l^2+B*cq[4]^2+C*sq[4]^2)/2 0 0;
        b*l*sq[2]*cq31 0 (A*cq[2]^2+sq[2]^2*(m*l^2+B*sq[4]^2+C*cq[4]^2))/2 0;
        0 0 A*cq[2] A/2
    ]
    M = Mt + Mt'

    dMt1 = @SMatrix [
        0 0 0 0;
        -b*l*cq[2]*sq31 0 0 0;
        b*l*sq[2]*sq31 0 0 0;
        0 0 0 0
    ]
    dM1 = dMt1 + dMt1'

    dMt2 = @SMatrix [
        0 0 0 0;
        -b*l*sq[2]*sq31 0 0 0;
        b*l*cq[2]*cq31 0 cq[2]*sq[2]*(-A+m*l^2+B*sq[4]^2+C*cq[4]^2) 0;
        0 0 -A*sq[2] 0
    ]
    dM2 = dMt2 + dMt2'

    dMt3 = @SMatrix [
        0 0 0 0;
        b*l*cq[2]*cq31 0 0 0;
        -b*l*sq[2]*sq31 0 0 0;
        0 0 0 0
    ]
    dM3 = dMt3 + dMt3'

    dMt4 = @SMatrix [
        0 0 0 0;
        0 cq[4]*sq[4]*(C-B) 0 0;
        0 0 -sq[2]^2*cq[4]*sq[4]*(C-B) 0;
        0 0 0 0
    ]
    dM4 = dMt4 + dMt4'

    # NOTE: The original MATLAB code erroneously used positions q(1)*dM1 + q(2)*dM2...
    # The time derivative of the mass matrix should be constructed using velocities (dq).
    # Corrected this to dq[1]*dM1 + ... to ensure physically accurate mass-matrix derivation.
    dM = dq[1] * dM1 + dq[2] * dM2 + dq[3] * dM3 + dq[4] * dM4

    Ihat = SVector(sq[2] * cq[3], sq[2] * sq[3], cq[2])
    Jhat = SVector(
        -cq[4] * sq[3] - sq[4] * cq[2] * cq[3],
        cq[4] * cq[3] - sq[4] * cq[2] * sq[3],
        sq[4] * sq[2],
    )
    Khat = SVector(
        sq[4] * sq[3] - cq[4] * cq[2] * cq[3],
        -sq[4] * cq[3] - cq[4] * cq[2] * sq[3],
        cq[4] * sq[2],
    )

    Rref2body = vcat(Ihat', Jhat', Khat')

    deltaRcg =
        Rref2body * @SMatrix [
            -b*sq[1] l*cq[2]*cq[3] -l*sq[2]*sq[3] 0.0;
            b*cq[1] l*cq[2]*sq[3] l*sq[2]*cq[3] 0.0;
            0.0 -l*sq[2] 0.0 0.0
        ]

    deltaW = @SMatrix [
        0.0 0.0 cq[2] 1.0;
        0.0 cq[4] sq[4]*sq[2] 0.0;
        0.0 -sq[4] cq[4]*sq[2] 0.0
    ]

    aero = aerodynamics(X, u, pars)
    Fa = aero.F
    Ma = aero.M

    Qaero = deltaRcg' * Fa + deltaW' * Ma

    Qgenerator = SVector(generator(q[1], dq[1], pars), 0.0, 0.0, 0.0)
    Qlines = SVector(0.0, 0.0, 0.0, -klines * q[4])

    Q = Qaero + Qgenerator + Qlines

    dV = SVector(0.0, -m * g * l * sq[2], 0.0, 0.0)
    dLdq =
        SVector(
            dot(dq, dM1 * dq), dot(dq, dM2 * dq), dot(dq, dM3 * dq), dot(dq, dM4 * dq)
        ) / 2.0 - dV

    # Calculate states derivative avoiding explicit inversion operations over allocations
    # @warn "Regularizing mass matrix M in eom" maxlog = 1
    M += 1e-6 * I  # one(eltype(M)) *
    # ddq = -M \ (dM*dq - dLdq - Q)
    # ddq = -inv(M) * (dM * dq - dLdq - Q)   # alternative for benchmarking: StaticArrays 4x4 inv is branch-free
    ddq = -cholesky_solve(M, dM * dq - dLdq - Q)

    return vcat(dq, ddq)
end

# =========================================================
# Simulation Initialization
# =========================================================

LU = 1.0
TU = 1.0

const pars = CA(;
    larm=2.0 / LU,
    lcg=20.0 / LU,
    mass=5.0,
    inertia=5.0 .* [5.0, 0.1, 0.5] ./ LU^2,
    inertiaarm=(2.0 / LU)^2 * 2.0,
    g=9.81 * TU^2 / LU,
    kgenerator=0.0 * TU^2 / LU,
    cgenerator=1.0 * TU / LU,
    klines=0.0 * TU^2 / LU,
    aero=(
        S=5.0 / LU^2,
        rho=1.225 / LU^3,
        panels=[
            0.1/LU 0.1/LU;
            1.25/LU -1.25/LU;
            0.0 0.0
        ],
        delta=30.0 * pi / 180.0,
    ),
)

x0 = @SVector [0.0, pi / 6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
const u_const = @SVector [-20.0 * pi / 180.0, -20.0 * pi / 180.0]
const u_const = @SVector [0, 0]

f_ode = (X, pars, t=0) -> eom(X, u_const, pars, t)
# @benchmark f_ode(x0, pars)
# @profview_allocs @benchmark f_ode(x0, pars)
# @profview @benchmark f_ode(x0, pars)

# Define simulation problem
tspan = (0.0, 20.0 / TU)
prob = ODEProblem(f_ode, x0, tspan, pars)

## Solve
sol = solve(prob, Tsit5(); reltol=1e-6, abstol=1e-6, progress=true)

# =========================================================
# Figures & Visualization
# =========================================================
t_plot = range(0.0, 20.0 / TU; length=101)
x_plot = sol.(t_plot)

# 1. Trajectory Animation (Replaces the drawnow interactive MATLAB loop)
calc_OA(x, pars) = pars.larm * SVector(cos(x[1]), sin(x[1]), 0.0)
function calc_OG(x, pars)
    return calc_OA(x, pars) +
           pars.lcg * SVector(sin(x[2]) * cos(x[3]), sin(x[2]) * sin(x[3]), cos(x[2]))
end

OA = calc_OA.(x_plot, Ref(pars))
OG = calc_OG.(x_plot, Ref(pars))

if false
    anim = @animate for j in 1:length(t_plot)
        plot(
            [0.0, pars.larm],
            [0.0, 0.0],
            [0.0, 0.0],
            color=:black,
            linewidth=2,
            label=false,
            aspect_ratio=:equal,
            showaxis=false,
            grid=false,
        )
        plot!(
            [0.0, 0.0], [0.0, pars.larm], [0.0, 0.0], color=:black, linewidth=2, label=false
        )
        plot!(
            [0.0, 0.0], [0.0, 0.0], [0.0, pars.larm], color=:black, linewidth=2, label=false
        )

        # Arm and Lines
        plot!(
            [0.0, OA[j][1]],
            [0.0, OA[j][2]],
            [0.0, OA[j][3]],
            color=:red,
            linewidth=2,
            label=false,
        )
        plot!(
            [OG[j][1], OA[j][1]],
            [OG[j][2], OA[j][2]],
            [OG[j][3], OA[j][3]],
            color=:blue,
            linewidth=2,
            label=false,
        )

        xlims!(-pars.lcg, pars.lcg)
        ylims!(-pars.lcg, pars.lcg)
        zlims!(0, pars.lcg)
    end
    g = gif(anim; fps=30) # Use this if you wish to dump a fluid video!
    display(g)

    # 2. Beta angle plot

    get_beta(x) = x[2] * (180.0 / pi)
    p_beta = plot(
        t_plot .* TU,
        get_beta.(x_plot);
        color=:black,
        linewidth=2,
        xlabel="Time [s]",
        ylabel="β [deg]",
        framestyle=:box,
        legend=false,
    )
    display(p_beta)

    # 3. Aerodynamic incidence plot
    get_i1(x, u, p) = aerodynamics(x, u, p).i1 * (180.0 / pi)
    get_i2(x, u, p) = aerodynamics(x, u, p).i2 * (180.0 / pi)

    p_inc = plot(
        t_plot .* TU,
        get_i1.(x_plot, Ref(u_const), Ref(pars));
        label="Panel 1",
        linewidth=2,
        xlabel="Time [s]",
        ylabel="Incidence [deg]",
        framestyle=:box,
    )
    plot!(
        p_inc,
        t_plot .* TU,
        get_i2.(x_plot, Ref(u_const), Ref(pars));
        label="Panel 2",
        linewidth=2,
    )
    display(p_inc)
end

#########################
## Generate a trajectory with simplified model
#########################

using KEEP.PointMassPara: build_vbpara
using KEEP.PointMass4: τ_to_θφ
using KEEP.LimitCycle: compute_limit_cycle

vbp = build_vbpara()
const lc = compute_limit_cycle(vbp; sense=(+), save_everystep=true)

fα(t) = lc(t)[1]
fθ(t) = τ_to_θφ(lc(t)[2], vbp)[1]
fφ(t) = τ_to_θφ(lc(t)[2], vbp)[2]

function fpos(t)
    LU = vbp.l
    a, θl, φl = fα(t), fθ(t), fφ(t)
    return LU .* [
        cos(a) + vbp.r * sin(θl) * cos(φl),
        sin(a) + vbp.r * sin(θl) * sin(φl),
        vbp.r * cos(θl),
    ]
end

function fβ(t)
    v = ForwardDiff.derivative(fpos, t)
    θl, φl = fθ(t), fφ(t)
    J0 = [-sin(φl), cos(φl), 0.0] # Side axis at β=0
    K0 = [-cos(θl) * cos(φl), -cos(θl) * sin(φl), sin(θl)] # Up axis at β=0
    return atan(
        -(v[1] * J0[1] + v[2] * J0[2] + v[3] * J0[3]),
        (v[1] * K0[1] + v[2] * K0[2] + v[3] * K0[3]),
    )
end

function fx(t)
    α = fα(t)
    θ = fθ(t)
    φ = fφ(t)
    β = fβ(t)
    dα = ForwardDiff.derivative(fα, t)
    dθ = ForwardDiff.derivative(fθ, t)
    dφ = ForwardDiff.derivative(fφ, t)
    dβ = ForwardDiff.derivative(fβ, t)
    return SA[α, θ, φ, β, dα, dθ, dφ, dβ]
end

"""
Map (α, τ, dα, dτ) to (α, θ, φ, β, dα, dθ, dφ, dβ)
LLM: correct and/or finish"""
function f2to4(x)
    α, τ, dα, dτ, = x
    θ, φ = τ_to_θφ(τ, vbp)
    return dθ, dφ = ForwardDiff.derivative(τ -> τ_to_θφ(τ, vbp), τ) * dτ
end

"""Find closest (α, τ, dα, dτ) with L2 norm
LLM: finish"""
function f4to2(x) end

# =========================================================
# State Mapping Functions (2D Reduced State <-> 4D Full State)
# =========================================================

"""
Compute system position relative to origin given arm angle α and path coordinate τ.
"""
function position_from_α_τ(α, τ, vbp)
    LU = vbp.l
    r = vbp.r
    θ, φ = τ_to_θφ(τ, vbp)
    return LU .*
           SVector(cos(α) + r * sin(θ) * cos(φ), sin(α) + r * sin(θ) * sin(φ), r * cos(θ))
end

"""
Compute roll angle β given 2D reduced state (α, τ, dα, dτ).
"""
function calc_beta(α, τ, dα, dτ, vbp)
    θ, φ = τ_to_θφ(τ, vbp)

    # Velocity vector v = d/dt position_from_α_τ
    v = ForwardDiff.derivative(t -> position_from_α_τ(α + dα * t, τ + dτ * t, vbp), 0.0)

    J0 = SVector(-sin(φ), cos(φ), 0.0)
    K0 = SVector(-cos(θ) * cos(φ), -cos(θ) * sin(φ), sin(θ))

    return atan(-dot(v, J0), dot(v, K0))
end

"""
Map 2D reduced state vector (α, τ, dα, dτ) to full 8-element state vector (α, θ, φ, β, dα, dθ, dφ, dβ).
"""
function f2to4(x)
    α, τ, dα, dτ = x[1], x[2], x[3], x[4]

    # 1. Map τ -> (θ, φ) and dτ -> (dθ, dφ)
    θ, φ = τ_to_θφ(τ, vbp)
    dθ_dτ, dφ_dτ = ForwardDiff.derivative(t -> SVector(τ_to_θφ(t, vbp)...), τ)
    dθ = dθ_dτ * dτ
    dφ = dφ_dτ * dτ

    # 2. Compute β and dβ = dβ/dt
    β = calc_beta(α, τ, dα, dτ, vbp)
    dβ = ForwardDiff.derivative(t -> calc_beta(α + dα * t, τ + dτ * t, dα, dτ, vbp), 0.0)

    return SA[α, θ, φ, β, dα, dθ, dφ, dβ]
end

"""
Find tau that minimizes L2 distance between τ_to_θφ(τ, vbp) and target angles (θ, φ).
Uses grid search followed by Newton-Raphson refinement.
"""
function find_closest_tau(θ, φ, vbp; τ_range=range(-2π, 2π; length=400))
    E(τ_val) = sum(abs2, SVector(τ_to_θφ(τ_val, vbp)...) - SA[θ, φ])

    # Grid search for robust global initial guess
    best_τ = first(τ_range)
    min_E = E(best_τ)
    for τ_val in τ_range
        val = E(τ_val)
        if val < min_E
            min_E = val
            best_τ = τ_val
        end
    end

    # Newton-Raphson quadratic refinement
    τ = best_τ
    for _ in 1:10
        dE = ForwardDiff.derivative(E, τ)
        ddE = ForwardDiff.derivative(t -> ForwardDiff.derivative(E, t), τ)
        if abs(ddE) < 1e-12
            break
        end
        step = dE / ddE
        τ -= step
        if abs(step) < 1e-10
            break
        end
    end
    return τ
end

"""
Map full 8-element state vector (α, θ, φ, β, dα, dθ, dφ, dβ) to 2D reduced state vector (α, τ, dα, dτ)
using L2 distance minimization for position and least-squares projection for velocity.
"""
function f4to2(x)
    α, θ, φ, β, dα, dθ, dφ, dβ = x

    # 1. Find optimal path coordinate τ* minimizing L2 norm ||(θ(τ), φ(τ)) - (θ, φ)||²
    τ = find_closest_tau(θ, φ, vbp)

    # 2. Least-squares projection of (dθ, dφ) onto tangent direction J = d(θ, φ)/dτ
    J = ForwardDiff.derivative(t -> SVector(τ_to_θφ(t, vbp)...), τ)
    dτ = dot(J, SA[dθ, dφ]) / (sum(abs2, J) + 1e-12)

    return SA[α, τ, dα, dτ]
end

#########################
## Reproducing the trajectory
#########################

x0 = fx(0)
prob = ODEProblem(f_ode, x0, tspan, pars)
sol = solve(prob, Tsit5(); reltol=1e-6, abstol=1e-6, progress=true)
plot(sol)

#########################
## Optimal Control
#########################
TF_MIN, TF_MAX = (1/2, 2) .* lc.t[end]

generated_power(x) = -x[5] * generator(x[1], x[5], pars)

f(x, u) = eom(x, u, pars)
f(rand(8), rand(2))

# Suivi de trajectoire
# ocp_traj = @def begin
#     tf ∈ R, variable
#     t ∈ [0, tf], time
#     x = (α, θ, φ, β, α̇, θ̇, φ̇, β̇) ∈ R⁸, state
#     u = (uₗ, uᵣ) ∈ R², control

#     α(0) == 0
#     α̇(0) >= 0
#     x(0) - x(tf) == zeros(8)

#     ẋ(t) == f(x(t), u(t))

#     ∫() → min
# end

# Orbite périodique minimisant le contrôle
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

# Initial guess
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

if false
    init = nothing
end
# default :optimized ADNLP backend uses ReverseDiff-over-ForwardDiff Hessian
# (SparseReverseADHessian), which fails with ForwardDiff "Cannot determine ordering
# of Dual tags" on this StaticArrays-heavy eom. The :default backend uses pure
# ForwardDiff (SparseADHessian) and works.
options = (; init=init, grid_size=20, backend=:default)
solve(ocp; options..., max_wall_time=1.0)
@profview solve(ocp; options..., max_wall_time=10.0)
# @profview solve(ocp; options..., max_walltime = 10)
# @profview_allocs solve(ocp; grid_size=10, options...)
# solve(ocp, init=init)

# =========================================================
# BUG REPORT (2026-08-11)
# =========================================================
# Two independent autodiff bugs were found and fixed in this file:
#
# BUG 1 — `M \ rhs` breaks OptimalControl's sparsity detection.
#   `eom` solved the mass-matrix system with `ddq = -M \ (dM*dq - dLdq - Q)`.
#   OptimalControl's ADNLP modeler computes the Jacobian/Hessian sparsity pattern
#   by pushing SparseConnectivityTracer.GradientTracer numbers through `eom`.
#   StaticArrays' PIVOTED `lu` does `if absi > amax` (lu.jl:169), and `>` is not
#   Bool-safe for GradientTracer -> "TypeError: non-boolean used in boolean context".
#   (ForwardDiff itself was never the problem: its Duals compare like Float64s.)
#   FIX: branch-free 4x4 Cholesky solve `cholesky_solve` (pure + - * / sqrt).
#   `M` is symmetric positive-definite after the `+1e-6*I` regularization, so
#   Cholesky is valid. An `inv(M)` alternative is kept commented for benchmarking
#   (StaticArrays' 4x4 `inv` is also branch-free). Reproduced & asserted in
#   test_ForwardDiff_LinearAlgebra.jl.
#
# BUG 2 — default `:optimized` ADNLP backend fails in the Hessian.
#   The default backend uses SparseReverseADHessian (ReverseDiff-over-ForwardDiff),
#   which throws "Cannot determine ordering of Dual tags" on this StaticArrays-heavy
#   `eom`. FIX: `modeler=OptimalControl.ADNLP(backend=bypass(:default))`, which uses
#   the pure-ForwardDiff SparseADHessian.
#
# After both fixes the pipeline runs end-to-end (sparsity, Jacobian, Hessian all
# compute). Note: Ipopt itself may still fail to converge ("Restoration Failed",
# dual infeasibility blows up ~1e23 after ~600 iters) — that is a separate NUMERICAL
# issue of the OCP/initial guess, not an autodiff bug.
#
# ⚡ TLDR: `cholesky_solve` (branch-free, replaces `M \ rhs`) + `backend=bypass(:default)`
#   fix both ForwardDiff/LinearAlgebra AD failures. Ipopt non-convergence is a separate
#   numerics issue. See test_ForwardDiff_LinearAlgebra.jl for the MRE.

X = SA[1:8...]
u = SA[9, 10]

# @profview @btime eom(X, u, pars)
