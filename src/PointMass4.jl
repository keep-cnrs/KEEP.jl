module PointMass4

import LinearAlgebra: norm, normalize, ⋅, ×, dot
using StaticArrays: SA, SMatrix, @SMatrix, StaticArrays
import NaNMath as nm
import OrdinaryDiffEqVerner: Vern9
import ForwardDiff

import KEEP.PointMass4 as PM4
using KEEP.PointMassPara
using KEEP.TorqueFunction
using KEEP.Integrate: _integrate
import KEEP: DEFAULT_TOLERANCE, TAU0

export dynamics, integrate, integrate_verified, build_u, kite_acceleration

@doc raw"""
Compute the jacobian of `f` at `x`.
"""
function jacobian(f, x)
    ## Integration using ForwardDiff is 25% faster compared to AbstractDifferentiation with ForwardDiff backend
    # return AD.jacobian(ab, f, x)[1]
    return ForwardDiff.jacobian(f, x)
end

@doc raw"""
Compute ``D²f(v, v)`` where
 - `f` is a vector valued function,
 - `v1` is the first sampling vector,
 - `v2` is the second sampling vector.

The output is a vector of length `length(f(x))`.
"""
function scnd_deriv_sample(f, x, v1, v2)
    ## Integration using ForwardDiff is 25% faster compared to AbstractDifferentiation with ForwardDiff backend
    ## jacobian-vector product should be done this way
    # return AD.jacobian(ab, τ2 -> AD.jacobian(ab, τ1 -> f(x + τ1 * v1 + τ2 * v2), 0.0)[1], 0.0)[1]
    return ForwardDiff.derivative(τ2 -> ForwardDiff.derivative(τ1 -> f(x + τ1 * v1 + τ2 * v2), 0.0), 0.0)
end

"""
Compute the polar angle θ and azimuthal angle φ given τ. Follows ISO 80000-2.
"""
function τ_to_θφ(τ, p)
    return p.θ0 + p.Δθ * sin(2τ), p.φ0 + p.Δφ * sin(τ)
end

@doc raw"""
Compute the coordinates of ``\hat α``.
"""
function compute_αhat(α)
    sin_α, cos_α = sincos(α)
    return SA[cos_α, sin_α, zero(α)]
end

@doc raw"""
Compute the coordinates of ``\hat τ``.
"""
function compute_τhat(τ, p)
    θ, φ = τ_to_θφ(τ, p)
    sin_θ, cos_θ = sincos(θ)
    sin_φ, cos_φ = sincos(φ)
    return SA[
        sin_θ*cos_φ,
        sin_θ*sin_φ,
        cos_θ
    ]
end


@doc raw"""
(α, τ) -> (R, τ)
"""
function compute_Rτ(q, p)
    L, _, _ = lmt(p)
    α, τ = q
    dot_prod = compute_αhat(α) ⋅ compute_τhat(τ, p)
    # ε = 1e-5
    # my_sqrt(x) = sqrt(sqrt(x^2 + ε^2))
    return SA[L*(dot_prod+nm.sqrt(p.l^2-(1-dot_prod^2))), q[2]]
end

@doc raw"""
Compute the Cartesian coordinates using R and τ.
"""
function compute_OK(Rτ, p)
    R, τ = Rτ
    return R * compute_τhat(τ, p)
end


@doc raw"""
Compute du with automatic differentiation: the reference implementation.

Same model as [`dynamics`](@ref), but every kinematic derivative is obtained by
`ForwardDiff` in the right-hand side (six AD calls per evaluation). It is kept as
the independent check of the closed form, and as the `rhs = :verified` backend of
`KEEP.Optimization.optimize`. Do not delete it: it is what the equivalence test
compares against.

u = (α, τ, dα, dτ)
du = (dα, dτ, ddα, ddτ)
"""
function _dynamics_verified(u, p, t=0)
    L, M, T = lmt(p)
    α, _ = q = u[SA[1, 2]]
    dα, _ = dq = u[SA[3, 4]]

    # Compute values and derivatives
    Rτ = compute_Rτ(q, p)
    OK = compute_OK(Rτ, p)

    D_Rτ = jacobian(q -> compute_Rτ(q, p), q)
    dRτ = D_Rτ * dq
    D_OK = jacobian(Rτ -> compute_OK(Rτ, p), Rτ)

    ddRτ = scnd_deriv_sample(q_ -> compute_Rτ(q_, p), q, dq, dq)
    ddOK = scnd_deriv_sample(Rτ_ -> compute_OK(Rτ_, p), Rτ, dRτ, dRτ)

    # Compute unit vectors
    αhat = compute_αhat(α)
    OA = αhat * L
    rhat = normalize(OK - OA)

    # Compute gravity force
    Fgrav = SA[0, 0, -p.g-p.m_l*p.l] * M * L * T^-2

    # Compute aerodynamical force
    wind = SA[abs(OK[3]/(p.z0*L))^(1/p.n_wind), 0, 0] * L * T^-1
    kite_speed = D_OK * dRτ
    app_wind = wind - kite_speed  # apparent wind, We

    # Apparent wind orthogonal to the lines
    app_wind_perp_to_lines = app_wind - (rhat ⋅ app_wind) * rhat
    ew = normalize(app_wind_perp_to_lines)

    xw = -normalize(app_wind)
    yw = rhat × ew
    zw = xw × yw

    lift_k = norm(app_wind)^2 * (p.c_L * M * L^-1)
    drag_k = lift_k / p.f
    drag_l = (p.c_D_l * M * L^-2) * (p.l * L) * norm(app_wind_perp_to_lines)^2
    Faero_on_kite = -drag_k * xw - lift_k * zw
    Faero_on_lines = drag_l * ew
    Faero = Faero_on_kite + Faero_on_lines

    # torque_function takes a dimensionless dα
    torque = torque_function(dα * T, p) * M * L^2 * T^-2
    l_tension = rhat ⋅ (SA[0, 0, 1] × OA) + 1e-100  # Avoid division by zero for q = (0, 0)
    A_tension = -(p.I * M * L^2) / l_tension * rhat * SA[1, 0]'
    b_tension = -torque / l_tension * rhat

    D_OK_D_Rτ = D_OK * D_Rτ
    A = D_OK_D_Rτ .* M .- A_tension
    b = (ddOK * M +
        D_OK * ddRτ * M +
        -b_tension - Fgrav - Faero
    )
    instantaneous_power = torque * dα
    ddq = (D_OK' * A) \ (-D_OK' * b)
    return SA[dq[1], dq[2], ddq[1], ddq[2], instantaneous_power]
end

@doc raw"""
θ, φ and their first two τ-derivatives (eq. 4 and 10 of `pm2-model-derivation-reworked.tex`).
"""
function _sph_derivs(τ, p)
    θ = p.θ0 + p.Δθ * sin(2τ)
    φ = p.φ0 + p.Δφ * sin(τ)
    θp = 2p.Δθ * cos(2τ)
    θpp = -4p.Δθ * sin(2τ)
    φp = p.Δφ * cos(τ)
    φpp = -p.Δφ * sin(τ)
    return θ, φ, θp, θpp, φp, φpp
end

@doc raw"""
Explicit kinematics of the cone constraint, eq. (6)–(30): `R = r(α,τ)`, its first
derivatives `R_α, R_τ`, its Hessian `HR`, and the lemniscate frame
`rhat, rhatp, rhatpp`. Closed form, no automatic differentiation.
"""
function _kinematics(α, τ, p)
    L, _, _ = lmt(p)
    A = L
    ELL = L * p.l
    θ, φ, θp, θpp, φp, φpp = _sph_derivs(τ, p)
    sψ, cψ = sincos(φ - α)
    sθ, cθ = sincos(θ)
    sφ, cφ = sincos(φ)

    p_ = A * sθ * cψ
    Δ = nm.sqrt(p_^2 + ELL^2 - A^2)
    R = p_ + Δ

    p_α = A * sθ * sψ
    p_τ = A * (θp * cθ * cψ - φp * sθ * sψ)
    p_αα = -p_
    p_ατ = A * (θp * cθ * sψ + φp * sθ * cψ)
    p_ττ = A * ((θpp * cθ - (θp^2 + φp^2) * sθ) * cψ -
                2θp * φp * cθ * sψ - φpp * sθ * sψ)

    Rp = R / Δ
    R_α = Rp * p_α
    R_τ = Rp * p_τ
    Rpp = (ELL^2 - A^2) / Δ^3
    HR = @SMatrix [Rp*p_αα+Rpp*p_α^2  Rp*p_ατ+Rpp*p_α*p_τ
                   Rp*p_ατ+Rpp*p_α*p_τ Rp*p_ττ+Rpp*p_τ^2]

    rhat = SA[sθ * cφ, sθ * sφ, cθ]
    ahat = SA[cos(α), sin(α), zero(α)]
    rhatp = SA[θp * cθ * cφ - φp * sθ * sφ,
        θp * cθ * sφ + φp * sθ * cφ,
        -θp * sθ]
    nrhatp2 = θp^2 + sθ^2 * φp^2
    eθ = SA[cθ * cφ, cθ * sφ, -sθ]
    eφ = SA[-sφ, cφ, zero(θ)]
    rhatpp = -nrhatp2 * rhat + (θpp - φp^2 * sθ * cθ) * eθ +
             (φpp * sθ + 2θp * φp * cθ) * eφ

    return (; R, Δ, HR, R_α, R_τ, rhat, ahat, rhatp, rhatpp, nrhatp2,
        θ, φ, θp, φp, A, ELL)
end

"Forces: aerodynamics on kite and lines, weight, and the generator torque."
function _forces(k, dq, p)
    L, M, T = lmt(p)
    OK = k.R * k.rhat
    OA = k.A * k.ahat
    lhat = (OK - OA) / k.ELL                 # loop closure: ‖OK - OA‖ = ℓ exactly

    F_g = SA[0, 0, -p.g - p.m_l * p.l] * M * L * T^-2

    Rdot = k.R_α * dq[1] + k.R_τ * dq[2]
    D_OK = @SMatrix [k.rhat[1] k.R*k.rhatp[1]
                     k.rhat[2] k.R*k.rhatp[2]
                     k.rhat[3] k.R*k.rhatp[3]]
    kite_speed = D_OK * SA[Rdot, dq[2]]

    wind = SA[abs(OK[3] / (p.z0 * L))^(1 / p.n_wind), 0, 0] * L * T^-1
    app_wind = wind - kite_speed
    perp = app_wind - (lhat ⋅ app_wind) * lhat
    ew = perp / norm(perp)
    xw = -app_wind / norm(app_wind)
    yw = lhat × ew
    zw = xw × yw

    lift_k = norm(app_wind)^2 * (p.c_L * M * L^-1)
    drag_k = lift_k / p.f
    drag_l = (p.c_D_l * M * L^-2) * (p.l * L) * norm(perp)^2
    F_a = (-drag_k * xw - lift_k * zw) + drag_l * ew

    torque = torque_function(dq[1] * T, p) * M * L^2 * T^-2
    return F_a, F_g, torque, Rdot
end

@doc raw"""
Compute du, the evolution of the state — closed-form implementation (default).

Same model as [`_dynamics_verified`](@ref), with the explicit formulas of
`applications/jtcam/scratch/pm2-model-derivation-reworked.tex` instead of six
`ForwardDiff` calls per evaluation: **2.33× cheaper per evaluation**, 2.1× on a
limit-cycle solve, identical to the verified version to 1e-15 (`test/04a`, run by
`verify_fast_rhs.jl` in the benchmark).

u = (α, τ, dα, dτ)
du = (dα, dτ, ddα, ddτ)
"""
function dynamics(u, p, t=0)
    L, M, _ = lmt(p)
    α, τ = u[1], u[2]
    dα, dτ = u[3], u[4]
    dq = SA[dα, dτ]

    k = _kinematics(α, τ, p)
    F_a, F_g, torque, Rdot = _forces(k, dq, p)

    # + 1e-100 mirrors the guard on `l_tension` in `_dynamics_verified`: μ vanishes on
    # the section seed (α = 0, τ = TAU0), where the arm balance is a removable singularity
    mu = k.A * k.R / k.ELL * sin(k.θ) * sin(k.φ - α) + 1e-100
    κ = p.I * M * L^2 / mu
    Mmat = @SMatrix [M*k.R_α+κ*k.Δ/k.ELL  M*k.R_τ
                     -κ*k.Δ/k.ELL*k.R_τ   M*k.R^2*k.nrhatp2]

    h1 = -M * (dot(dq, k.HR * dq) - k.R * dτ^2 * k.nrhatp2) +
         k.rhat ⋅ F_a + k.rhat ⋅ F_g - torque * k.Δ / (mu * k.ELL)
    h2 = -M * (k.R * (2Rdot * dτ * k.nrhatp2 + k.R * dτ^2 * (k.rhatp ⋅ k.rhatpp))) +
         (k.R * k.rhatp) ⋅ F_a + (k.R * k.rhatp) ⋅ F_g +
         torque * k.Δ * k.R_τ / (mu * k.ELL)

    ddq = Mmat \ SA[h1, h2]
    return SA[dα, dτ, ddq[1], ddq[2], torque * dα]
end


"""
[For post-processing]

Given u, p and possibly du, compute the line tension
"""
function compute_line_tension(u, p; du=PM4.dynamics(u, p, 0))
    L, M, T = lmt(p)
    α = u[1]
    dα = u[3]
    ddα = du[3]
    αhat = PM4.compute_αhat(α)
    OA = p.a * αhat
    OK = PM4.compute_OK(PM4.compute_Rτ(u, p), p)
    rhat = normalize(OK - OA)
    torque = torque_function(dα * T, p) * M * L^2 * T^-2
    l_tension = rhat ⋅ (SA[0, 0, 1] × OA) + 1e-100
    A_tension = -p.I * M * L^2 / l_tension * rhat * SA[1, 0]'
    b_tension = -torque / l_tension * rhat
    F_tension = A_tension[:, 1] * ddα + b_tension
    return -F_tension ⋅ rhat
end

@doc raw"""
Cartesian acceleration ``\ddot{x}`` of the kite, with ``x = OK(τ, α)``.

`u = (α, τ, dα, dτ, W)` is a state and `p` a `VBPara`. The angular
accelerations come from `dynamics(u, p, 0)`. Positions from `compute_OK` are in
metres and time is in seconds — `build_vbpara` keeps `a`, `m` and `w0`
physical — so the result is in m s⁻².
"""
function kite_acceleration(u, p)
    q = u[SA[1, 2]]
    dq = u[SA[3, 4]]
    ddq = dynamics(u, p, 0)[SA[3, 4]]
    pos(q_) = compute_OK(compute_Rτ(q_, p), p)
    return scnd_deriv_sample(pos, q, dq, dq) + jacobian(pos, q) * ddq
end

function build_u(α, τ, dα, dτ, W)
    return SA[α, τ, dα, dτ, W]
end

function build_u(; α=0, dα=0, dτ=0, τ=TAU0, W=0.0)
    return build_u(α, τ, dα, dτ, W)
end

function integrate(u0, tf, p, alg=Vern9(); save_everystep=false, tol=DEFAULT_TOLERANCE, kwargs...)
    return _integrate(dynamics, u0, tf, p, alg; save_everystep=save_everystep, tol=tol, kwargs...)
end

"Integrate with the AD reference dynamics ([`_dynamics_verified`](@ref))."
function integrate_verified(u0, tf, p, alg=Vern9(); save_everystep=false, tol=DEFAULT_TOLERANCE, kwargs...)
    return _integrate(_dynamics_verified, u0, tf, p, alg; save_everystep=save_everystep, tol=tol, kwargs...)
end



end  # module