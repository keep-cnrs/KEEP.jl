using Test
using LinearAlgebra

using KEEP: PointMass4 as PM4
using KEEP.PointMassPara

p = build_vbpara(build_para())
u0 = PM4.build_u(α=0.3, τ=0.5, dα=0.4, dτ=1.0)
h = 1e-4

# Central difference of the Cartesian kite position along the true flow.
sol = PM4.integrate(u0, 2h, p; save_everystep=true, tol=1e-10)
pos(u) = PM4.compute_OK(PM4.compute_Rτ(u[1:2], p), p)
a = PM4.kite_acceleration(sol(h), p)
a_fd = (pos(sol(2h)) - 2pos(sol(h)) + pos(sol(0.0))) / h^2

@test a ≈ a_fd rtol = 1e-6

## The closed-form `dynamics` must reproduce the AD reference `_dynamics_verified`.
## Appended here because this file already validates the same right-hand side.
import KEEP: TAU0
@testset "closed-form dynamics == _dynamics_verified" begin
    states = vec([PM4.build_u(α, τ, dα, dτ, 0.0)
                  for α in (-π, -1.0, 0.0, 1.0, 2.5), τ in (0.0, 0.7, 1.9, 3.5),
                  dα in (-2.0, -0.3, 0.8), dτ in (-1.5, 0.4, 2.0)])
    # Both forms solve the arm balance for ddα, never for the line tension, so they are
    # regular on the fold curves α = φ(τ), φ(τ) + π, where the lever arm of the lines
    # vanishes (and where all the equilibria lie): include states exactly on them.
    φ(τ) = p.φ0 + p.Δφ * sin(τ)
    append!(states, [PM4.build_u(φ(τ) + s, τ, dα, dτ, 0.0)
                     for τ in (0.0, 0.7, 1.9, 3.5), s in (0.0, π), dα in (-0.3, 0.0, 0.8), dτ in (-1.5, 0.0, 2.0)])
    push!(states, PM4.build_u(α=0.0, τ=TAU0, dα=0.0, dτ=1.0))
    err = 0.0
    for u in states
        y_ref = PM4._dynamics_verified(u, p)
        y_new = PM4.dynamics(u, p)
        @test all(isfinite, y_new)
        err = max(err, maximum(abs.(y_new .- y_ref)) / max(maximum(abs, y_ref), eps()))
    end
    @test err < 1e-12

    # and the same trajectory over a nontrivial arc
    u0 = PM4.build_u(α=0.4, τ=0.7, dα=0.3, dτ=-0.2)
    a = PM4.integrate(u0, 0.8, p; tol=1e-12).u[end]
    b = PM4.integrate_verified(u0, 0.8, p; tol=1e-12).u[end]
    @test isapprox(a, b; rtol=1e-9)
end

## On the fold curves the arm is free (the lines exert no torque): at rest, ddα = 0, and
## the dynamics and the tension are continuous across the fold.
@testset "regularity on the fold curves" begin
    φ(τ) = p.φ0 + p.Δφ * sin(τ)
    for τ in range(-π, π, 37), s in (0.0, π)
        @test abs(PM4.dynamics(PM4.build_u(φ(τ) + s, τ, 0.0, 0.0, 0.0), p)[3]) < 1e-12
        u0 = PM4.build_u(φ(τ) + s, τ, 0.7, 1.8, 0.0)
        u1 = PM4.build_u(φ(τ) + s + 1e-7, τ, 0.7, 1.8, 0.0)
        y0, y1 = PM4.dynamics(u0, p), PM4.dynamics(u1, p)
        @test all(isfinite, y0)
        @test maximum(abs, y0 - y1) < 1e-4 * max(1, maximum(abs, y0))
        F0, F1 = PM4.compute_line_tension(u0, p), PM4.compute_line_tension(u1, p)
        @test isfinite(F0) && abs(F0 - F1) < 1e-4 * max(1, abs(F0))
    end
end
