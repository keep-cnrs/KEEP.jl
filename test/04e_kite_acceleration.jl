using Test
using LinearAlgebra

using KEEP: PointMass4 as PM4
using KEEP.PointMassPara

p = build_vbpara(build_para())
u0 = PM4.build_u(; α=0.3, τ=0.5, dα=0.4, dτ=1.0)
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
    states = vec([
        PM4.build_u(α, τ, dα, dτ, 0.0) for α in (-π, -1.0, 0.0, 1.0, 2.5),
        τ in (0.0, 0.7, 1.9, 3.5), dα in (-2.0, -0.3, 0.8),
        dτ in (-1.5, 0.4, 2.0)
    ])
    # μ = 0 (sin(φ−α) = 0) is a removable singularity of the arm balance: there the AD
    # reference divides by `l_tension + 1e-100` and returns O(1e83) while the closed form
    # returns 0. No trajectory reaches that line, so equivalence is asserted off it, plus
    # at the bootstrap seed, where the singular term carries a factor dα and vanishes.
    μline(u) = abs(sin(p.φ0 + p.Δφ * sin(u[2]) - u[1])) < 1e-9
    states = filter(u -> !μline(u), states)
    u_seed = PM4.build_u(α=0.0, τ=TAU0, dα=0.0, dτ=1.0)
    push!(states, u_seed)
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
