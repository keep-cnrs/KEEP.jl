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
