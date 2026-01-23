using Test
using StaticArrays
using StrFormat

import KEEP.PointMass4 as PM4
using KEEP.PointMassPara
using KEEP.TorqueFunction

vbp = build_vbpara()
u0 = SA{Float64}[0, 1, 2, 3, 0]
tf = 10

reset_torque_function!()
sim_default = PM4.integrate(u0, tf, vbp)

set_torque_function!(LINEAR_TORQUE)
sim_linear = PM4.integrate(u0, tf, vbp)

set_torque_function!(RATIONAL_TORQUE)
sim_rational = PM4.integrate(u0, tf, vbp)

reset_torque_function!()
sim_reset = PM4.integrate(u0, tf, vbp)

@info f"[Avg power before/after]
Default torque function : \%.0f(sim_default.u[end][5] / tf) W
Linear torque function  : \%.0f(sim_linear.u[end][5] / tf) W
Rational torque function: \%.0f(sim_rational.u[end][5] / tf) W
Reset torque function   : \%.0f(sim_reset.u[end][5] / tf) W"

# We used a different torque function that produces a different solution
@test sim_linear.u != sim_rational.u

# Resetting the torque function switches back to default
@test sim_default.u == sim_reset.u

# torque_function is type-stable (not type stable for Integer dα because of zero(dα).)
@inferred torque_function(1., vbp)
