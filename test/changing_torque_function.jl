using Test
using StaticArrays
using StrFormat

import KEEP.PointMass4 as PM4
using KEEP.PointMassPara
using KEEP.TorqueFunction

vbp = build_vbpara()
u0 = SA{Float64}[0, 1, 2, 3, 0]
tf = 10

reset_torque_function()
sim_pre = PM4.integrate(u0, tf, vbp)

set_torque_function!_function(LINEAR_TORQUE)
sim_post = PM4.integrate(u0, tf, vbp)

@info f"[Avg power before/after]
Default torque function : \%.0f(sim_pre.u[end][5] / tf) W
Linear torque function  : \%.0f(sim_post.u[end][5] / tf) W"

# We used a different torque function that produces a different solution
@test sim_pre.u != sim_post.u

# Resetting the torque function switches back to default = Rational
reset_torque_function()
@test torque_function(1., vbp) == TorqueFunction.torque_function_rational(1., vbp)

# torque_function is type-stable (not type stable for Integer dα because of zero(dα).)
@inferred torque_function(1., vbp)
