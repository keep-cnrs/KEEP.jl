
using Plots
using StaticArrays
using Test
using NonlinearSolve
using SplitApplyCombine
using Statistics
using Distances: pairwise, Euclidean
using LinearAlgebra: norm, I
using QuasiMonteCarlo: sample, SobolSample
using Setfield: @set

import KEEP.PointMass4 as PM4
using KEEP.PointMassPara
import KEEP.Visualisation as VIS
using KEEP.LimitCycle
import KEEP.LimitCycle: TAU0
using KEEP.TorqueFunction
using KEEP.Optimization

# Include functions for steady state calculations
include("05.4_functions.jl")


# 1. Build parameters
p = build_para()
vbp = build_vbpara(p)

# 2. Initial conditions
τ0, dτ0 = 1e-1, 0
u0 = SA[0.0, τ0, 0.0, dτ0, 0.0]
tf = 60.0

# 3. Integrate the model
sol = PM4.integrate(u0, tf, vbp; save_everystep=true)

# 4. Plot the trajectory
VIS.plot_trajectory_4D(sol; tspan=tf)
plot!(title="4D model, t=0..60s")

# 5. Save the plot
# savefig("my_trajectory_plot.png")


#=
# SECTION: EQUILIBRIA AND LIMIT CYCLES
# Simulations convergent toujours vers one of two cycles limites, par conséquent stables (ou attractifs)
# IC1 donne cycle limite 1, IC2 donne cycle limite 2 -> il existe t tq (1-t) IC1 + t IC2 donne autre chose : un équilibre dans notre cas
# Cherche q tq ddq(q, dq=0) = 0
# Petite vitesse en tau dans la direction du cycle limite souhaité, intégration sur ~10 cycles + callback -> cycle limite
=#

println("\n--- EQUILIBRIA AND LIMIT CYCLES ---")

# Find two limit cycles from different initial conditions
tol = 1e-12
lc_neg = compute_limit_cycle(SA[0.0, TAU0, 0.0, -10, 0.0], vbp; tol=tol, save_everystep=true)
lc_pos = compute_limit_cycle(SA[0.0, TAU0, 0.0, 10, 0.0], vbp; tol=tol, save_everystep=true)

println("Found two limit cycles:")
println("  - Negative cycle (good one) average power: ", lc_neg.u[end][end] / lc_neg.t[end])
println("  - Positive cycle (bad one) average power: ", lc_pos.u[end][end] / lc_pos.t[end])

# Find equilibria
τ_init = 0.5
ss = aligned_and_opposite_steady_states(τ_init, vbp; tol=tol)
println("\nFound equilibrium points.")
@test all(norm.(ddq.(ss, Ref(vbp))) .< tol)

# Start near an equilibrium with a small velocity to find a limit cycle
u0_from_ss = SA[ss[1]..., 0.0, 1e-3, 0.0]
lc_from_ss = compute_limit_cycle(u0_from_ss, vbp; tol=tol, save_everystep=true)
println("\nFound a limit cycle starting from an equilibrium point.")
println("  - Average power: ", lc_from_ss.u[end][end] / lc_from_ss.t[end])


#=
# SECTION: NUMERICAL RESULTS
# ipopt et opti sur une plage de vent de (l, Δθ, Δφ): géométrie et (I, b): dynamique du bras
=#

println("\n--- NUMERICAL RESULTS ---")

# Set up optimization problem
set_torque_function!_function(LINEAR_TORQUE)
# The user mentioned parameters (l, Δθ, Δφ) and (I, b).
# In our model, these correspond to:
# l: tether length, which is `l` in `p`
# Δθ, Δφ: not directly parameters, but related to the motion of the kite.
# I: moment of inertia of the arm, which is `I_eq`
# b: damping of the arm, which is related to `torque_slope`

# We will optimize for `r`, `I_eq`, and `torque_slope` as in `test/07_LimitCycleOptim.jl`.
# `r` is the radius of the drum.
syms = [:r, :I_eq, :torque_slope]
p0 = build_para(torque_slope=100)
lower, upper = make_bounds(p0, syms)
lower[3] = 1
upper[3] = 1e6
tol_optim = 1e-6

println("Running optimization for parameters: ", syms)
stats, model = optimize(p0, syms, lower, upper; tol=tol_optim)

if stats.status == :first_order
    println("Optimization successful!")
    L, M, T = lmt(p0)
    para_dims, state_dims, iterate_dims, power_dim = compute_dims(p0, syms)
    (; solution, objective) = stats
    val = solution[5] * power_dim
    println("  - Optimal average power: ", val, " W")
    sol_params = solution[end-length(syms)+1:end] .* para_dims
    println("  - Optimal parameters: ")
    for (s, v) in zip(syms, sol_params)
        println("    - ", s, ": ", v)
    end
else
    println("Optimization did not converge. Status: ", stats.status)
end

reset_torque_function()
