using StaticArrays
using Plots
using OrdinaryDiffEqTsit5

using KEEP.PointMassPara
using KEEP.PointMass4
using KEEP.LimitCycle

let
    include("02_EightPlot.jl")
end

# Simulating the cost function without the limit cycle
function build_f(x0, T0, dT0, Pbar, ddPbar)
    function f(x, t)
        T = T0 + dT0 * (x - x0)
        return (Pbar + (x - x0)^2 * ddPbar / 2) + Pbar * T * sin(π * t / T) / (π * t / T)
    end
    return f
end

x0, T0, dT0, pbar, ddpbar = (45, 3, -.2, 50, -1)
f = build_f(x0, T0, dT0, pbar, ddpbar)

x = range(40, 50, 501)
t = range(30, 300, step=.1)

display(plot(t, f.(45, t)))
plot(x, f.(x, 20))
plot!(x, f.(x, 200))
plot!(x, f.(x, 1000))
display(plot!(x, f.(x, 10000)))

x = range(44, 46, 1001)
plot(x, f.(x, 1000))
plot!(x, f.(x, 10000))
display(plot!(x, f.(x, 100000)))


## Optimise using the limit cycle
function gain(x, reference_p)
    # p_dict = Dict(fieldnames(Para) .=> getfield.(Ref(reference_p), fieldnames(Para)))
    # p_dict[:Δφ] = x[1]
    # p = build_para(;p_dict...)
    p = build_para(reference_p; r=x[1], Δθ=x[2], Δφ=x[3], I_eq=x[4])
    ode_prob = ODEProblem(dynamics, SVector(last_initial_state), 1e9, p)
    limit_cycle = compute_limit_cycle(ode_prob)
    last_initial_state .= limit_cycle[begin]
    limit_cycle[end][5]
end

const last_initial_state = MVector(0, TAU0, 0, 0, 0)
const reference_p = build_para()

x0 = [reference_p.r, reference_p.Δθ, reference_p.Δφ, reference_p.I_eq]
gain(x0, reference_p)

using Optimization
using OptimizationOptimJL
using OptimizationNLopt

lb = [15.0, deg2rad(5), deg2rad(5), 1300]
ub = [100.0, deg2rad(70), deg2rad(70), 10000]
optim_prob = OptimizationProblem(gain, x0, reference_p, lb=lb, ub=ub; sense=MaxSense)
sol = solve(optim_prob, NLopt.LN_COBYLA(); maxiters=200)
sol.objective
rad2deg.([sol.minimizer[2], sol.minimizer[3]])
