using Test
using StrFormat
using Setfield: @set
using Plots
import ComponentArrays: ComponentArray as CA, getaxes
using Printf

import KEEP.PointMass4 as PM4
using KEEP.PointMassPara
using KEEP.LimitCycle: shoot, endpoint_residuals
using KEEP.TorqueFunction: torque_function, set_torque_function!, reset_torque_function!, LINEAR_TORQUE
using KEEP.Optimization: compute_dims, make_bounds, optimize
using KEEP.Visualization

reset_torque_function!()
include("utils.jl")
tol = 1e-6

default(lw=3, formatter=:plain, label="")

set_torque_function!(LINEAR_TORQUE)

# set_torque_function!(RATIONAL_TORQUE)
# syms = [:r, :I_eq, :Cmax, :Ωmax]

p0 = build_para()

## Step 1
syms0 = [:r, :I_eq, :torque_slope]
lower0, upper0 = make_bounds(p0, syms0)
lower0[3] = 1
upper0[3] = 1e6
solution1, stats, model = optimize(p0, syms0, lower0, upper0; tol=tol, max_wall_time=500.);  # High wall_time for github CI
p1 = CA(p0; solution1.params...)

syms = syms0
lower = lower0
upper = upper0
solution = solution1
p = p1

## Step 2
if false
    syms1 = [:r, :I_eq, :torque_slope, :Δθ, :Δφ]
    lower1, upper1 = make_bounds(p0, syms1)
    lower1[3] = 1
    upper1[3] = 1e6
    solution2, stats, model = optimize(p1, syms1, lower1, upper1; tol=tol, max_wall_time=500., initial_guess=solution1)  # High wall_time for github CI
    p2 = CA(p1; solution2.params...)

    syms = syms1
    lower = lower1
    upper = upper1
    solution = solution2
    p = p2
end

## Avant/Après : texte

# solution  = [α0, dα0, dτ0, tf, Wf / tf, optim_paras...]
# L, M, T = lmt(p)
para_dims, state_dims, iterate_dims, power_dim = compute_dims(p0, syms)
x0 = CA(model.meta.x0 .* iterate_dims, getaxes(solution))
(; α0, dα0, dτ0, tf, P, params) = solution

_cons_buffer = similar(model.meta.ucon)
constraints_neg_log_error = -log10(sqrt(sum(model.c!(_cons_buffer, stats.solution) .^ 2)))

init_sol = x0.params
sol = solution.params

init_val = x0.P
val = solution.P
@test stats.status == :first_order

# Build the compact padded table
# 1. Gather all table data into an array of Tuples: (Symbol, Low, High, Init, Sol)
table_data = [
    ("tf", 1e-3, 100.0, x0.tf, solution.tf),
    ((string(s), lower[i], upper[i], x0.params[s], solution.params[s]) for (i, s) in enumerate(syms))...
]

# 2. Compile row formatter and build the table as a single array of strings
row_fmt = Printf.Format("%-14s | %10.3g | %10.3g | %10.3g | %10.3g")
table_str = join([
        @sprintf("%-14s | %10s | %10s | %10s | %10s", "Symbol", "Low", "High", "Init", "Sol"),
        join(["-"^w for w in (14, 10, 10, 10, 10)], "-+-"),
        (Printf.format(row_fmt, row...) for row in table_data)...
    ], "\n")

# Log the result
@info f"""Status: \%(string(stats.status))
Constraints neglog error (should be 0 < e < 16): \%.1f(constraints_neg_log_error) (~correct digits)
dτ0 (should be > 0): \%.3f(dτ0)

Average power:
\tInit: \%.1f(x0.P) W
\tSol : \%.1f(solution.P) W

Variables:
\%(table_str)
"""



# Two ways to construct the vbpara
vbp0 = build_vbpara(p0)
vbp = build_vbpara(p)

## Compare initial guess and solution

initial_guess_shooting = x0[1:4]
initial_guess_sim = shoot(initial_guess_shooting, vbp0, tol=tol, save_everystep=true)
@test sum(abs, endpoint_residuals(initial_guess_sim)) < 8tol  # 4 is length of residuals, 2 is a "safety" factor, 2*4=8


# Construct the solution's limit cycle (from the shooting + the parameters), check that we get a limit cycle and that the power is correct
shooting = solution[1:4]
solution_sim = shoot(shooting, vbp, tol=tol, save_everystep=true)
@test sum(abs, endpoint_residuals(solution_sim)) < 8tol

# Parameters for plotting
N = 1000
t_init = range(extrema(initial_guess_sim.t)..., length=N)
t_sol = range(extrema(solution_sim.t)..., length=N)

## Power plotting

function power(u, vbp)
    L, M, T = lmt(vbp)
    dα = u[3]
    return dα * torque_function(dα * T, vbp) * (M * L^2 * T^-2)
end

initial_guess_power = power.(initial_guess_sim.(t_init), Ref(vbp0))
initial_guess_avg_power = initial_guess_sim.u[end][end] / initial_guess_sim.t[end]
solution_power = power.(solution_sim.(t_sol), Ref(vbp))
solution_avg_power = solution_sim.u[end][end] / solution_sim.t[end]
max_t = max(initial_guess_sim.t[end], solution_sim.t[end])
max_P = max(maximum(initial_guess_power), maximum(solution_power))

plot(xlabel="t (s)", ylabel="Power (W)", legend=:outerbottom, title="Power")
plot!(t_init, initial_guess_power, c=:grey75)
plot!(t_sol, solution_power, c=:orange)
hline!([initial_guess_avg_power], c=:grey50, ls=:dash)
hline!([solution_avg_power], c=:red, ls=:dash)
annotate!(max_t, initial_guess_avg_power - max_P * 0.05, text(f"\%.0f(initial_guess_avg_power) W", :right))
annotate!(max_t, solution_avg_power - max_P * 0.05, text(f"\%.0f(solution_avg_power) W", :right))
plot!(1:0, c=:grey, label="Initial guess")
plot!(1:0, c=:red, label="After optimising")
display(plot!())


pow_init_fig = plot_avg_power_4D(initial_guess_sim)
pow_sol_fig = plot_avg_power_4D(solution_sim)
max_y = maximum(fig.subplots[1].attr[:yaxis].plotattributes[:extrema].emax for fig in (pow_init_fig, pow_sol_fig))
for fig in (pow_init_fig, pow_sol_fig)
    fig.subplots[1].attr[:yaxis].plotattributes[:extrema].emax = max_y
end
plot!(pow_init_fig, title="Initial guess")
plot!(pow_sol_fig, title="After optimising")
plots = [pow_init_fig pow_sol_fig]
display(plot(plots..., size=500 .* size(plots'), layout=size(plots)))
# savefig(plot!(), "docs/media/optim_initial_guess_power.pdf")
# savefig(plot!(), "docs/media/optim_solution_power.pdf")

## Arm position and velocity checks
α_init = [initial_guess_sim(t, idxs=1) for t in t_init]
dα_init = [initial_guess_sim(t, idxs=3) for t in t_init]
α_sol = [solution_sim(t, idxs=1) for t in t_sol]
dα_sol = [solution_sim(t, idxs=3) for t in t_sol]

plot(xlabel="α (rad)", ylabel="dα (rad/s)", legend=:outerbottom, legend_columns=2, title="Arm position and velocity")
plot!(α_init, dα_init, c=:grey75, label="Initial guess")
plot!(α_sol, dα_sol, idxs=(1, 3), c=:orange, label="After optimisation")
# hline!([-vbp0.Ωmax, vbp0.Ωmax], c=:grey50, ls=:dash)
hline!([-vbp.Ωmax, vbp.Ωmax], c=:red, ls=:dash)
vline!([-π / 2, π / 2], c=:black, label="-π/2 < α < π/2")
display(plot!())

using ForwardDiff: derivative

_, _, T = lmt(vbp)

radpersec2rpm = 60 / 2π
α_init_max = rad2deg(maximum(α_init))  # °
α_sol_max = rad2deg(maximum(α_sol))  # °
dα_init_max = maximum(dα_init) * radpersec2rpm
dα_sol_max = maximum(dα_sol) * radpersec2rpm
Ωmax_init = vbp0.Ωmax / T * radpersec2rpm
Ωmax_sol = vbp.Ωmax / T * radpersec2rpm

dtorque_dω_dim = power_dim * T^2
dtorque_dω0_init = derivative(ω -> torque_function(ω, vbp0), 0) * dtorque_dω_dim
dtorque_dω0_sol = derivative(ω -> torque_function(ω, vbp), 0) * dtorque_dω_dim

@info f"""Arm position and velocity
Initial guess:
\tmax(α) = \%.2f(α_init_max)° < 90° = \%s(α_init_max < 90)
\tmax(dα) (\%.2f(dα_init_max) rpm) < Ωmax (\%.2f(Ωmax_init) rpm) : \%s(dα_init_max < Ωmax_init)
\td(torque)/dω at ω = 0 : \%i(dtorque_dω0_init) kg m^2 s^-1

Solution:
\tmax(α) = \%.2f(α_sol_max)° < 90° = \%s(α_sol_max < 90)
\tmax(dα) (\%.2f(dα_sol_max) rpm) < Ωmax (\%.2f(Ωmax_sol) rpm) : \%s(dα_sol_max < Ωmax_sol)
\td(torque)/dω at ω = 0 : \%i(dtorque_dω0_sol) kg m^2 s^-1
"""

## Line tension comparison
tension_init = PM4.compute_line_tension.(initial_guess_sim.u, Ref(vbp0))
tension_sol = PM4.compute_line_tension.(solution_sim.u, Ref(vbp))

plot(title="Line tension", xlabel="Time (s)", ylabel="Tension (N)")
plot!(initial_guess_sim.t, tension_init, c=:grey75, label="Initial guess")
plot!(solution_sim.t, tension_sol, c=:orange, label="After optimisation")
# plot!(lc.t, tension_lc, c=:lightgrey, label="Other limit cycle")
plot!([0.], [0.], alpha=0, label="", yformatter=:plain)
display(plot!())

@info f"""Line tension comparison
Before optim:
\tMin tension: \%.2f(minimum(tension_init)) N
\tMax tension: \%.2f(maximum(tension_init)) N

After optim:
\tMin tension: \%.2f(minimum(tension_sol)) N
\tMax tension: \%.2f(maximum(tension_sol)) N
"""

## Phase evolution comparison

function compute_state_scale(sim)
    return maximum.(abs, eachrow(reduce(hcat, sim.u)))
end

function plot_trajectory(sim; time_scale=1, state_scale=compute_state_scale(sim), N=100, fig=plot(), kwargs...)
    t = range(sim.t[1], sim.t[end], length=N)
    u = eachrow(reduce(hcat, sim.(t)))
    plot!(fig, t / time_scale, u ./ state_scale; kwargs...)
end

state_scale = compute_state_scale(solution_sim)
plot_trajectory(solution_sim, color=(1:5)', N=1000, labels=["α" "τ" "dα" "dτ" "W"])
plot_trajectory(initial_guess_sim, state_scale=state_scale, fig=plot!(), color=(1:5)', ls=:dash, N=1000, label="")
plot!(1:0, ls=:dash, label="Initial guess", c=:black, title="Normalized phase vs time (s)")
display(plot!())


## Animation
function anim_2d(sim; fps=20, kwargs...)
    vbp = sim.prob.p
    return @animate for t in range(0, sim.t[end], length=round(Int, fps * sim.t[end]))[1:end-1]
        plot_eight_circle([sim(t)], vbp; plot_kwargs=(label="", color=:blue, kwargs...))
        # p_plot = init_plot_eight_circle(vbp)
        # add_point_plot_eight_circle!(sim(t), p_plot, label="", color=:blue, kwargs...)
    end
end

fps = 15
if true  # creating the animations somehow opens a window at each frame
    flat_anim_init = anim_2d(initial_guess_sim, fps=fps)
    flat_anim_sol = anim_2d(solution_sim)
    display(gif(flat_anim_init, fps=fps))
    display(gif(flat_anim_sol, fps=fps))
end

# gif(flat_anim_init, "docs/media/optim_initial_guess_flat_anim.gif", fps=fps)
# gif(flat_anim_sol, "docs/media/optim_solution_flat_anim.gif")

reset_torque_function!()


## Sensitivity analysis
# Local
#=
using ForwardDiff: JacobianConfig

if false
    vbp = build_vbpara()
    cfg = JacobianConfig(p -> optimize(p, args), vbp)
    vbp_ = eltype(cfg.duals).(vbp, cfg.seeds)  # Why not same length ?

    @time stats, model = optimize(vbp, args...);
end
=#

# Landscape around the optimal point
#=
using GlobalSensitivity

f = (x) -> .5*(x[1]^2 - x[2]^2/10) + x[1]^2 * x[2]

bounds = [(-1, 1), (-1, 1)]
N = 100

m = gsa(f, Sobol(), bounds, samples=N)
bar(m.S1, ylim = (0, 1))
bar(m.ST, ylim = (0, 1))
=#


# Benchmarking
#=
import ForwardDiff
using BenchmarkTools

x = copy(stats.solution)
x[4] = 1000
println("#in = ", length(x), ", #out = ", length(_cons_buffer))
@btime model.c!(_cons_buffer, x);  # value call
_Jcons_buffer = zeros(length(x), length(_cons_buffer))
@time @btime ForwardDiff.jacobian!(_Jcons_buffer, model.c!, _cons_buffer, x);  # Full jacobian

function jvp!(jvp_buffer, f!, x, v)
    g!(jvp_buffer, t) = f!(jvp_buffer, x + t * v)
    ForwardDiff.derivative!(jvp_buffer, g!, jvp_buffer, 0)
end

@btime jvp!(_cons_buffer, model.c!, x, ones(length(x)));
=#

## Finding the phase lag of the arm, close to π/4
# s = solution_sim
# φ = π/4
# plot(s.t, s[1:4, :]', label=["α" "τ" "dα" "dτ"])
# plot!(s.t, π/2 * sin.(s[2, :] .+ φ), label="sin(τ + φ)", lc=:black)

# function resid(φ, s)
#     t = s.t
#     α = s(t, idxs=1)
#     τ = s(t, idxs=2)
#     α_approx = maximum(α) * sin.(τ .+ φ)
#     return sum(abs, α - α_approx) / length(t)
# end

# plot(φ -> resid(φ, s), -π, π, label="Residual")

# argmin(φ -> resid(φ, s), LinRange(-π, π, 1000))
