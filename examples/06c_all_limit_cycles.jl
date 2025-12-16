using Logging
using Test

using Setfield
using ProgressLogging
using Plots
using ComponentArrays
using OrdinaryDiffEqTsit5
using DiffEqCallbacks
using StaticArrays
using Statistics
using LinearAlgebra: norm
using Distances: pairwise, Euclidean
using QuasiMonteCarlo: sample, SobolSample

import KEEP.PointMass4 as PM4
using KEEP.PointMassPara
using KEEP.LimitCycle: TAU0, compute_limit_cycle
import KEEP: DEFAULT_TOLERANCE

include("00_general_functions.jl")
include("05.4_functions.jl")

function nan_norm(x, p=2)
    length(x) == 1 && return abs(x[1])
    return nm.pow(nm.sum(nm.pow.(nan_norm.(x, p), p)), 1/p)
end

function compute_limit_cycle_2callbacks(u0, vbp; nb_cycles=5)
    dτ_sign = ifelse(u0[4] > 0, +, -)
    τstart = 2π * nb_cycles
    τstop = τstart + 2π
    cb1 = ContinuousCallback((u, t, integrator) -> abs(u[2] - TAU0) - τstart, identity, save_positions=(true, false))
    cb2 = ContinuousCallback((u, t, integrator) -> abs(u[2] - TAU0) - τstop, terminate!, save_positions=(true, false))
    cb = CallbackSet(cb1, cb2)

    sol = PM4.integrate(u0, 1e4, vbp, callback=cb, save_start=false, save_end=false)

    if length(sol.u) < 2
        return Inf, SA[NaN, NaN, NaN, NaN, NaN]
    end

    yₙ = sol.u[end]
    yₙ₋₁ = sol.u[end-1]
    T = sol.t[end] - sol.t[end-1]
    y_T = SA[yₙ[1], TAU0, yₙ[3], yₙ[4], yₙ[5] - yₙ₋₁[5]]
    return T, y_T
end

function compute_limit_cycle_2callbacks(vbp, dτ_sign::Union{typeof(+),typeof(-)}; τ_init=1)
    q0 = steady_state(τ_init, vbp)
    u0 = SA[q0..., 0, dτ_sign(.1), 0]
    return compute_limit_cycle_2callbacks(u0, vbp)
end

tol = DEFAULT_TOLERANCE
vbp = build_vbpara()

# function observe(u)
#     α, τ, dα, dτ, _ = u
#     return [cos(α), sin(α), cos(τ), sin(τ), dα, dτ]
# end

## Method 1: first compute the steady states, then add a small speed and integrate
dτ_sign = -
τmin, τmax = 0, π
αmin, αmax = -π, π
N = 100
qs = eachcol(sample(N, [αmin, τmin], Float64[αmax, τmax], SobolSample()))
ss = reduce_expand_reduce(combinedims([steady_state(α, τ, vbp; tol=1e-15) for (α, τ) in qs]), tol)
[PM4.dynamics(SA[q..., 0, 0, 0], vbp, 0) for q in eachcol(ss)]

u0s = [ss; repeat([0, dτ_sign(.1), 0], inner=(1, size(ss, 2)))]
u0s
Ts, y_Ts = invert([compute_limit_cycle_2callbacks(SA[u0...], vbp) for u0 in eachcol(u0s)])
@test norm(pairwise(Euclidean(), Ts)) < sqrt(tol)  # If passed, all equilibriums end on the same cycle, for a given sign of dτ

## Method 2: Sample the phase space and directly integrate
# Pour les bornes, cf. https://github.com/ljad-cnrs/keep/issues/78#issuecomment-2714484140
s = 4  # Security factor, leeway around the approximated speeds
Δα = π
T = 2 * vbp.r * vbp.Δφ / vbp.v_ref  # Approximation grossière de la période
dαmin, dαmax = s * Δα / T * [-1, 1]
dτmin, dτmax = 2π / T * [1 / s, s]
q_dqs = let
    samples = sample(N, [αmin, τmin, dαmin, dτmin], [αmax, τmax, dαmax, dτmax], SobolSample())
    samples = [samples [1, 1, 1, -1] .* samples]
    eachcol(samples)
end
Ts, y_Ts = invert(@progress x = [compute_limit_cycle_2callbacks(SA[q_dq..., 0], vbp, nb_cycles = 5) for q_dq in q_dqs])
# @test norm(pairwise(Euclidean(), Ts)) < sqrt(tol) # If passed, there is only one attracting limit cycle

# TODO: plots = α/τ, équilibres et trajectoires (depuis équilibres pour un plot, depuis partout pour l'autre)

PM4.integrate(SA[q_dqs[1]..., 0], -.1, vbp)

function drift(q, t, vbp)
    return q - PM4.integrate()
end

u0 = SA[ss[:, 1]..., 0, 0, 0]
ε = 1e-5
tf = 15

for sense in (+, -)
    u0_ = @set u0[4] = sense(ε)
    u = PM4.integrate(u0_, tf, vbp; save_everystep=true)
    plot(u.t, invert(u.u)[1:4], title="steady state with dτ = $(sense(1e-5))", label=["α" "τ" "dα" "dτ"])
    display(plot!())
end

# TODO: A réécrire sans dynamicalsystems si besoin
# using Attractors, DynamicalSystems

# function observe(state)
#     @assert length(state) == 5
#     α, τ, dα, dτ, W = state
#     # return SA[cos(α), cos(τ)]
#     return SA[cos(α), cos(τ), sin(α), sin(τ), dα, dτ]
# end

# function unobserve(observed_state)
#     @assert length(observed_state) == 6
#     # cα, cτ = observed_state
#     # return SA[acos(cα), acos(cτ), 0, 0, 0]
#     cα, cτ, sα, sτ, dα, dτ = observed_state
#     return SA[atan(cα, sα), atan(cτ, sτ), dα, dτ, 0]
# end

# diffeq = (alg = Tsit5(), abstol = tol, reltol = tol)
# ds = CoupledODEs(PM4.dynamics, u0, vbp; diffeq)
# projected_ds = ProjectedDynamicalSystem(ds, observe, unobserve)


# # Plot cos(τ) vs cos(α) pour avoir un phase space torique
# for sense in (+, -)
#     reinit!(ds, @set u0[4] = sense(ε))  # Selon ε>0 ou ε<0, on voit que la spirale part vers le haut ou le bas et on tombe dans le bassin correspondant
#     U, t = trajectory(ds, 5, Δt = .01)  # t assez petit pour voir le comportement proche de l'équilibre
#     plot(columns(U)[1:2], xlabel="α", ylabel="τ", title="(α, τ) with dτ = $(sense(1e-5))")
#     scatter!(eachcol(U[1:1, 1:2])..., marker=:x, label="start")
#     display(plot!())

#     reinit!(projected_ds)
#     U, t = trajectory(projected_ds, tf, Δt = .01)
#     plot(columns(U)[1:2], xlabel="cos(α)", ylabel="cos(τ)", title="'torus visualisation' of the phase space with dτ = $(sense(1e-5))")
#     display(plot!())
# end


# featurizer(U, t) = SA[mean.(columns(U)[end-1:end])...]

# cαmin, cαmax = cτmin, cτmax = sαmin, sαmax = sτmin, sτmax = -1, 1
# dαmin, dαmax = dτmin, dτmax = -10, 10
# N = 101
# grid = (
#     range(cαmin, cαmax, length=N),
#     range(cτmin, cτmax, length=N),
#     range(sαmin, sαmax, length=N),
#     range(sτmin, sτmax, length=N),
#     range(dαmin, dαmax, length=N),
#     range(dτmin, dτmax, length=N),
# )
# sampler, = statespace_sampler(grid)

# mapper = AttractorsViaFeaturizing(projected_ds, featurizer, T=5, Ttr=5, Δt = .01)
# fs = basins_fractions(mapper, sampler, show_progress=should_verbose())  # Almost exactly split in half (49.6% and 50.4%, cover whole domain)
# attractors = extract_attractors(mapper)

# using SplitApplyCombine

# plot(columns(attractors[1])[1:2])
# plot!(columns(attractors[2])[1:2])

# # plot(eachrow(combinedims(unobserve.(attractors[1])))[1:2])

# plot(eachcol(Matrix(attractors[1])))


# # Same analysis at the optimum
# include("07_optim_utils.jl")

# p0 = build_para(Ωmax=π/4)  # this initial guess works, others don't...
# syms = [:r, :I_eq, :Cmax, :Ωmax]

# args, vbp0, dims, (L, M, T) = make_bounds(p0, syms)
# stats, model = optimize(p0, args...);
# solution = stats.solution

# solution_sim = simulation_from_iterate(solution, vbp0, syms)
# solution_u0 = solution_sim.u[1]
# solution_tf = solution_sim.t[end]
# solution_vbp = @set vbp0[syms] = solution[end-length(syms)+1:end]

# ds2 = CoupledODEs(PM4.dynamics, solution_u0, solution_vbp; diffeq)

# function featurizer2(U, t)
#     f1 = minimum(cos.(columns(U)[1]))
#     f2 = mean(columns(U)[4])
#     f3 = mean(columns(U)[3] .^2)
#     return SA[f1, f2, f3]
# end

# # featurizer2(U, t) = SA[minimum(cos.(columns(U)[1])), mean(columns(U)[4])]
# mapper2 = AttractorsViaFeaturizing(ds2, featurizer2, T=10, Ttr=10, Δt = .01)
# sampler2(r=10) = 2r*SA[rand(5, 1)...].-r

# # Plot du cycle limite trouvé par l'opti et des autres cycles limites par sampling
# plot()
# reinit!(ds2)
# Δt1, Δt2 = 100, 20
# U, t = trajectory(ds2, Δt1)
# plot!(cos.(columns(U)[1]), cos.(columns(U)[2]), xlabel="cos(α)", ylabel="cos(τ)", c=:red)
# @progress for _ in 1:1
#     reinit!(ds2, sampler2())
#     reinit!(ds2, SA[0, π/2, 0, -1000, 0])
#     trajectory(ds2, Δt1)  # Wait some time
#     local U, t = trajectory(ds2, Δt2)
#     plot!(cos.(columns(U)[1]), cos.(columns(U)[2]), xlabel="cos(α)", ylabel="cos(τ)", label="")
# end
# display(plot!())

# # Interestingly, setting dτ = -1000 or +1000 lead to the same cycle with dτ < 0

# # Refaire à la main
# ics = [sampler2() for _ in 1:200]
# features = extract_features(mapper2, ics, show_progress=should_verbose())
# group_labels = group_features(features, mapper2.group_config)
# fs = basins_fractions(group_labels)
# scatter(invert(features)...)

# #=
# Analyse des équilibres :
# 1. plot les nullclines pour dα = dτ = 0 et en couleur le cycle limite dans lequel on tombe
#   - Combien de points d'équlibre ?
#   - Quand y a-t-il une singularité, et quelle en est la conséquence pour les cycles limites ?
#   - Quand observe-t-on 3 cycles limites, pourquoi le troisième apparait-il ? Un qui se sépart en deux suite à une variation de paramètres
#   - Le troisième cycle est-il bien atteignable depuis un équilibre ?

# 2. Autour d'un équilibre, plot les nullclines en dα et dτ, et en couleur le cycle limite dans lequel on tombe
#   - Est-ce bien linéaire, ou le paysage est-il plus compliqué ? (-> Si depuis un équilibre je donne un dτ = +ε, est-ce que je tombe toujours dans le cycle dτ > 0 ?)
# =#


# function eigenplot_2d(mat)
#     plot()
#     circ = [[cos(θ), sin(θ)] for θ in range(0, 2π, length=100)]
#     transformed_circ = [mat * circ_i for circ_i in circ]

# end

# ## Stabilité des équilibres
# import ForwardDiff: jacobian as J
# using LinearAlgebra

# ss
# q = ss[:, 3]

# # underscore to avoid redefining ddq from 05.4_functions.jl
# function ddq_(q, dq, p)
#     PM4.dynamics(SA[q..., dq..., 0], p, 0)[3:4]
# end

# function J_Δq(q, p)
#     J(Δq -> ddq_(q .+ Δq, [0, 0], p), [0, 0])
# end

# function J_dq(q, p)
#     J(dq -> ddq_(q, dq, p), [0, 0])
# end

# scatter([eigvals(J_dq(q, vbp)) for q in eachcol(ss)])
# scatter([eigvals(J_Δq(q, vbp)) for q in eachcol(ss)])
# scatter(eachrow(ss)...)


# mat = J_Δq(q, vbp)
# N = 100
# circ = [[cos(θ), sin(θ)] for θ in range(0, 2π, length=N+1)[1:end-1]]
# transformed_circ = [mat * circ_i for circ_i in circ]

# plot(ratio=:equal)
# scatter!(invert(circ)...)
# scatter!(invert(transformed_circ)...)
# quiver!(invert(circ)..., quiver=Tuple.(transformed_circ .- circ))


# m1 = mat^(1/N)
# mat_pows = [m1^i for i in range(0, N)]
# us = [m[:, 1] for m in mat_pows]
# vs = [m[:, 2] for m in mat_pows]

# scatter([0], [0], ratio=:equal)
# plot!(invert(us)..., line_z=range(0, 1, length=N))
# plot!(invert(vs)..., line_z=range(0, 1, length=N))