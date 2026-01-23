using ADNLPModels, NLPModelsIpopt

using StaticArrays
using SplitApplyCombine: flatten

using KEEP: TAU0
using KEEP.PointMass4
using KEEP.PointMassPara
using KEEP.LimitCycle

# function flatten(x)
#     return reduce(vcat, x)
# end

# function unflatten(x)
#     # return [@view x[5i+1:5i+5] for i in 0:N]  # slower alternative
#     return eachcol(reshape(x, 5, :))
# end

#=
1. Calculer un cycle limite
2. Discrétiser le temps en N+1 points t_i = i * tf / N, i=0:N. On note h = tf/N.
3. Inconnues :
  - 1 : tf
  - 5(N+1) : α_i, τ_i, dα_i, dτ_i, W_i pour i = 0:N
  - M paramètres
3. Définir 5N contraintes avec un schéma Euler/Runge-Kutta/... (à partir d'un tableau de Butcher? des fonctions de OrdinaryDiffEq?) x_i = f(x_(i-1), h, p), i=1:N
4. Ajouter la contrainte du cycle :
    1. α_0 = α_N
    2. τ_0 = TAU0
    3. dα_0 = dα_N
    4. dτ_0 = dτ_N
    5. W_0 = 0
    6. τ_N = τ_0 - 2π
5. Maximiser W_N pour p appartenant à P. 5N+6+M inconnues avec 5N+6 contraintes et M paramètres à déterminer
=#

"""
N: Nombre d'intervalles de la discrétisation
optim_para_syms: vecteur des symboles des paramètres à optimiser
optim_para_lower, optim_para_upper: bornes des paramètres à optimiser
"""
function optimize(N, vbp, optim_para_syms, optim_para_lower, optim_para_upper)
    # Initialisation d'un état sur le cycle limite
    α0, dα0, dτ0 = 0, 0, -100
    u0 = SA[α0, TAU0, dα0, dτ0, 0]
    limit_cycle = compute_limit_cycle(u0, vbp; tol=1e-12, save_everystep=true)
    tf = limit_cycle.t[end]
    Xs = limit_cycle(0:tf/N:tf)
    optim_paras = vbp[optim_para_syms]
    dyn_pure = x -> dynamics(x, vbp, 0)

    # Construction du problème d'optimisation
    optim_vars0 = [flatten(Xs); tf; optim_paras]

    αmin, αmax = -π, π
    τmin, τmax = TAU0 - 2π, TAU0
    dαmax = 10maximum(stack(Xs)[3, :])
    dαmin = -dαmax
    dτmin, dτmax = 10minimum(stack(Xs)[4, :]), 0
    Wmin, Wmax = 0, 10maximum(stack(Xs)[5, :])

    ε = 1e-1
    lb_state = [αmin, τmin, dαmin, dτmin, Wmin] .- ε
    lb = [repeat(lb_state, outer=N+1); 0; optim_para_lower]

    ub_state = [αmax, τmax, dαmax, dτmax, Wmax] .+ ε
    ub = [repeat(ub_state, outer=N+1); 10tf; optim_para_upper]

    function obj(x) x[5N+5] / x[1] end

    function cons!(c, optim_vars)
        h = optim_vars[end] / N
        # Xs = unflatten(optim_vars[1:5N+5])  # size = 5N + 5
        Xs = eachcol(reshape(optim_vars[1:5N+5], 5, :))
        dXs = dyn_pure.(Xs)  # size = 5N + 5
        # Formule des trapèzes : x_i+1 = x_i + h/2 * (f(x_i) + f(x_i+1))
        for i in 1:N c[5i-4:5i] .= Xs[i+1] - Xs[i] - h / 2 * (dXs[i] + dXs[i+1]) end
        # residuals = Xs[2:end] - Xs[1:end-1] - h/2 * (dXs[1:end-1] + dXs[2:end])  # size = (N, 5)
        # c[1:5N] .= flatten(residuals)            # size = 5N
        c[5N+1] = Xs[1][1] - Xs[end][1]         # α_0 = α_N
        c[5N+2] = Xs[1][2] - TAU0               # τ_0 = TAU0
        c[5N+3] = Xs[1][3] - Xs[end][3]         # dα_0 = dα_N
        c[5N+4] = Xs[1][4] - Xs[end][4]         # dτ_0 = dτ_N
        c[5N+5] = Xs[1][5]                      # W_0 = 0
        c[5N+6] = Xs[end][2] - TAU0 + 2π        # τ_N = τ_0 - 2π
        return c
    end

    # StackOverflowError, see https://github.com/SciML/DiffEqCallbacks.jl/issues/29 ?
    # model = ADNLPModel!(optim_vars -> optim_vars[5N+5], optim_vars0, lb, ub, cons!, zeros(5N+6), zeros(5N+6); minimize=false)  # sparsity detection error
    model = ADNLPModel!(optim_vars -> optim_vars[5N+5], optim_vars0, lb, ub, cons!, zeros(5N+6), zeros(5N+6); minimize=false, backend=:generic)

    return stats = ipopt(model; tol=1e-3, max_wall_time=120.)  # 5 secondes
end

vbp = build_vbpara()
optim_para_syms = [:r, :I_eq]
optim_para_lower = Float64[10, 100]
optim_para_upper = Float64[100, 10000]

stats = optimize(10, vbp, optim_para_syms, optim_para_lower, optim_para_upper)