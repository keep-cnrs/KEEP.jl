using StaticArrays: SA
using OrdinaryDiffEqTsit5
using LinearAlgebra: norm
using ADNLPModels, NLPModelsIpopt
using Setfield: @set
using Logging
import ComponentArrays: ComponentArray as CA, getaxes

using KEEP.PointMassPara
import KEEP: DEFAULT_TOLERANCE
import KEEP.LimitCycle: compute_limit_cycle
using KEEP.LimitCycle: compute_limit_cycle

include("00_general_functions.jl")

const MIN_TF = 1e-3
const MAX_TF = 1e2
const MIN_POWER = 1e-3
const MAX_POWER = 1e5
const MIN_ALPHA = -π
const MAX_ALPHA = π
const MIN_D_ALPHA = -100.0
const MAX_D_ALPHA = 100.0
const MIN_D_TAU = -100.0
const MAX_D_TAU = 0.0  # Assuming dτ0 is negative.
const CONSTRAINT_TOL_MULTIPLIER = 1.0

#=
p
p -> vbp
## make_bounds
p -> [center, lower, upper]
## optimize
=#

"""
Automatically builds the lower and upper bounds
"""
function make_bounds(p, choosen_syms; multiplier=3)
    center = collect(p[choosen_syms])
    lower = center ./ multiplier
    upper = center .* multiplier
    return lower, upper
end

function compute_dims(p, choosen_syms)
    vbp0 = build_vbpara(p)
    L, M, T = lmt(vbp0)
    p_adim = build_para(normalize_vbpara(vbp0))
    para_dims = collect((p ./ p_adim)[choosen_syms])
    state_dims = [1, 1, 1/T, 1/T, M * L^2 * T^-2]
    iterate_dims = [state_dims[[1, 3, 4]]..., T, state_dims[5], para_dims...]
    return para_dims, state_dims, iterate_dims
end

"""
tol: tolerance (abs/rel) of the ODE solver, taking sqrt(tol) for the optimization problem
"""
function optimize(p , optim_para_syms, optim_para_lower, optim_para_upper; tol=DEFAULT_TOLERANCE, max_wall_time=120.)
    function integrate(optim_vars::CA, vbp::VBPara)
        DType = eltype(optim_vars)
        vbp = DType.(vbp)
        (; α0, dα0, dτ0, tf, optim_params) = optim_vars
        vbp[keys(optim_params)] .= optim_params
        u0 = SA[α0, TAU0, dα0, dτ0, 0]
        return PM4.integrate(u0, tf, vbp, tol=tol)
    end

    # convert optim_vars to CA
    integrate(optim_vars, vbp::VBPara) = integrate(CA(optim_vars, OPT_VARS_AXES), vbp)

    # Trivial objective function
    obj(optim_vars::CA, vbp) = optim_vars.P
    obj(optim_vars, vbp) = obj(CA(optim_vars, OPT_VARS_AXES), vbp)

    function cons!(c, optim_vars, vbp)
        ode_sol = integrate(optim_vars, vbp)
        α0, dα0, dτ0, tf, P, _ = optim_vars
        αf, τf, dαf, dτf, Wf = ode_sol.u[end]
        c[1] = αf - α0
        c[2] = τf - TAU0 + 2π
        c[3] = dαf - dα0
        c[4] = dτf - dτ0
        c[5] = 2(Wf - P * tf) / (Wf + P * tf + eps(Wf))
        return c
    end
    cons(optim_vars::AbstractArray{T}, vbp) where {T} = cons!(Array{T}(undef, 5), optim_vars, vbp)  # wrapper used O(1) times

    vbp0 = build_vbpara(p)

    ## Setup
    # Variable d'optimisation : α0, dα0, dτ0, tf, P, [paramètres]...
    vbp = normalize_vbpara(vbp0)  # without dimensions
    
    # Initialisation d'un état sur le cycle limite
    @info "Computing limit cycle for the initial guess"
    limit_cycle = compute_limit_cycle(vbp; sense=-, tol=tol)
    tf = limit_cycle.t[end]
    α0, _, dα0, dτ0, Wf = limit_cycle.u[end]
    
    optim_paras = vbp[optim_para_syms]
    optim_vars = CA(α0=α0, dα0=dα0, dτ0=dτ0, tf=tf, P=Wf / tf, optim_params=optim_paras)
    OPT_VARS_AXES = getaxes(optim_vars)

    @info "Checking that constraints are satisfied for the initial guess"
    @assert all(abs.(cons(optim_vars, vbp)) .< sqrt(tol)) "Initial guess breaks non-linear constraints"

    lower_vbp = build_vbpara(@set p[optim_para_syms] = optim_para_lower)
    upper_vbp = build_vbpara(@set p[optim_para_syms] = optim_para_upper)

    ## Contraintes en boite
    # Boite arbitraire pour α, dα, dτ, tf, P
    lb = CA([MIN_ALPHA, MIN_D_ALPHA, MIN_D_TAU, MIN_TF, MIN_POWER, lower_vbp[optim_para_syms]...], OPT_VARS_AXES)
    ub = CA([MAX_ALPHA, MAX_D_ALPHA, MAX_D_TAU, MAX_TF, MAX_POWER, upper_vbp[optim_para_syms]...], OPT_VARS_AXES)

    # Boite à ±10% de l'Initialisation
    # ε = .1
    # lb = min.((1 - ε) * optim_vars, 1 / (1 - ε) * optim_vars)
    # ub = max.((1 + ε) * optim_vars, 1 / (1 + ε) * optim_vars)

    ## Contraintes de raccordement en bout de cycle
    lc = uc = zero(cons(optim_vars, vbp))

    @assert all(lb .<= optim_vars .<= ub) "Box does not contain initial guess"
    @assert all(lc .<= uc) "Non-linear constraints are not feasible"

    @info "Creating model"
    model = ADNLPModel!(
        x -> obj(x, vbp), collect(optim_vars), collect(lb), collect(ub), (c, x) -> cons!(c, x, vbp), collect(lc), collect(uc); minimize=false, backend=:generic)
    @info "Starting solve"
    stats = ipopt(model; tol=sqrt(tol), max_wall_time=max_wall_time, print_level=ifelse(should_verbose(), 5, 0))

    # sim = integrate(stats.solution, vbp)
    # @show sim.prob.p
    # @show getaxes(sim.prob.p)
    return stats, model
end


"""iterate : [α0, dα0, dτ0, tf, avg_P, params...]
vbp0 : reference vbpara
syms : ComponentArrays.getaxes(params)"""
function simulation_from_iterate(iterate, vbp0, syms; tol=DEFAULT_TOLERANCE)
    α0, dα0, dτ0, tf, _, params... = iterate
    vbp = @set vbp0[syms] = params
    u0 = SA[α0, TAU0, dα0, dτ0, 0.0]
    return PM4.integrate(u0, tf, vbp; save_everystep=true, tol=tol)
end

"""build an iterate from a simulation"""
function iterate_from_simulation(sim, syms)
    vbp = sim.prob.p
    L, M, T = lmt(vbp)
    α0, dα0, dτ0 = sim.u[1][[1, 3, 4]] ./ SA[1, 1/T, 1/T]
    Wf = sim.u[end][end] / (M * L^2 * T^-2)
    tf = sim.t[end] / T
    return [α0, dα0, dτ0, tf, Wf / tf, vbp[syms]...]
end

#=

function make_bounds(p, choosen_syms)
    vbp0 = build_vbpara(p)
    L, M, T = lmt(vbp0)

    # [:r, :I_eq, :Δθ, :Δφ, :Ωmin, :Ωmax, :Cmax]
    DIMS = [L, M * L^2, 1, 1, T^-1, T^-1, M*L^2*T^-2]

    inds = findall(s -> s ∈ choosen_syms, SYMS)
    # inds = findfirst.(isequal.(choosen_syms), Ref(SYMS))
    optim_para_syms = SYMS[inds]
    optim_para_lower = LOWER[inds]
    optim_para_upper = UPPER[inds]
    optim_para_dims = DIMS[inds]
    return (optim_para_syms, optim_para_lower, optim_para_upper), vbp0, optim_para_dims, (L, M, T)
end

=#
