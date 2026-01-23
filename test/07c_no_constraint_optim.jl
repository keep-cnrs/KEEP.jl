using Test
import ForwardDiff as FD
using StaticArrays
using Setfield
using ComponentArrays
using NonlinearSolve
using Optimization, OptimizationGCMAES

using KEEP: TAU0
import KEEP.PointMass4 as PM4
using KEEP.PointMass4
using KEEP.PointMassPara
using KEEP.LimitCycle


function as_component_array(optim_vars)
    return ComponentArray(collect(optim_vars), OPT_VARS_AXES)
end


"""shoot_vars -> (u0, tf, Wf0)"""
function unpack_shoot_vars(shoot_vars)
    α0, dα0, dτ0, tf, Wf0 = shoot_vars
    return SA[α0, TAU0, dα0, dτ0, 0], tf, Wf0
end


"""(α0, dα0, dτ0, tf, Wf0) -> (αf - α0, dαf - dα0, dτf - dτ0, Wf0 - Wf)"""
function shoot(shoot_vars, vbp)
    u0, tf, Wf0 = unpack_shoot_vars(shoot_vars)
    sol = PM4.integrate(u0, tf, vbp; reltol=tol, abstol=tol, save_everystep=true)
    α0, _, dα0, dτ0, _ = u0
    αf, τf, dαf, dτf, Wf = sol.u[end]
    return SA[αf-α0, τf-TAU0+2π, dαf-dα0, dτf-dτ0, Wf0-Wf]
end

function initialize_shoot_init(vbp)
    # what did I want to do ???
end

# TODO: Test that all updating methods give the same SHOOT_INIT
"""Update SHOOT_INIT using TrustRegion on the shooting function"""
function update_shoot_init_rootfind(optim_vars, vbp)
    vbp[keys(optim_vars)] .= optim_vars
    shoot_root = solve(NonlinearProblem(shoot, SA[SHOOT_INIT...], vbp), TrustRegion(); abstol=10tol, reltol=10tol, maxiters=1000).u
    SHOOT_INIT .= shoot_root
    nothing
end
function update_shoot_init_rootfind(optim_vars::ComponentArray{<:FD.Dual}, vbp)
    update_shoot_init_rootfind(getfield.(optim_vars, :value), vbp)
end

function update_shoot_init_compute_limit_cycle(optim_vars, vbp)
    vbp[keys(optim_vars)] .= optim_vars
    u0, tf, _ = unpack_shoot_vars(SHOOT_INIT)
    sol = PM4.integrate(u0, 1000tf, vbp)
    limit_cycle = compute_limit_cycle(ode_prob; tol=tol)
    SHOOT_INIT[[1, 2, 3, 5]] .= limit_cycle.u[end][[1, 3, 4, 5]]
    SHOOT_INIT[4] = limit_cycle.t[end]
    nothing
end
function update_shoot_init_compute_limit_cycle(optim_vars::ComponentArray{<:FD.Dual}, vbp)
    update_shoot_init_compute_limit_cycle(getfield.(optim_vars, :value), vbp)
end

function update_shoot_init_integration(optim_vars, vbp)
    # # We get 2 significant digits each cycle, and we want -log10(tol) significant digits
    # nb_cycles = 2round(-log10(tol)/2, RoundUp)
    # cond_save = u -> sin((u[2] - TAU0) / 2)
    # cond_save = u -> (4*2π - abs(u[2] - TAU0)) * (3*2π - abs(u[2] - TAU0))
    # cond_stop = u -> 4.5*2π - abs(u[2] - TAU0)  # 8 cycles + a half to make sure we hit the last save at cycle n°8. save at stop is not save thanks to `save_end=false`

    # cb_save = ContinuousCallback((u, t, integrator) -> cond_save(u), identity, save_positions=(true, false))
    # # cb_save = ContinuousCallback((u, t, integrator) -> cond_save(u), int -> display("$int.t"), save_positions=(true, false))  # debug
    # # cb_stop = ContinuousCallback((u, t, integrator) -> cond_stop(u), terminate!, save_positions=(false, false))
    # cb_stop = DiscreteCallback((u, t, integrator) -> cond_stop(u) < 0, terminate!, save_positions=(false, false))
    # cb = CallbackSet(cb_save, cb_stop)

    cb1 = ContinuousCallback((u, t, integrator) -> abs(u[2] - TAU0) - 8π, identity, save_positions=(true, false))
    cb2 = ContinuousCallback((u, t, integrator) -> abs(u[2] - TAU0) - 10π, terminate!, save_positions=(true, false))
    cb = CallbackSet(cb1, cb2)

    vbp[keys(optim_vars)] .= optim_vars
    u0, tf, _ = unpack_shoot_vars(SHOOT_INIT)
    sol = PM4.integrate(u0, 100tf, vbp)
    sol = solve(ode_prob; callback=cb, save_everystep=false, save_start=false, save_end=false, reltol=tol, abstol=tol)
    yₙ = sol.u[end]
    yₙ₋₁ = sol.u[end-1]
    SHOOT_INIT[1] = rem2pi(yₙ[1], RoundNearest)  # α0
    SHOOT_INIT[[2, 3]] = yₙ[[3, 4]]  # dα0, dτ0
    SHOOT_INIT[5] = yₙ[5] - yₙ₋₁[5]  # Wf0
    SHOOT_INIT[4] = sol.t[end] - sol.t[end-1]  # tf
    nothing
end
function update_shoot_init_integration(optim_vars::ComponentArray{<:FD.Dual}, vbp)
    update_shoot_init_integration(getfield.(optim_vars, :value), vbp)
end

function obj(optim_vars::ComponentArray, vbp)
    # 1. update SHOOT_INIT with params cast to Float64
    # 2. integrate over one limit cycle using SHOOT_INIT as initial conditions (vbp can be dual if derivatives are requested)
    # 3. return the average power Wf/tf. d(Wf)/d(params) OK, what about d(tf)/d(params)?
    SHOOT_INIT[3] = -10  # ensure that the kite cycles in the right direction
    update_shoot_init_integration(optim_vars, vbp)
    vbp = convert(typeof(optim_vars), vbp)
    vbp[keys(optim_vars)] .= optim_vars
    u0, tf, _ = unpack_shoot_vars(SHOOT_INIT)
    sol = PM4.integrate(u0, 2tf, vbp)
    one_cycle_cb = ContinuousCallback((u, t, integrator) -> abs(u[2] - TAU0) - 2π, terminate!, save_positions=(true, false))
    sol = solve(ode_problem; callback=one_cycle_cb, abstol=tol, reltol=tol, save_everystep=false, save_start=false, save_end=false)
    if SHOOT_INIT[3] > 0
        @show "Le kite tourne dans le mauvais sens"
    end
    return sol.u[end][5] / sol.t[end]
end
obj(optim_vars, vbp) = obj(as_component_array(optim_vars), vbp)

function obj_one_integration(dτ_sign::Union{typeof(+),typeof(-)})
    function obj_one_integration(optim_vars::ComponentArray, vbp)
        @show 1
        vbp = convert(typeof(optim_vars), vbp)
        vbp[keys(optim_vars)] .= optim_vars

        q0 = steady_state(1, vbp; tol=10tol)
        u0 = SA[q0..., 0, dτ_sign(1), 0]
        @show 2

        cb1 = ContinuousCallback((u, t, integrator) -> u[2] - dτ_sign(8π), identity, save_positions=(true, false))
        cb2 = ContinuousCallback((u, t, integrator) -> u[2] - dτ_sign(10π), terminate!, save_positions=(true, false))
        cb = CallbackSet(cb1, cb2)
        sol = PM4.integrate(u0, 1e4, vbp)
        sol = solve(ode_prob; callback=cb, save_everystep=false, save_start=false, save_end=false, reltol=tol, abstol=tol)
        @show 3

        yₙ = sol.u[end]
        yₙ₋₁ = sol.u[end-1]
        T = sol.t[end] - sol.t[end-1]
        y_T = SA[yₙ[1], TAU0, yₙ[3], yₙ[4], yₙ[5]-yₙ₋₁[5]]
        return y_T[5] / T
    end
    obj_one_integration(optim_vars, vbp) = obj_one_integration(as_component_array(optim_vars), vbp)
    return obj_one_integration
end

#= Optimisation parameters

para: paramètres dimensionnels
α0, dα0, dτ0: condition initiale pour la recherche de cycle limite, TODO: Point d'équilibre + spécifier cycle gauche/droite
optim_para_syms: vecteur des symboles des paramètres à optimiser

=#
const SHOOT_INIT = Array{Float64,1}(undef, 5)
para1 = build_para()
para2 = build_para(C_D_l=0)
para3 = build_para(ρ_l=0)
para4 = build_para(n_wind=0)
para5 = build_para(C_D_l=0, ρ_l=0)
para6 = build_para(C_D_l=0, n_wind=0)
para7 = build_para(ρ_l=0, n_wind=0)
para8 = build_para(C_D_l=0, ρ_l=0, n_wind=0)

begin
    para = para1
    vbp0 = build_vbpara(para)
    L, M, T = lmt(vbp0)
    α0, dα0, dτ0 = 0, 0, 0
    tol = 1e-9

    # optim_para_syms = [:I_eq]
    # optim_para_lower = Float64[500]
    # optim_para_upper = Float64[2000]
    # optim_para_dims = [M * L^2]

    # optim_para_syms = [:r, :I_eq, :Δθ, :Δφ]
    # optim_para_lower = Float64[10, 100, deg2rad(5), deg2rad(5)]
    # optim_para_upper = Float64[100, 10000, deg2rad(20), deg2rad(90)]
    # optim_para_dims = [L, M * L^2, 1, 1]

    optim_para_syms = [:r, :I_eq, :Δφ]
    optim_para_lower = Float64[10, 100, deg2rad(5)]
    optim_para_upper = Float64[200, 40000, deg2rad(90)]
    optim_para_dims = [L, M * L^2, 1]

    ## Setup
    OPT_VARS_AXES = getaxes(ComponentArray(; (optim_para_syms .=> 0.)...))
    char_time = get_char_time(vbp0)
    SHOOT_INIT .= Float64[α0, dα0, dτ0, 1, 0]  # (α0, dα0, dτ0, tf, Wf0)
    vbp = normalize_vbpara(vbp0)  # without dimensions

    # Initialisation d'un état sur le cycle limite
    sol = PM4.integrate(unpack_shoot_vars(SHOOT_INIT)[1], 1000 * char_time / T, vbp)
    limit_cycle = compute_limit_cycle(ode_prob; tol=tol)
    SHOOT_INIT[[1, 2, 3]] .= limit_cycle[end][[1, 3, 4]]
    SHOOT_INIT[4] = limit_cycle.t[end]
    SHOOT_INIT[5] = limit_cycle[end][5]

    ## Contraintes en boite
    # Boites à la main
    optim_vars = vbp[optim_para_syms]
    lb = optim_para_lower ./ optim_para_dims
    ub = optim_para_upper ./ optim_para_dims

    @test all(lb .<= optim_vars .<= ub)

    opt_func = OptimizationFunction(obj, AutoForwardDiff())
    opt_prob = OptimizationProblem(opt_func, optim_vars, vbp; lb=lb, ub=ub, sense=MaxSense)
    opt_alg = OptimizationGCMAES.GCMAESOpt()
    @time sol = solve(opt_prob, opt_alg; maxtime=60)  # maxtime, maxiters
    # vérifier que l'on tourne dans le bon sens et que l'on a un cycle limite

    println([sol.u[end-length(optim_para_syms)+1:end] .* optim_para_dims; sol.objective * L * M^2 * T^-2])
    # Quelles sont les bornes saturées
    println([sol.u .* optim_para_dims; lb .* optim_para_dims; ub .* optim_para_dims])
    @test all(lb .<= sol.u .<= ub)  # Parameters inside the box
    @test SHOOT_INIT[3] < 0  # dτ0 < 0
end



sol.objective
# Quelles bornes sont saturées, quel %gain par rapport aux paramètres initiaux
sol.u .* optim_para_dims
lb .* optim_para_dims
ub .* optim_para_dims
val_before = opt_func(optim_vars, vbp) * M * L^2 * T^-2
val_after = opt_func(sol.u, vbp) * M * L^2 * T^-2
gain = (val_after - val_before) / val_before
# injecter les paramètres dans 05_LimitCycle.jl pour avoir le scaling correct de la puissance


using Plots
make_plots = false
if make_plots
    default(label="")
    opt_adim = as_component_array(sol.u)
    lb_adim = as_component_array(lb)
    ub_adim = as_component_array(ub)
    dims = as_component_array(optim_para_dims)
    opt = dims .* opt_adim
    dim_obj = M * L^2 * T^-2
    obj_val = sol.objective * dim_obj

    for para_sym in optim_para_syms
        partial = x -> obj((@set opt_adim[para_sym] = x), vbp)

        x_adim = range(lb_adim[para_sym], ub_adim[para_sym], length=51)
        y_adim = partial.(x_adim)

        x = x_adim .* dims[para_sym]
        y = y_adim .* dim_obj

        plot(x, y, label="obj($para_sym)")
        vline!(opt[[para_sym]])
        plot!(title="$(para_sym) = $(opt[para_sym])", xlabel="$(para_sym)")
        display(plot!())
    end
end



function benchmark_obj(obj, optim_vars, vbp, N)
    Dobj(optim_vars, vbp) = ForwardDiff.gradient(x -> obj(x, vbp), optim_vars)
    DDobj(optim_vars, vbp) = ForwardDiff.hessian(x -> obj(x, vbp), optim_vars)
    DDobj(optim_vars, vbp)
    return nothing

    TOL .*= 1e6
    for f in (obj, Dobj, DDobj)
        @time f(optim_vars, vbp)
    end
    TOL ./= 1e6

    for f in (obj, Dobj, DDobj)
        @time for i in 1:N
            f(optim_vars, vbp)
        end
    end
end

TOL = 1e-11
N = 1
benchmark = false
if benchmark
    benchmark_obj(obj, optim_vars, vbp, N)
end
nothing
#=
Bug avec Rodas5P -> Je ne l'ai plus... mais ça charge loooongtemps
=#
#=
0.778822 seconds (2.14 M allocations: 142.579 MiB, 3.77% gc time, 97.25% compilation time)
  2.674582 seconds (10.57 M allocations: 683.877 MiB, 10.08% gc time, 99.13% compilation time)
597.711314 seconds (17.95 M allocations: 9.102 GiB, 87.41% gc time, 1.20% compilation time)
  0.108412 seconds (2.22 k allocations: 953.656 KiB)
  0.082930 seconds (2.70 k allocations: 2.904 MiB)
515.908538 seconds (1.00 M allocations: 8.133 GiB, 88.07% gc time)
=#
