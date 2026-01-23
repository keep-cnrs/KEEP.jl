#=
The `optimize` function(s) must have a standard signature:
 - `optimize(vbp, symbols, lower, upper)`

A function `make_bounds` which, given a `p` (not VBPara), a vector `symbols` and a `multiplier` factor, computes `lower`, `upper` and `init`.
=#

module Optimization

using StaticArrays: SA
# using OrdinaryDiffEqTsit5
using ADNLPModels, NLPModelsIpopt
using Setfield: @set
using Logging
using LinearAlgebra: norm
import ComponentArrays: ComponentArray as CA, getaxes

using KEEP.PointMassPara
import KEEP: DEFAULT_TOLERANCE, TAU0
import KEEP.PointMass4 as PM4
using KEEP.LimitCycle: compute_limit_cycle

export make_bounds, compute_dims, optimize


const MIN_TF = 1e-3
const MAX_TF = 1e2
const MIN_POWER = 1e-3
const MAX_POWER = 1e5
const MIN_ALPHA = -π
const MAX_ALPHA = π
const MIN_D_ALPHA = -100.0
const MAX_D_ALPHA = 100.0
const MIN_ABS_D_TAU = 0.0
const MAX_ABS_D_TAU = 100.0
const CONSTRAINT_TOL_MULTIPLIER = 1.0

"""
Check verbosity of Logging
"""
function should_verbose()
  return Logging.min_enabled_level(current_logger()) <= Info
end

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
  para_dims = collect((p./p_adim)[choosen_syms])
  state_dims = [1, 1, 1 / T, 1 / T, M * L^2 * T^-2]
  power_dim = state_dims[5] / T
  iterate_dims = [state_dims[[1, 3, 4]]..., T, power_dim, para_dims...]
  return para_dims, state_dims, iterate_dims, power_dim
end

"""
Performs optimization to find optimal limit cycle parameters using Ipopt.

# Arguments
- `p`: An initial `VBPara` (physical units) used to establish base parameters
  and characteristic dimensions.
- `optim_para_syms`: A vector of `Symbol`s for the `VBPara` parameters that will
  be optimized.
- `optim_para_lower`: A vector of `Float64` specifying the lower bounds for the
  `optim_para_syms` *in their physical units*.
- `optim_para_upper`: A vector of `Float64` specifying the upper bounds for the
  `optim_para_syms` *in their physical units*.
- `tol`: Absolute/relative tolerance for the ODE solver. The optimization
  solver's tolerance will be `sqrt(tol)`.
- `max_wall_time`: Maximum wall-clock time in seconds for the Ipopt solver.
- `linear_solver`: The linear solver to use for the Ipopt solver, cf. https://coin-or.github.io/Ipopt/OPTIONS.html#OPT_Linear_Solver. Defaults to "mumps".

# Returns
- `stats`: The statistics object from the Ipopt solver.
- `model`: The `ADNLPModel` object used for the optimization.

# Notes
The optimization is performed in a normalized (dimensionless) space.
The `optim_para_lower` and `optim_para_upper` inputs are in physical units,
but they are converted to normalized units internally for the optimization.
The `stats.solution` will be in normalized units.
"""
function optimize(p, optim_para_syms, optim_para_lower, optim_para_upper; sense=+, tol=10DEFAULT_TOLERANCE, max_wall_time=120., linear_solver="mumps")

  function integrate(optim_vars::CA, vbp::VBPara)
    DType = eltype(optim_vars)
    vbp = DType.(vbp)
    (; α0, dα0, dτ0, tf, optim_params) = optim_vars
    vbp[keys(optim_params)] .= optim_params
    u0 = SA[α0, TAU0, dα0, dτ0, 0]
    return PM4.integrate(u0, tf, vbp, tol=tol / 10)
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
    c[2] = τf - TAU0 - sense(2π)
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
  limit_cycle = compute_limit_cycle(vbp; sense=sense, tol=tol / 10)
  tf = limit_cycle.t[end]
  α0, _, dα0, dτ0, Wf = limit_cycle.u[end]

  optim_paras = vbp[optim_para_syms]
  optim_vars = CA(α0=α0, dα0=dα0, dτ0=dτ0, tf=tf, P=Wf / tf, optim_params=optim_paras)
  OPT_VARS_AXES = getaxes(optim_vars)

  @info "Checking that constraints are satisfied for the initial guess"
  guess_error = norm(cons(optim_vars, vbp))
  if guess_error > tol
    @warn "optimize: Initial guess breaks non-linear constraints : residuals = $guess_error > tol = $tol"
  end

  lower_vbp = build_vbpara(@set p[optim_para_syms] = optim_para_lower)
  upper_vbp = build_vbpara(@set p[optim_para_syms] = optim_para_upper)

  ## Contraintes en boite
  # Boite arbitraire pour α, dα, dτ, tf, P
  min_d_tau, max_d_tau = minmax(sense(MIN_ABS_D_TAU), sense(MAX_ABS_D_TAU))
  lb = CA([MIN_ALPHA, MIN_D_ALPHA, min_d_tau, MIN_TF, MIN_POWER, lower_vbp[optim_para_syms]...], OPT_VARS_AXES)
  ub = CA([MAX_ALPHA, MAX_D_ALPHA, max_d_tau, MAX_TF, MAX_POWER, upper_vbp[optim_para_syms]...], OPT_VARS_AXES)

  # Boite à ±10% de l'Initialisation
  # ε = .1
  # lb = min.((1 - ε) * optim_vars, 1 / (1 - ε) * optim_vars)
  # ub = max.((1 + ε) * optim_vars, 1 / (1 + ε) * optim_vars)

  ## Contraintes de raccordement en bout de cycle
  lc = uc = zero(cons(optim_vars, vbp))

  symbols_not_in_bounds = [s for s in keys(optim_vars) if !(lb[s] < optim_vars[s] < ub[s])]
  if length(symbols_not_in_bounds) > 0
    s = first(symbols_not_in_bounds)
    @assert false "Optimization variable $s with value $(optim_vars[s]) is not in bounds [$(lb[s]), $(ub[s])]."
  end
  @assert all(lc .<= uc) "Non-linear constraints are not feasible"

  @info "Creating model"
  model = ADNLPModel!(
    x -> obj(x, vbp), collect(optim_vars), collect(lb), collect(ub), (c, x) -> cons!(c, x, vbp), collect(lc), collect(uc); minimize=false, backend=:generic)
  @info "Starting solve"
  stats = ipopt(model; tol=tol, max_wall_time=max_wall_time, linear_solver=linear_solver, print_level=ifelse(should_verbose(), 5, 0))
  @info "To remove logging, use a Logging level below Info"

  # sim = integrate(stats.solution, vbp)
  # @show sim.prob.p
  # @show getaxes(sim.prob.p)
  return stats, model
end


end  # module