#=
The `optimize` function(s) must have a standard signature:
 - `optimize(vbp, symbols, lower, upper)`

A function `make_bounds` which, given a `p` (not VBPara), a vector `symbols` and a `multiplier` factor, computes `lower`, `upper` and `init`.
=#

module Optimization

using StaticArrays: SA
# using OrdinaryDiffEqVerner
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
    para_dims = collect((p ./ p_adim)[choosen_syms])
    state_dims = [1, 1, 1 / T, 1 / T, M * L^2 * T^-2]
    power_dim = state_dims[5] / T
    iterate_dims = [state_dims[[1, 3, 4]]..., T, power_dim, para_dims...]
    return para_dims, state_dims, iterate_dims, power_dim
end

"""
Recursively updates `target` (dimensionless) with the keys present in `guess` (dimensional).
Automatically scales the guess using `scales` to make it dimensionless.
"""
function update_initial_guess!(target, guess, scales)
    isempty(guess) && return

    for k in keys(guess)
        if k in keys(target)
            # Support Dicts, NamedTuples, and ComponentArrays natively
            val = guess isa AbstractDict ? guess[k] : getproperty(guess, k)

            # If both target and guess are nested (e.g., params), recurse
            if val isa Union{NamedTuple,AbstractDict,CA} && target[k] isa CA
                update_initial_guess!(target[k], val, scales[k])
            else
                # User provided a dimensional value, convert to dimensionless
                target[k] = val / scales[k]
            end
        else
            @warn "Initial guess key `$k` not found in optimization variables. Ignoring."
        end
    end
end

"""
Solve `model` with `solver` (`:ipopt`, or `:madnlp` when `using MadNLP`).

`hessian = :auto` keeps the solver's default (exact Hessian, differentiated through
the ODE). `:lbfgs` uses a limited-memory quasi-Newton Hessian, which costs one
constraint and one Jacobian evaluation per iteration instead of a full Hessian;
`:exact` forces the exact Hessian.
"""
function solve_nlp(model, solver::Symbol; tol, max_wall_time, linear_solver="mumps",
    hessian=:auto, print_level=0)
    solver === :ipopt && return _solve_ipopt(model; tol, max_wall_time, linear_solver,
        hessian, print_level)
    solver === :madnlp && return _solve_madnlp(model; tol, max_wall_time, hessian, print_level)
    throw(ArgumentError("unknown solver `$solver`: use :ipopt or :madnlp"))
end

function _solve_ipopt(model; tol, max_wall_time, linear_solver="mumps", hessian=:auto, print_level=0)
    kw = Dict{Symbol,Any}(:tol => tol, :max_wall_time => max_wall_time,
        :linear_solver => linear_solver, :print_level => print_level)
    if hessian === :lbfgs
        kw[:hessian_approximation] = "limited-memory"
        # Ipopt's default (6) is too short for this problem: the secant pairs are
        # skipped on curvature failures and the iteration stalls. Two pairs per
        # variable is the smallest setting that converges on the reference solve.
        kw[:limited_memory_max_history] = max(20, 2 * length(model.meta.x0))
    end
    hessian === :exact && (kw[:hessian_approximation] = "exact")
    return ipopt(model; kw...)
end

# The MadNLP backend registers itself here when its package extension loads; a plain
# method definition would be an illegal method overwrite during precompilation.
const MADNLP_BACKEND = Ref{Any}(nothing)

"Did the solver report convergence? (Ipopt returns a `Symbol`, MadNLP a SolverCore status.)"
_solved(status::Symbol) = status === :first_order
_solved(status) = occursin("SOLVE_SUCCEEDED", string(status))

function _solve_madnlp(model; kwargs...)
    MADNLP_BACKEND[] === nothing && error("solver=:madnlp requires the MadNLP extension: " *
        "run `using MadNLP` (and `import KEEP`) first.")
    return MADNLP_BACKEND[](model; kwargs...)
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
- `sense`: + for rising on the sides of the lemniscate, - for kite dropping on the sides
- `tol`: Absolute/relative tolerance for the ODE solver. The dynamics inside the
  non-linear constraints is integrated at `tol / 10`, and Ipopt is called with
  tolerance `tol`.
- `max_wall_time`: Maximum wall-clock time in seconds for the Ipopt solver.
- `linear_solver`: The linear solver to use for the Ipopt solver, cf. https://coin-or.github.io/Ipopt/OPTIONS.html#OPT_Linear_Solver. Defaults to "mumps".
- `solver`: `:ipopt` (default) or `:madnlp`; the latter needs `using MadNLP` (KEEP ships
  it as a package extension, so Ipopt-only users pay nothing).
- `hessian`: `:exact` (default) uses the exact Hessian. `:lbfgs` uses a limited-memory
  quasi-Newton Hessian — cheaper per iteration but it wanders on bad starts, so it retries
  with the exact Hessian when the cheap attempt fails or throws. `:auto` leaves the solver
  default (exact). On the whole `paper.jl` the exact Hessian is 1.9x faster (154 s against
  287 s), because failed L-BFGS attempts dominate the sweeps.
- `rhs`: `:auto` (default, the closed-form `PointMass4.dynamics`) or `:verified`
  (the AD reference `_dynamics_verified`). `:fast`/`:ad` are accepted as aliases.
- `symmetry`: exploit the invariance `S(α,τ,dα,dτ) = (-α, τ+π, -dα, dτ)` and integrate
  only half a period, imposing `x(tf) = S x(0)`. Returns the full period in `tf`.
  Faster, but fragile: the symmetric formulation has more than one stationary point and a
  seed a fraction off the symmetric manifold can settle on the wrong one (measured: 4.2 kW
  against 12.25 kW). Check the returned power before trusting it.
- `initial_guess`: A `NamedTuple`, `Dict`, or `ComponentArray` specifying initial starting values in **physical units**. This can partially or fully overwrite the baseline limit cycle guess. Because it accepts and returns physical dimensions, the returned `solution` can be passed directly into this argument for iterative refinements. Defaults to an empty NamedTuple `(;)`.

# Returns
- `solution`: A `ComponentArray` containing the optimized variables (states, time, power, and optimized parameters) mapped to their **physical dimensions**.
- `stats`: The statistics object from the Ipopt solver.
- `model`: The `ADNLPModel` object used for the optimization.

# Notes
The optimization is performed in a normalized (dimensionless) space to improve solver conditioning.
The `optim_para_lower`, `optim_para_upper`, and `initial_guess` inputs are expected in physical units and are converted to normalized units internally. While the raw `stats.solution` remains in normalized units, the returned `solution` ComponentArray is automatically scaled back to physical units.
"""
function optimize(p, optim_para_syms, optim_para_lower, optim_para_upper; sense=(+), tol=10DEFAULT_TOLERANCE, max_wall_time=120., linear_solver="mumps", solver=:ipopt, hessian=:exact,
    rhs=:auto, symmetry=false, initial_guess=(;))
    @info "To remove logging, use a Logging level below Info"

    function integrate(optim_vars::CA, vbp::VBPara)
        DType = eltype(optim_vars)
        vbp = DType.(vbp)
        (; α0, dα0, dτ0, tf, params) = optim_vars
        vbp[keys(params)] .= params
        u0 = SA[α0, TAU0, dα0, dτ0, 0]
        return (rhs === :verified || rhs === :ad) ?
               PM4.integrate_verified(u0, tf, vbp; tol=tol / 10) :
               PM4.integrate(u0, tf, vbp, tol=tol / 10)   # closed form (default)
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
        if symmetry   # half period, x(tf) = S x(0) with S = (-α, τ+π, -dα, dτ)
            c[1] = αf + α0
            c[2] = τf - TAU0 - sense(π)
            c[3] = dαf + dα0
            c[4] = dτf - dτ0
        else
            c[1] = αf - α0
            c[2] = τf - TAU0 - sense(2π)
            c[3] = dαf - dα0
            c[4] = dτf - dτ0
        end
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
    # `tf` is the integration horizon: the half period when the symmetry is used
    tf = symmetry ? limit_cycle.t[end] / 2 : limit_cycle.t[end]
    α0, _, dα0, dτ0, Wf = limit_cycle.u[end]

    params = vbp[optim_para_syms]
    # average power over the *full* cycle, whatever the integration horizon
    optim_vars = CA(α0=α0, dα0=dα0, dτ0=dτ0, tf=tf, P=Wf / limit_cycle.t[end], params=params)
    OPT_VARS_AXES = getaxes(optim_vars)

    # Use initial guess provided by user
    _, _, scales_vec, _ = compute_dims(p, optim_para_syms)
    scales = CA(scales_vec, OPT_VARS_AXES)
    update_initial_guess!(optim_vars, initial_guess, scales)
    # a user guess gives the period; convert it to the integration horizon
    symmetry && (:tf in keys(initial_guess)) && (optim_vars.tf /= 2)

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

    # Ipopt needs the starting point strictly inside the box, but a warm start taken
    # from a bound-active solution sits exactly on the bound. Compare the flat vectors:
    # `keys(optim_vars)` collapses the nested `params` field to one key, so a
    # per-symbol check cannot see a parameter sitting on its bound.
    x0, lo, hi = collect(optim_vars), collect(lb), collect(ub)
    outside = .!((lo .< x0) .& (x0 .< hi))
    if any(outside)
        @warn "optimize: initial guess clamped into the box" idxs = findall(outside) value = x0[outside] lower = lo[outside] upper = hi[outside]
        optim_vars = CA(clamp.(x0, nextfloat.(lo), prevfloat.(hi)), OPT_VARS_AXES)
    end
    @assert all(lc .<= uc) "Non-linear constraints are not feasible"

    @info "Creating model"
    model = ADNLPModel!(
        x -> obj(x, vbp), collect(optim_vars), collect(lb), collect(ub), (c, x) -> cons!(c, x, vbp), collect(lc), collect(uc); minimize=false, backend=:generic)
    @info "Starting solve"
    solve = h -> solve_nlp(model, solver; tol=tol, max_wall_time=max_wall_time,
        linear_solver=linear_solver, hessian=h, print_level=ifelse(should_verbose(), 5, 0))

    # A limited-memory Hessian is ~20x cheaper per iteration but can wander on badly
    # scaled starts (the sweep re-optimizations of `paper.jl` reach 3000 iterations and
    # non-physical parameters). Fall back to the exact Hessian when the cheap attempt
    # does not reach a feasible point, so `hessian = :lbfgs` can never return something
    # worse than what the exact Hessian would have given from the same start. The attempt
    # is also wrapped: a non-physical iterate can reach the ODE and throw (NaN horizon),
    # and that exception would escape callers such as `paper.jl`'s wind sweeps.
    stats = if hessian === :lbfgs
        fast = try
            solve(:lbfgs)
        catch err
            @warn "optimize: limited-memory Hessian attempt failed, retrying with the exact Hessian" exception = sprint(showerror, err)
            nothing
        end
        residual = fast === nothing ? Inf : norm(cons(CA(fast.solution, OPT_VARS_AXES), vbp))
        if fast !== nothing && _solved(fast.status) && residual <= max(100tol, 1e-3)
            fast
        else
            fast === nothing || @warn "optimize: limited-memory Hessian did not reach a feasible point, retrying with the exact Hessian" status = fast.status residual
            solve(:exact)
        end
    else
        solve(hessian)
    end

    ## Post-processing dimensions
    @info "Post-processing the solution into a ComponentArray"
    solution_dimless = CA(stats.solution, OPT_VARS_AXES)
    solution = solution_dimless .* scales
    symmetry && (solution.tf *= 2)   # report the full period, as without the symmetry

    return solution, stats, model
end


end  # module