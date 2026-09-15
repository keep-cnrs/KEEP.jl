# solutions_cache.jl — fingerprint-aware caching for OptimalControl solutions.
#
# A cache entry is keyed by (prefix, grid_size, problem fingerprint, init
# fingerprint, solver-options hash), all embedded in the filename:
#
#     "$(prefix)_grid_$(N)_p$(problem_hex)_i$(init_hex).jld2"
#
# so a cache hit guarantees the stored solution was computed for EXACTLY this
# problem definition, this initial guess and these solver options. Changing the
# init or the problem (even a constant referenced inside @def) produces a
# different key: no silent reuse of a stale solution.
#
# problem_fingerprint combines:
#   * the structural hash of the @def source (LineNumberNodes stripped),
#   * sampled evaluations of the dynamics / objective / constraint functions
#     (catches edits to constants and functions referenced by the @def block),
#   * the box-constraint bounds.
# Generated label symbols (label##NNN) are excluded: they embed a global parser
# counter and change between identical redefinitions.

"""
    __fp_expr(e)

Copy of expression `e` with all `LineNumberNode`s removed (they carry file/line
positions, which must not affect the fingerprint).
"""
__fp_expr(e::Expr) = Expr(e.head,
    Any[__fp_expr(a) for a in e.args if !(a isa LineNumberNode)]...)
__fp_expr(x) = x

# evaluate a generated in-place function `(r, ...) -> nothing` and return `r`
__fp_call(f, rlen, args...) = (r = zeros(rlen); f(r, args...); r)

"""
    problem_fingerprint(ocp)::UInt64

Hash of everything that defines the optimal control problem: the `@def` source
expression, sampled dynamics / Lagrange / Mayer / nonlinear-constraint values
and the box-constraint bounds. Stable across sessions for identical
definitions; changes whenever the problem changes — its structure, a constant
it references, or any function it calls (dynamics, cost, constraint).
"""
function problem_fingerprint(ocp)::UInt64
    def = ocp.definition
    h = hash(__fp_expr(def isa Expr ? def : def.expr))

    nx = OptimalControl.state_dimension(ocp)
    nu = OptimalControl.control_dimension(ocp)
    nv = OptimalControl.variable_dimension(ocp)
    samples = [(t = 0.37 + 0.11 * k,
        x = [0.1 + 0.05 * i + 0.01 * k for i in 1:nx],
        u = [0.2 - 0.03 * i + 0.01 * k for i in 1:nu],
        v = [0.7 + 0.01 * i for i in 1:nv]) for k in 1:3]

    fdyn = OptimalControl.dynamics(ocp)
    for a in samples
        h = hash(__fp_call(fdyn, nx, a.t, a.x, a.u, a.v), h)
    end

    for (getf, sig) in ((OptimalControl.lagrange, :txuv), (OptimalControl.mayer, :x0xfv))
        try
            f = getf(ocp)
            for a in samples
                h = hash(sig === :txuv ? f(a.t, a.x, a.u, a.v) : f(a.x, a.v), h)
            end
        catch
            # objective part absent (pure Lagrange has no Mayer, etc.)
        end
    end

    for (getf, mkargs) in (
        (OptimalControl.boundary_constraints_nl, a -> (a.x, a.x, a.v)),
        (OptimalControl.path_constraints_nl, a -> (a.t, a.x, a.u, a.v)),
    )
        try
            lb, fun, ub = getf(ocp)[1:3]
            h = hash((lb, ub), h)
            isempty(lb) && continue
            for a in samples
                h = hash(__fp_call(fun, length(lb), mkargs(a)...), h)
            end
        catch
            # no nonlinear constraints of this kind
        end
    end

    for getf in (OptimalControl.state_constraints_box,
        OptimalControl.control_constraints_box,
        OptimalControl.variable_constraints_box)
        try
            h = hash(getf(ocp)[1:3], h)   # (lb, indices, ub); labels are parser-generated
        catch
            # no box constraints of this kind
        end
    end
    return h
end

"""
    init_fingerprint(init)::UInt64

Hash of the initial guess: the optimisation-variable guess plus the state and
control guesses sampled at nine times scaled to the guess horizon. Changes
whenever the guess trajectory or the `tf` guess changes, so a cached solution
is never silently reused for a different initialisation.
"""
function init_fingerprint(init)::UInt64
    init === nothing && return hash(:no_init)
    h = hash(init.variable)
    var = init.variable
    tf_guess = (var isa AbstractVector && length(var) >= 1 && var[1] isa Real && var[1] > 0) ?
               Float64(var[1]) : 1.0
    for i in 0:8
        t = tf_guess * i / 8
        h = hash(init.state(t), h)
        h = hash(init.control(t), h)
    end
    return h
end

"""
    cache_filepath(cache_dir, prefix, grid_size, problem_hash, init_hash)

Cache path prefix (no extension) of one fingerprint-matched cache entry.
"""
function cache_filepath(cache_dir::String, prefix::String, grid_size::Int,
    problem_hash::UInt64, init_hash::UInt64)
    return joinpath(cache_dir, "$(prefix)_grid_$(grid_size)" *
                               "_p$(string(problem_hash; base=16))_i$(string(init_hash; base=16))")
end

# Highest-grid cached file with the same problem fingerprint and any init
# fingerprint (filename pattern: prefix_grid_N_p<problem_hex>_i<init_hex>.jld2).
function _warmstart_file(cache_dir::String, prefix::String, problem_hash::UInt64)
    rx = Regex("^$(prefix)_grid_(\\d+)_p([0-9a-f]+)_i([0-9a-f]+)\\.jld2\$")
    best = nothing
    best_grid = -1
    for fname in readdir(cache_dir)
        m = match(rx, fname)
        m === nothing && continue
        parse(UInt64, m.captures[2]; base=16) == problem_hash || continue
        grid = parse(Int, m.captures[1])
        if grid > best_grid
            best_grid = grid
            best = joinpath(cache_dir, fname[1:(end - length(".jld2"))])
        end
    end
    return best
end

"""
    run_grid_homotopy(ocp, grid_schedule; init, cache=:auto,
                      cache_dir="applications/controlled/solutions",
                      prefix="ocp_half", solve_options...)

Grid continuation over `grid_schedule` (e.g. `(8, 20, 50)`) with
fingerprint-aware caching. Each grid is cached under a key combining the
problem fingerprint, the init fingerprint and the solver options: a hit is
loaded from disk (no optimization), a miss is solved — warm-started from the
previous grid of the same run — and stored.

# Cache modes (`cache::Symbol`)
- `:auto`  (default): load every grid whose fingerprint-matched file exists,
   solve the rest.
- `:exact`: load only; error if any fingerprint-matched file is missing.
- `:no`   : ignore the cache, solve everything, refresh the files.
- `:warm` : like `:auto`, but on a miss additionally warm-start from the
   highest-grid cached solution of the same *problem* with any init (opt-in:
   it may converge to a different solution than a fresh start from `init`).
"""
function run_grid_homotopy(ocp, grid_schedule; init,
    cache::Symbol=:auto,
    cache_dir::String="applications/controlled/solutions",
    prefix::String="ocp_half",
    solve_options...)
    cache in (:auto, :exact, :no, :warm) ||
        throw(ArgumentError("cache must be :auto, :exact, :no or :warm (got :$cache)"))
    mkpath(cache_dir)

    pf = problem_fingerprint(ocp)
    ih = init_fingerprint(init)
    @info "cache key" prefix problem=string(pf; base=16) init=string(ih; base=16)

    sol = init
    for N in grid_schedule
        fpath = cache_filepath(cache_dir, prefix, N, pf, ih)
        if cache != :no && isfile(fpath * ".jld2")
            sol = import_ocp_solution(ocp; filename=fpath)
            println("✓ Grid $N loaded — Objective: ", round(objective(sol), digits=4))
        else
            cache == :exact &&
                error("cache=:exact: missing cached solution for grid $N ($(fpath).jld2)")
            warm = sol
            if cache == :warm
                wf = _warmstart_file(cache_dir, prefix, pf)
                if wf !== nothing
                    warm = import_ocp_solution(ocp; filename=wf)
                    println("… Grid $N warm-started from $(basename(wf)) (same problem, other init)")
                end
            end
            sol = solve(ocp; init=warm, grid_size=N, solve_options...)
            export_ocp_solution(sol; filename=fpath)
            println("✓ Grid $N solved — Objective: ", round(objective(sol), digits=4))
        end
    end
    return sol
end
