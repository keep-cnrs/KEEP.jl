# Three coexisting prograde limit cycles at v_ref = 2.48 (inside the window
# between the two folds), for the comparison figures.
#
# The three cycles are three crossings of ONE continuation path (short sheet ->
# fold 1 -> long sheet -> fold 2 -> long-long sheet), so a branch continuation
# that KEEPS its solutions (`ContResult.sol`) yields all three orbits. This is
# the only reliable way to get the long-long state: it is a saddle the Poincare
# sampler never reaches and there is no intermediate archive to seed Newton.
#
# Serializes scratch/three_cycles_shooting.jls in the SAME layout as
# scratch/two_cycles_shooting.jls: (; short, long, longlong, v_ref), each
#   (; x, t, u, p, tf)  with x = shooting unknown, t = PHYSICAL time (s),
#   u = 5×N BVP orbit.
#
# Run: MSTAR=40 julia --project=applications/sprint2026 scratch/_three_cycles_248.jl
using Serialization, Printf, LinearAlgebra, StaticArrays
include(joinpath(@__DIR__, "..", "BK_tests_0910.jl"))

const OPT = (
    params_opt=(r=47.44987027243979, I_eq=3541.2653832051565, torque_slope=3862.561181135744),
    shooting=[-1.2304408672910867, 1.3582650527729334, 1.3900409853869808, 3.8474750196009815],
)
const VREF = parse(Float64, get(ENV, "VREF", "2.48"))
const MSTAR = parse(Int, get(ENV, "MSTAR", "40"))
const MAXSTEPS = parse(Int, get(ENV, "MAXSTEPS", "400"))

cycle_res(u) = maximum(abs.((u[1, end] - u[1, 1], u[2, end] - u[2, 1] - 2π,
    u[3, end] - u[3, 1], u[4, end] - u[4, 1])))

function main()
    setup = make_setup(opt=OPT)
    method = BVP.Shooting(MSTAR, ODE_ALG, true)
    bvp = make_bvp(make_model(method, setup), method)
    x0 = shooting_warm_start(bvp, method, setup)
    prob = BVP.BVPBifProblem(bvp, x0, setup.nt_p0, (@optic _.v_ref);
        jacobian=make_jac(method), record_from_solution=record_period, plot_solution=plot_solution)
    x0, res_pre = damped_newton(prob, x0, setup.nt_p0; make_presolve(method)...)
    prob = BifurcationKit.re_make(prob; u0=x0)
    sol = BifurcationKit.solve(prob, Newton(), NewtonPar(tol=1e-10, verbose=false, linesearch=true))
    @assert BifurcationKit.converged(sol) "warm-start Newton did not converge"
    optn = NewtonPar(tol=1e-10, verbose=false, linesearch=true)
    optc = ContinuationPar(p_min=0.1, p_max=50.05, dsmax=0.1, ds=0.01,
        detect_bifurcation=0, newton_options=optn, max_steps=MAXSTEPS, n_inversion=6,
        save_sol_every_step=1)
    t0 = time()
    br = continuation(prob, PALC(), optc; plot=false, verbosity=0, normC=norminf, bothside=true)
    p = collect(br.branch.param); tf = collect(br.branch.tf); ss = br.sol
    @printf("branch: n=%d  sol=%d  p∈[%.4f,%.4f]  tf∈[%.3f,%.3f]  (%.1fs)\n",
        length(p), length(ss), minimum(p), maximum(p), minimum(tf), maximum(tf), time() - t0)
    @assert length(ss) == length(p) "sol/branch length mismatch: $(length(ss)) vs $(length(p))"

    pv = merge(setup.nt_p0, (v_ref=VREF,))
    jac = make_jac(method)
    mkprob(x) = BVP.BVPBifProblem(bvp, x, pv, (@optic _.v_ref);
        jacobian=jac, record_from_solution=record_period, plot_solution=plot_solution)
    bands = ((0.0, 20.0, :short), (20.0, 75.0, :long), (75.0, Inf, :longlong))
    out = Dict{Symbol,NamedTuple}()
    for (lo, hi, name) in bands
        cand = Int[]
        for i in 1:(length(p) - 1)
            (p[i] - VREF) * (p[i+1] - VREF) <= 0 || continue
            tfi = tf[i] + (VREF - p[i]) * (tf[i+1] - tf[i]) / (p[i+1] - p[i])
            (isfinite(tfi) && lo <= tfi < hi) && push!(cand, i)
        end
        if isempty(cand)
            @printf("  %-8s : NO crossing\n", name); continue
        end
        i = first(cand)
        tfi = tf[i] + (VREF - p[i]) * (tf[i+1] - tf[i]) / (p[i+1] - p[i])
        j = abs(tf[i] - tfi) <= abs(tf[i+1] - tfi) ? i : i + 1
        x = ss[j] isa NamedTuple ? ss[j].x : raw_x(ss[j])
        tfn = x[5]                                   # Newton-polish onto exactly VREF
        xp = try damped_newton(mkprob(x), x, pv; tol=1e-10, max_iter=60)[1] catch; nothing end
        xp === nothing || (x = xp)
        t, u = orbit(method, bvp, x, pv)
        U = Array(u)
        @printf("  %-8s : node tf=%.4f → polished tf=%.4f s (crossing %.4f)  cycle res=%.2e\n",
            name, tfn, x[5], tfi, cycle_res(U))
        out[name] = (; x=x, t=t .* x[5], u=u, p=VREF, tf=x[5])   # t: PHYSICAL seconds (0→tf)
        flush(stdout)
    end
    @assert length(out) == 3 "expected 3 cycles, got $(length(out))"
    serialize(joinpath(@__DIR__, "three_cycles_shooting.jls"),
        (; short=out[:short], long=out[:long], longlong=out[:longlong], v_ref=VREF, M=MSTAR))
    @printf("saved three_cycles_shooting.jls (v_ref=%.2f, M=%d)\n", VREF, MSTAR)
end

main()
