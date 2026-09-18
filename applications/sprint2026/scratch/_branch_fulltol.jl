# Full-tolerance limit-cycle branch(es) for the figure, at the multiple-shooting
# arc counts given on the command line (default: 40 80). Saves
# scratch/brS_shoot_M<M>_fulltol.jls = (; p, tf, M, secs).
#
# Run: julia --project=applications/sprint2026 scratch/_branch_fulltol.jl [M...]
using Serialization, Printf
include(joinpath(@__DIR__, "..", "BK_tests_0910.jl"))

const OPT = (
    params_opt=(r=47.44987027243979, I_eq=3541.2653832051565, torque_slope=3862.561181135744),
    shooting=[-1.2304408672910867, 1.3582650527729334, 1.3900409853869808, 3.8474750196009815],
)
const MS = isempty(ARGS) ? [40, 80] : parse.(Int, ARGS)
const MAXSTEPS = parse(Int, get(ENV, "MAXSTEPS", "400"))

function branch_fulltol(M, setup)
    method = BVP.Shooting(M, ODE_ALG, true)
    bvp = make_bvp(make_model(method, setup), method)
    x0 = shooting_warm_start(bvp, method, setup)
    prob = BVP.BVPBifProblem(bvp, x0, setup.nt_p0, (@optic _.v_ref);
        jacobian=make_jac(method), record_from_solution=record_period, plot_solution=plot_solution)
    t0 = time()
    x0, res_pre = damped_newton(prob, x0, setup.nt_p0; make_presolve(method)...)
    prob = BifurcationKit.re_make(prob; u0=x0)
    sol = BifurcationKit.solve(prob, Newton(), NewtonPar(tol=1e-10, verbose=false, linesearch=true))
    @assert BifurcationKit.converged(sol) "M=$M Newton did not converge"
    optn = NewtonPar(tol=1e-10, verbose=false, linesearch=true)
    optc = ContinuationPar(p_min=0.1, p_max=50.05, dsmax=0.1, ds=0.01,
        detect_bifurcation=0, newton_options=optn, max_steps=MAXSTEPS, n_inversion=6)
    br = continuation(prob, PALC(), optc; plot=false, verbosity=0, normC=norminf, bothside=true)
    p = collect(br.branch.param); tf = collect(br.branch.tf)
    i0 = argmin(p)
    serialize(joinpath(@__DIR__, "brS_shoot_M$(M)_fulltol" *
        (MAXSTEPS == 400 ? "" : "_max$(MAXSTEPS)") * ".jls"),
        (; p=p, tf=tf, M=M, secs=time() - t0, res_pre=res_pre, maxsteps=MAXSTEPS))
    @printf("M=%3d  n=%3d  pmin=%.8f (tf=%.3f)  long-leg tfmax=%.2f  short-leg tfmin=%.4f  %.1fs  rss=%.0fMiB  maxsteps=%d\n",
        M, length(p), p[i0], tf[i0], maximum(tf[1:i0]), minimum(tf[i0:end]), time() - t0,
        Sys.maxrss() / 2^20, MAXSTEPS)
    flush(stdout)
    return p, tf
end

function main()
    setup = make_setup(opt=OPT)
    for M in MS
        branch_fulltol(M, setup)
    end
end

main()
