# Deflated enumeration of coexisting prograde limit cycles of the physical-time
# BVP at one wind. The Poincare sampler (all_limit_cycles) only finds ATTRACTING
# cycles; deflation penalises a known root so a Newton step can land on a saddle
# branches (or on a disconnected cycle) instead.
#
# Two anchors:
#   VREF=2.47  -> the archived short + long (saddle) pair is a valid solution;
#                 demonstrates deflation recovering the NON-attracting long cycle
#                 from the short one. Result: exactly 2 cycles.
#   VREF=9.0   -> the reference/optimization wind; one attracting cycle expected.
#                 Result: exactly 1 cycle.
#
# A perturbed-guess search for a cycle beyond the seeds is NOT attempted: a
# divergent shooting guess makes the adaptive integrator crawl and cannot be
# interrupted (with_timeout does not fire mid-integration), so it is unbounded
# in time. A full enumeration must be seeded from continuation-derived states.
#
# UNITS: BVP unknown is PHYSICAL/SI (dα,dτ rad/s, tf s). The sampler is NORMALIZED
# (T0 = l/v_ref s): dα_SI = dα_norm/T0, tf_SI = T_norm*T0.
#
# Run: VREF=2.47 MSTAR=20 julia --project=applications/sprint2026 scratch/_deflate_cycles.jl
#      VREF=9.0  MSTAR=20 julia --project=applications/sprint2026 scratch/_deflate_cycles.jl
using Serialization, Printf, LinearAlgebra, StaticArrays
include(joinpath(@__DIR__, "..", "BK_tests_0910.jl"))
import KEEP
import KEEP.PointMass4 as PM4
import KEEP.LimitCycle as LC
using KEEP.PointMassPara: build_vbpara

const OPT = (
    params_opt=(r=47.44987027243979, I_eq=3541.2653832051565, torque_slope=3862.561181135744),
    shooting=[-1.2304408672910867, 1.3582650527729334, 1.3900409853869808, 3.8474750196009815],
)
const VREF = parse(Float64, get(ENV, "VREF", "2.47"))
const MSTAR = parse(Int, get(ENV, "MSTAR", "20"))

nearest_poincare(v) = (rows = deserialize(joinpath(@__DIR__, "bvp_cycles_poincare.jls"));
    rows[argmin(abs.([r.v_ref for r in rows] .- v))])

"BVP warm start (SI): integrate the nearest normalized cycle, sample M arcs."
function warm_from_sampler(bvp, setup)
    row = nearest_poincare(VREF)
    T0 = setup.vbp.l / row.v_ref
    vbp = build_vbpara(merge(setup.nt_p0, (v_ref=row.v_ref,)))
    u0, Tn = LC.unpack_shooting(SA[row.α0, row.dα0, row.dτ0, row.T])
    sol = PM4.integrate(u0, Tn, vbp; tol=1e-11, save_everystep=true)
    tt = collect(sol.t) ./ sol.t[end]
    U = reduce(hcat, sol.u)
    interp = t -> (i = clamp(searchsortedfirst(tt, t), 1, length(tt));
        SA[U[1, i], U[2, i], U[3, i] / T0, U[4, i] / T0, Tn * T0])
    x = BVP.generate_solution(bvp, interp); x[end] = 0.0
    @printf("warm start from sampler v_ref=%.4f (tf=%.4f s)\n", row.v_ref, Tn * T0); flush(stdout)
    return x
end

"Resample an archived 5×N BVP orbit (t = PHYSICAL time 0→tf) onto the M arcs."
function bvp_from_orbit(bvp, u, t)
    tt = collect(t); tt = (tt .- tt[1]) ./ (tt[end] - tt[1])
    interp = s -> (i = clamp(searchsortedfirst(tt, s), 1, length(tt));
        SA[u[1, i], u[2, i], u[3, i], u[4, i], u[5, i]])
    x = BVP.generate_solution(bvp, interp); x[end] = 0.0
    return x
end

cycle_res(u) = maximum(abs.((u[1, end] - u[1, 1], u[2, end] - u[2, 1] - 2π,
    u[3, end] - u[3, 1], u[4, end] - u[4, 1])))

function avg_power(u, pv, setup)
    T0 = setup.vbp.l / VREF
    pw = [PM4.dynamics(SA[u[1, k], u[2, k], T0 * u[3, k], T0 * u[4, k], 0.0], get_fast_para(pv))[5]
          for k in 1:size(u, 2)]
    return sum(pw) / length(pw)
end

function main()
    setup = make_setup(opt=OPT)
    method = BVP.Shooting(MSTAR, ODE_ALG, true)
    bvp = make_bvp(make_model(method, setup), method)
    pv = merge(setup.nt_p0, (v_ref=VREF,))
    jac = make_jac(method)
    mkprob(x) = BVP.BVPBifProblem(bvp, x, pv, (@optic _.v_ref);
        jacobian=jac, record_from_solution=record_period, plot_solution=plot_solution)

    roots = Vector{Vector{Float64}}()
    function add_root!(xr, tag)
        xr === nothing && (@printf("  %-10s Newton failed\n", tag); flush(stdout); return false)
        cr = cycle_res(orbit(method, bvp, raw_x(xr), pv)[2])
        tfr = raw_x(xr)[5]
        cr > 1e-6 && (@printf("  %-10s rejected (cycle res %.1e)\n", tag, cr); flush(stdout); return false)
        any(abs(tfr - raw_x(q)[5]) < 1e-3 for q in roots) &&
            (@printf("  %-10s duplicate (tf=%.4f)\n", tag, tfr); flush(stdout); return false)
        push!(roots, collect(xr)); @printf("  %-10s root tf=%.4f s\n", tag, tfr); flush(stdout); return true
    end

    ## roots: archived pair at 2.47 (already solutions) and/or sampler-Newton.
    ## The archived orbit is resampled AND Newton-polished — node placement in
    ## the shooting mesh need not coincide with a uniform grid.
    newton(xg) = try damped_newton(mkprob(xg), xg, pv; tol=1e-9, max_iter=100)[1] catch; nothing end
    archf = joinpath(@__DIR__, "two_cycles_shooting.jls")
    if isfile(archf) && abs(VREF - 2.47) < 1e-9
        cyc = deserialize(archf)
        add_root!(newton(bvp_from_orbit(bvp, Array(cyc.short.u), cyc.short.t)), "arch short")
        add_root!(newton(bvp_from_orbit(bvp, Array(cyc.long.u), cyc.long.t)), "arch long")
    else
        add_root!(newton(warm_from_sampler(bvp, setup)), "sampler")
    end
    @assert !isempty(roots) "no root found"

    ## deflation demo: deflate root1 only, recover each other root from its state
    if length(roots) >= 2
        M1 = DeflationOperator(2.0, 1.0, [copy(roots[1])])
        for k in 2:length(roots)
            s = try BifurcationKit.solve(mkprob(roots[k]), M1,
                    NewtonPar(tol=1e-8, verbose=false, linesearch=true, max_iterations=20)) catch; nothing end
            @printf("  deflate root1 → from root%d: %s\n", k,
                (s !== nothing && BifurcationKit.converged(s)) ?
                "recovered tf=$(round(raw_x(collect(s.u))[5], digits=4)) s" : "no convergence")
            flush(stdout)
        end
    end

    ## Beyond the seeds: a perturbed-guess search with the M-arc manual-Jacobian
    ## Newton is unbounded in time here — a divergent shooting guess makes the
    ## adaptive integrator crawl, and with_timeout cannot interrupt it. The
    ## deflation DEMO above is the usable result (it recovers the non-attracting
    ## saddle); a full enumeration would seed from continuation-derived states.
    found = 0

    rows = NamedTuple[]
    for (k, xr) in enumerate(roots)
        o = orbit(method, bvp, raw_x(xr), pv)
        push!(rows, (tf=raw_x(xr)[5], α0=o[2][1, 1], mu=monodromy_max(o[2], o[1], raw_x(xr)[5], pv),
            power=avg_power(o[2], pv, setup), cres=cycle_res(o[2]), tag=k))
    end
    sort!(rows, by=r -> r.tf)
    @printf("\nv_ref=%.4f : %d distinct prograde BVP cycle(s) (%d beyond the seeds)\n", VREF, length(rows), found)
    @printf("%-4s %-10s %-10s %-10s %-11s %-10s\n", "#", "tf (s)", "α0 (rad)", "|μ|max", "power (W)", "cycle res")
    for r in rows
        @printf("%-4d %-10.4f %-+10.5f %-10.3e %-11.1f %-10.2e\n", r.tag, r.tf, r.α0, r.mu, r.power, r.cres)
    end
    serialize(joinpath(@__DIR__, "deflated_cycles_V$(VREF)_M$(MSTAR).jls"),
        (v_ref=VREF, rows=rows, n_extra=found))
end

main()
