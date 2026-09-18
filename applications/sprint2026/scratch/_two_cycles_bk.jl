# Extract TWO limit cycles at one v_ref from a BifurcationKit BVP continuation:
# the "short" (stable) branch and the "long" (saddle) branch past the fold.
# Run for BOTH discretizations:  DISC=shooting  and  DISC=collocation.
#
#   include("_two_cycles_bk.jl")   # reads ENV["DISC"]
#
# Uses save_sol_every_step (BK ContinuationPar default = 1), so br.sol holds
# every step's orbit. Warm start comes from a KEEP cycle at v_ref = 9 (default
# params), avoiding the optimization pipeline.
using Pkg
Pkg.activate(@__DIR__)

import OrdinaryDiffEqTsit5 as ODE
using BifurcationKit, LinearAlgebra, StaticArrays, Printf, Serialization
using Plots
import KEEP
using KEEP: TAU0
using KEEP.PointMass4: dynamics, integrate
using KEEP.PointMassPara: build_vbpara, build_para
using KEEP.LimitCycle: unpack_shooting

const BVP = BifurcationKit.BVP
const MODE = get(ENV, "DISC", "shooting")
const VRT = 2.47

const VBP0 = build_vbpara()
const P0 = build_para(VBP0)
const NT_P0 = NamedTuple(P0)
const NTPARA = typeof(NamedTuple(build_vbpara(P0)))
const PARAM_CACHE = Ref{Union{Nothing, Tuple{Float64, NTPARA}}}(nothing)
function get_fast_para(params)
    c = PARAM_CACHE[]
    if c === nothing || c[1] != params.v_ref
        c = (params.v_ref, NamedTuple(build_vbpara(params)))
        PARAM_CACHE[] = c
    end
    return c[2]
end
function F_fast(u, params, t=0)
    α, τ, dα, dτ, tf = u
    T = params.l / params.v_ref
    _, _, ddα, ddτ, _ = dynamics(SA[α, τ, T*dα, T*dτ, 0], get_fast_para(params))
    out = tf .* SA[dα, dτ, ddα/T^2, ddτ/T^2, 0]
    return u isa SVector ? out : Vector(out)
end
g(u0, uT, p) = SA[uT[1]-u0[1], uT[3]-u0[3], uT[4]-u0[4], u0[2]-TAU0, uT[2]-TAU0-2π]
const STATE_SIZE = 5

const DISC = MODE == "shooting" ? BVP.Shooting(40, ODE.Tsit5(), true) :
                                  BVP.Collocation(Ntst=30, m=5, meshadapt=true)

function make_model()
    if MODE == "shooting"
        odeprob = ODE.ODEProblem(F_fast, SA[0.0, TAU0, 0.0, 0.0, 0.6], (0, 1), NT_P0)
        return BVP.BVPModel(odeprob, g; n=STATE_SIZE)
    else
        return BVP.BVPModel(F_fast, g; n=STATE_SIZE)
    end
end
model = make_model()
bvp = MODE == "shooting" ? BVP.discretize(model, DISC; abstol=1e-12, reltol=1e-10) :
                            BVP.discretize(model, DISC)

# Warm start: KEEP cycle at v_ref = 9 (default params), velocities are normalized
# in KEEP -> multiply by 1/T0 to get the physical velocities F_fast expects.
s9 = [-0.7771, 1.3899, 1.8811, 2.82632]
u0c, Tc = unpack_shooting(s9)
solc = integrate(u0c, Tc, build_vbpara(NT_P0); tol=1e-11, save_everystep=true)
tt = collect(solc.t) ./ solc.t[end]
U = reduce(hcat, solc.u)
T0_9 = 2 / 9
interp(t) = (i = clamp(searchsortedfirst(tt, t), 1, length(tt));
             SA[U[1, i], U[2, i], U[3, i] / T0_9, U[4, i] / T0_9, Tc * T0_9])
x0 = BVP.generate_solution(bvp, interp)
MODE == "shooting" && (x0[end] = Tc * T0_9)

raw_x(x) = x isa BifurcationKit.BVPSavedSolutionAndState ? BifurcationKit.saved_solution(x) : x
record_period(x, p; kw...) = (tf=raw_x(x)[5],)

"Refine a cycle at fixed `v_ref` (the saved continuation point may sit up to dsmax away)."
function refine(x; vref=VRT)
    prb = BVP.BVPBifProblem(bvp, x, merge(NT_P0, (v_ref=vref,)), (@optic _.v_ref);
        jacobian=(MODE == "shooting" ? BifurcationKit.AutoDiffDense() : BifurcationKit.DenseAnalytical()),
        record_from_solution=record_period)
    return BifurcationKit.solve(prb, Newton(), optn).u
end

# Shooting residual pin (deterministic; see BK_tests_0910.jl header).
function BVP.bvp_residual(d_bvp::BVP.DiscretizedBVP{<:BVP.BVPModel, <:BVP.Shooting}, X, p)
    model_ = BVP.get_model(d_bvp); disc_ = BVP.get_discretizer(d_bvp)
    n = BVP.state_dimension(model_); t0, tf = BVP.get_time_interval(model_); M = BVP.mesh_size(disc_)
    Xm = reshape(@view(X[1:(n*M)]), n, M)
    out = similar(X)
    outm = reshape(@view(out[1:(n*M)]), n, M)
    BVP.bvp_residual_bare!(d_bvp, outm, Xm, p, tf - t0)
    out[end] = X[end]
    return out
end

optn = NewtonPar(tol=1e-10, verbose=false, linesearch=true)
prob = BVP.BVPBifProblem(bvp, x0, NT_P0, (@optic _.v_ref);
    jacobian=(MODE == "shooting" ? BifurcationKit.AutoDiffDense() : BifurcationKit.DenseAnalytical()),
    record_from_solution=record_period)
sol = BifurcationKit.solve(prob, Newton(), optn)
println("warm-start converged: ", BifurcationKit.converged(sol))
prob = BifurcationKit.re_make(prob; u0=sol.u)

optc = ContinuationPar(p_min=0.0, p_max=25.0, dsmax=0.1, ds=0.01,
    detect_bifurcation=0, newton_options=optn, max_steps=600)
br = continuation(prob, PALC(), optc; plot=false, verbosity=1, normC=norminf, bothside=true)
println(MODE, ": continuation steps = ", length(br.branch))

# Collect saved solutions near VRT.
tfof(x) = raw_x(x)[5]
cand = [(i, s.p, tfof(s.x)) for (i, s) in enumerate(br.sol) if abs(s.p - VRT) < 0.12]
sort!(cand, by=last)
if isempty(cand) || maximum(last, cand) / minimum(last, cand) < 1.3
    println(MODE, ": only ONE branch near v_ref=", VRT, " (did not pass the fold): ",
            isempty(cand) ? "no saved point" :
            "tf ≈ " * string(round(minimum(last, cand), digits=2)) * " s")
else
    xshort = refine(raw_x(br.sol[cand[1][1]].x))
    xlong = refine(raw_x(br.sol[cand[end][1]].x))
    tf1, tf2 = tfof(xshort), tfof(xlong)
    println(MODE, ": short/long at v_ref=", VRT, " -> tf = ",
            round(tf1, digits=3), " s / ", round(tf2, digits=3), " s")
    n = BVP.state_dimension(bvp)
    function orbit_of(x, p, tf)
        s = MODE == "shooting" ? BVP.get_solution_bvp(bvp, @view(x[1:(n * DISC.M)]), BifurcationKit.setparam(prob, p)) :
                                 BVP.get_solution_bvp(bvp, x, BifurcationKit.setparam(prob, p))
        return collect(s.t) ./ s.t[end] .* tf, s.u
    end
    pv = merge(NT_P0, (v_ref=VRT,))
    (t1, u1) = orbit_of(xshort, pv, tf1); (t2, u2) = orbit_of(xlong, pv, tf2)
    p = plot(layout=(2, 1), size=(850, 650))
    plot!(p[1], t1, u1[1, :]; label="short tf=$(round(tf1, digits=2))s", xlabel="t (s)", ylabel="α (rad)", title="$MODE: α(t)")
    plot!(p[1], t2, u2[1, :]; label="long  tf=$(round(tf2, digits=2))s", lw=2)
    plot!(p[2], u1[1, :], u1[3, :]; label="short", xlabel="α (rad)", ylabel="dα (rad/s)", title="phase portrait")
    plot!(p[2], u2[1, :], u2[3, :]; label="long", lw=2)
    savefig(p, joinpath(@__DIR__, "BK_tests_0910_two_cycles_$(MODE).png"))
    serialize(joinpath(@__DIR__, "two_cycles_$(MODE).jls"),
        (short=(t=t1, u=u1, p=VRT, tf=tf1), long=(t=t2, u=u2, p=VRT, tf=tf2)))
    println(MODE, ": saved BK_tests_0910_two_cycles_$(MODE).png")
end
