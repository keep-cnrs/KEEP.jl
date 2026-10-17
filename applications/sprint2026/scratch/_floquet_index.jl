# Floquet stability index of the coexisting prograde cycles, from their shooting
# orbits, by a reorthonormalised (QR) product of per-arc flow Jacobians.
#
# A raw product of the per-arc Jacobians is dominated by the unstable branches'
# growth (|mu| ~ 1e34 for long, 1e75 for long-long at v_ref=2.48), so the
# tangent map Phi is re-QR'd after EVERY arc and log|diag(R)| accumulated over
# `nper` periods -> Lyapunov exponents, whose signs give the number of stable
# directions (trivial exponent ~ 0). Per-arc renormalisation keeps every factor
# bounded, unlike assembling the monodromy in one block.
#
# Run: julia --project=applications/sprint2026 scratch/_floquet_index.jl
using LinearAlgebra, StaticArrays, Serialization, Printf
import OrdinaryDiffEqTsit5 as ODE
using KEEP: KEEP
using KEEP.PointMass4: dynamics
using KEEP.PointMassPara: build_vbpara, build_para

const NT_P0 = NamedTuple(build_para(build_vbpara()))
const CACHE = Ref{Union{Nothing,Tuple{Float64,Any}}}(nothing)
function fp(p)
    c = CACHE[]
    (c === nothing || c[1] != p.v_ref) &&
        (c=(p.v_ref, NamedTuple(build_vbpara(p))); CACHE[]=c)
    return c[2]
end

"Normalised BVP RHS on the 5-state (α, τ, dα, dτ, tf); tf held fixed."
function F_fast(u, params, t=0)
    α, τ, dα, dτ, tf = u
    T = params.l / params.v_ref
    _, _, ddα, ddτ, _ = dynamics(SA[α, τ, T * dα, T * dτ, 0.0], fp(params))
    out = tf .* SA[dα, dτ, ddα / T ^ 2, ddτ / T ^ 2, 0.0]
    return u isa SVector ? out : Vector(out)
end

function fd_jac(f, x; h=1e-7)
    n = length(x)
    J = Matrix{Float64}(undef, n, n)
    for i in 1:n
        e = zeros(n)
        e[i] = h
        J[:, i] = (f(x .+ e) .- f(x .- e)) ./ (2h)
    end
    return J
end

"Flow of the 4 mechanical states over a normalised arc `ds` (physical dt = tf·ds)."
function flowseg(v, tf, ds, p)
    pr = ODE.ODEProblem((uu, pp, t) -> F_fast(uu, pp, t), vcat(v, tf), (0.0, ds), p)
    return copy(Array(ODE.solve(pr, ODE.Tsit5(); abstol=1e-13, reltol=1e-13).u[end])[1:4])
end

"Per-arc flow Jacobians along the saved orbit `d` (t is PHYSICAL seconds)."
function arc_jacobians(d, pp)
    tf = d.tf
    Js = Matrix{Float64}[]
    for k in 1:(size(d.u, 2) - 1)
        ds = (d.t[k + 1] - d.t[k]) / tf
        ds <= 0 && continue
        push!(Js, fd_jac(x -> flowseg(x, tf, ds, pp), collect(d.u[1:4, k]); h=1e-7))
    end
    return Js
end

"Floquet/Lyapunov exponents of cycle `d` (per-arc QR over `nper` periods)."
function floquet_exponents(d, pp; nper=25)
    Js = arc_jacobians(d, pp)
    Φ = Matrix{Float64}(I, 4, 4)
    acc = zeros(4)
    for _ in 1:nper, J in Js
        Φ = J * Φ
        F = qr(Φ)
        acc .+= log.(abs.(diag(F.R)))
        Φ = Matrix(F.Q)
    end
    return acc ./ nper, length(Js)
end

const HERE = @__DIR__

function report_file(f, mstar)
    isfile(joinpath(HERE, f)) || return nothing
    cy = deserialize(joinpath(HERE, f))
    m = hasproperty(cy, :M) ? cy.M : mstar
    for tag in (:short, :long, :longlong)
        hasproperty(cy, tag) || continue
        d = getfield(cy, tag)
        pp = merge(NT_P0, (v_ref=d.p, tf=d.tf))
        e, narcs = floquet_exponents(d, pp)
        e = sort(e; rev=true)
        triv = e[argmin(abs.(e))]                  # trivial (≈0) direction
        # dead-zone separates the converged trivial from the genuine exponents
        stable = count(x -> x < -0.5, e)
        unstable = count(x -> x > 0.5, e)
        @printf(
            "%-32s %-9s M=%-3s tf=%8.3f arcs=%-4d  #stable=%d #unstable=%d  trivial=%+.3f  max|mu|=%.3g\n",
            f,
            tag,
            m,
            d.tf,
            narcs,
            stable,
            unstable,
            triv,
            maximum(abs, exp.(e))
        )
        @printf("    exponents=[%s]\n", join([@sprintf("%+7.3f", x) for x in e], ", "))
        flush(stdout)
    end
end

report_file("three_cycles_shooting.jls", 40)
report_file("two_cycles_shooting_M20_fulltol.jls", 20)
println("DONE")
