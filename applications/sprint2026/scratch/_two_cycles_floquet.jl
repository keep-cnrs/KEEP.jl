# Post-process the two cycle trajectories saved by _two_cycles_bk.jl:
# ground-truth Floquet multipliers (period-map FD) + a 2x2 figure (α, dτ, τ, kite path).
using Pkg
Pkg.activate(@__DIR__)

import OrdinaryDiffEqTsit5 as ODE
using LinearAlgebra, StaticArrays, Serialization, Printf, Plots
import KEEP
using KEEP.PointMass4: dynamics, compute_Rτ, compute_OK
using KEEP.PointMassPara: build_vbpara, build_para

const NT_P0 = NamedTuple(build_para(build_vbpara()))
const PARA_CACHE = Ref{Union{Nothing, Tuple{Float64, typeof(NamedTuple(build_vbpara(build_para(build_vbpara()))))}}}(nothing)
function get_fast_para(p)
    c = PARA_CACHE[]
    if c === nothing || c[1] != p.v_ref
        c = (p.v_ref, NamedTuple(build_vbpara(p))); PARA_CACHE[] = c
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
function fd_jac(f, x; h=1e-7)
    n = length(x); J = Matrix{Float64}(undef, n, n)
    for i in 1:n
        e = zeros(n); e[i] = h
        J[:, i] = (f(x .+ e) .- f(x .- e)) ./ (2h)
    end
    return J
end

"Flow of the 4-state over a segment `ds` of the NORMALISED BVP variable
(F_fast advances the normalised variable; physical dt = tf * ds), tf fixed."
function flowseg(v, tf, ds, p)
    pr = ODE.ODEProblem((uu, pp, t) -> F_fast(uu, pp, t), vcat(v, tf), (0.0, ds), p)
    copy(Array(ODE.solve(pr, ODE.Tsit5(), abstol=1e-13, reltol=1e-13).u[end])[1:4])
end

"""
Segment-wise monodromy: product of the flow Jacobians over the SAVED orbit
segments. Robust where a single-shooting period-map FD is not — the long
(saddle) cycle is unstable enough over one period that its integrated return
misses closure by O(1) and the naive FD returns garbage.
"""
function monodromy(u, t, p)
    tf = u[5, 1]
    M = Matrix{Float64}(I, 4, 4)
    for k in 1:(size(u, 2) - 1)
        ds = (t[k+1] - t[k]) / tf
        ds <= 0 && continue
        J = fd_jac(x -> flowseg(x, tf, ds, p), u[1:4, k]; h=1e-7)
        M = J * M
        all(isfinite, M) || break
    end
    return M
end
function report(tag, d)
    p = merge(NT_P0, (v_ref=d.p,))
    mu = eigvals(monodromy(d.u, d.t, p))
    j = argmin(abs.(mu .- 1))
    nontriv = [abs(mu[k]) for k in eachindex(mu) if k != j]
    @printf("%-7s v_ref=%.4f tf=%.3fs  |mu| = %s\n", tag, d.p, d.tf,
            join([@sprintf("%.4g", abs(m)) for m in mu], ", "))
    @printf("        max non-trivial |mu| = %.4g\n", maximum(nontriv))
    return mu
end

files = [f for f in ("two_cycles_shooting.jls", "two_cycles_collocation.jls")
         if isfile(joinpath(@__DIR__, f))]
println("=== present: ", files, " ===")

plt = plot(layout=(2, 2), size=(950, 700))
colors = Dict("short" => 1, "long" => 2)
for f in files
    d = deserialize(joinpath(@__DIR__, f))
    mode = replace(f, "two_cycles_" => "", ".jls" => "")
    for tag in ("short", "long")
        haskey(d, Symbol(tag)) || continue
        dd = getfield(d, Symbol(tag))
        mu = report(tag, dd)
        p = merge(NT_P0, (v_ref=dd.p,))
        vbp = build_vbpara(p)
        X = reduce(hcat, [compute_OK(compute_Rτ([dd.u[1, k], dd.u[2, k]], vbp), vbp) for k in axes(dd.u, 2)])
        lbl = "$mode $tag (tf=$(round(dd.tf, digits=1))s)"
        plot!(plt[1], dd.t, dd.u[1, :]; label=lbl, title="α (rad)", xlabel="t (s)")
        plot!(plt[2], dd.t, dd.u[4, :]; label=lbl, title="dτ (rad/s)", xlabel="t (s)")
        plot!(plt[3], dd.t, dd.u[2, :]; label=lbl, title="τ (rad)", xlabel="t (s)")
        plot!(plt[4], X[1, :], X[2, :], X[3, :]; label=lbl, title="kite path", xlabel="x", ylabel="y", zlabel="z")
    end
end
savefig(plt, joinpath(@__DIR__, "BK_tests_0910_two_cycles.png"))
println("saved BK_tests_0910_two_cycles.png")
