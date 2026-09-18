# Confirm the lower terminus is a fold (saddle-node of cycles) by marching
# v_ref downward with a damped 4-D shooting Newton, and watching the
# multipliers. At a fold, ∂R/∂s becomes singular ⇔ a nontrivial multiplier
# → +1, and Newton eventually fails (can't pass a limit point in v_ref).
#
# Shooting unknowns s = [α0, dα0, dτ0, T] (T = normalized period);
# R(s; v_ref) = endpoint_residuals(flow) ∈ R^4 (sense = +, dτ > 0).
using Pkg
Pkg.activate(@__DIR__)

using KEEP
using KEEP.PointMassPara: build_para, build_vbpara
using KEEP.PointMass4: integrate
using KEEP.LimitCycle: compute_limit_cycle, build_shooting, unpack_shooting,
    endpoint_residuals
using StaticArrays, LinearAlgebra, Printf

const NT_P0 = NamedTuple(build_para(build_vbpara()))

fd_jac(f, x; h=1e-6) = (n = length(x); J = Matrix{Float64}(undef, n, n);
    for i in 1:n
        e = zeros(n); e[i] = h
        @views J[:, i] .= (f(x .+ e) .- f(x .- e)) ./ (2h)
    end; J)

vbp_at(vr) = build_vbpara(merge(NT_P0, (v_ref=vr,)))

function resid(s, vbp)
    u0, T = unpack_shooting(SA[s[1], s[2], s[3], s[4]])
    sol = integrate(u0, T, vbp; tol=1e-11)
    return collect(endpoint_residuals(sol; sense=+))
end

function multipliers(s, vbp)
    u0, T = unpack_shooting(SA[s[1], s[2], s[3], s[4]])
    pm = v -> copy(Array(integrate(SA[v[1], v[2], v[3], v[4], 0.0], T, vbp; tol=1e-11).u[end])[1:4])
    return sort(abs.(eigvals(fd_jac(pm, collect(Array(u0)[1:4])))), rev=true)
end

function newton(s0, vbp)
    s = Vector{Float64}(s0)
    for _ in 1:60
        r = resid(s, vbp)
        nr = norm(r)
        nr < 1e-10 && return s, nr
        J = fd_jac(y -> resid(y, vbp), s)
        ds = try
            J \ (-r)
        catch
            return s, Inf
        end
        all(isfinite, ds) || return s, Inf
        λ = 1.0
        while λ > 1e-6 && norm(resid(s .+ λ .* ds, vbp)) > nr
            λ /= 2
        end
        s .+= λ .* ds
    end
    return s, norm(resid(s, vbp))
end

function main()
    # Warm-start down to v_ref = 2.5 (compute_limit_cycle tracks reliably).
    s = [-0.7771, 1.3899, 1.8811, 2.82632]
    for vr in (6.0, 4.0, 3.0, 2.75, 2.6, 2.55, 2.52, 2.5)
        s = build_shooting(compute_limit_cycle(unpack_shooting(s)[1], vbp_at(vr); tol=1e-11))
    end
    println("start (v_ref=2.5): T=", round(s[4], digits=5))

    println("\n v_ref      T        tf (s)     |mu| (sorted)                     res")
    vr = 2.5
    while vr > 2.44
        vbp = vbp_at(vr)
        s, res = newton(s, vbp)
        mu = multipliers(s, vbp)
        @printf("%7.4f  %9.5f  %8.4f  %s  %.1e%s\n", vr, s[4], s[4] * (NT_P0.l / vr),
            join([@sprintf("%.5f", m) for m in mu], "  "), res, res < 1e-8 ? "" : "  <-- FAIL")
        flush(stdout)
        res < 1e-8 || break
        vr -= 0.001
    end
end

main()
