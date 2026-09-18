# Two limit-cycle trajectories at the SAME v_ref, one on each side of the fold:
#   * "short" branch: the strongly stable cycle (tf ~ 13 s at v_ref = 2.47),
#     easily obtained by integrating / Newton from above;
#   * "long" branch: the saddle cycle past the fold (tf ~ 28 s), reachable only
#     by pseudo-arclength continuation THROUGH the limit point.
# We continue the 4-D shooting map in arclength past the fold, then refine each
# branch at v_ref = 2.47 and integrate one period.
#
# Run under the memory guard: runwatch -l 1500 -t 280 -- julia ...
using Pkg
Pkg.activate(@__DIR__)

using KEEP
using KEEP.PointMassPara: build_para, build_vbpara
using KEEP.PointMass4: integrate, compute_Rτ, compute_OK
using KEEP.LimitCycle: compute_limit_cycle, build_shooting, unpack_shooting,
    endpoint_residuals
using StaticArrays, LinearAlgebra, Printf, Logging
using Plots
Logging.disable_logging(Logging.Warn)

const NT_P0 = NamedTuple(build_para(build_vbpara()))
vbp_at(vr) = build_vbpara(merge(NT_P0, (v_ref=vr,)))
T0(vr) = NT_P0.l / vr

function R4(z, vr=z[5]; tol=1e-10)
    u0, T = unpack_shooting(SA[z[1], z[2], z[3], z[4]])
    sol = integrate(u0, T, vbp_at(vr); tol=tol)
    return collect(endpoint_residuals(sol; sense=+))
end

function tangent(z, τprev)
    J = Matrix{Float64}(undef, 4, 5)
    for i in 1:5
        e = zeros(5); e[i] = 1e-7
        @views J[:, i] .= (R4(z .+ e) .- R4(z .- e)) ./ (2e-7)
    end
    τ = svd(J).V[:, end]
    dot(τ, τprev) < 0 && (τ = -τ)
    return τ ./ norm(τ)
end

"FD multipliers of the cycle at shooting s with v_ref = vr."
function multipliers(s, vr)
    u0, T = unpack_shooting(SA[s[1], s[2], s[3], s[4]])
    vbp = vbp_at(vr); v0 = collect(Array(u0)[1:4])
    pm = v -> copy(Array(integrate(SA[v[1], v[2], v[3], v[4], 0.0], T, vbp; tol=1e-11).u[end])[1:4])
    J = Matrix{Float64}(undef, 4, 4)
    for i in 1:4
        e = zeros(4); e[i] = 1e-7
        @views J[:, i] .= (pm(v0 .+ e) .- pm(v0 .- e)) ./ (2e-7)
    end
    return sort(abs.(eigvals(J)), rev=true)
end

"Refine the cycle at fixed v_ref (4 unknowns: α0,dα0,dτ0,T)."
function solve_fixed(z0, vr)
    s = Vector{Float64}(z0[1:4]); vbp = vbp_at(vr)
    F = s -> begin
        u0, T = unpack_shooting(SA[s[1], s[2], s[3], s[4]])
        collect(endpoint_residuals(integrate(u0, T, vbp; tol=1e-11); sense=+))
    end
    for _ in 1:60
        r = F(s); norm(r) < 1e-11 && break
        J = Matrix{Float64}(undef, 4, 4)
        for i in 1:4
            e = zeros(4); e[i] = 1e-7
            @views J[:, i] .= (F(s .+ e) .- F(s .- e)) ./ (2e-7)
        end
        ds = try J \ (-r) catch; return s, Inf end
        all(isfinite, ds) || return s, Inf
        λ = 1.0
        while λ > 1e-6 && norm(F(s .+ λ .* ds)) > norm(r)
            λ /= 2
        end
        s .+= λ .* ds
    end
    return s, norm(F(s))
end

function traj(s, vr)
    u0, T = unpack_shooting(SA[s[1], s[2], s[3], s[4]])
    vbp = vbp_at(vr)
    sol = integrate(u0, T, vbp; tol=1e-11, save_everystep=true)
    t = collect(sol.t) .* T0(vr)
    U = reduce(hcat, sol.u)
    xy = [compute_OK(compute_Rτ([U[1, k], U[2, k]], vbp), vbp) for k in axes(U, 2)]
    X = reduce(hcat, xy)
    return t, U, X
end

function main()
    # seed on the stable branch at v_ref = 3
    s = [-0.7771, 1.3899, 1.8811, 2.82632]
    for vr in (6.0, 4.0, 3.0)
        s = build_shooting(compute_limit_cycle(unpack_shooting(s)[1], vbp_at(vr); tf=120, tol=1e-10))
    end
    z = [s[1], s[2], s[3], s[4], 3.0]
    τ = tangent(z, ones(5))
    τ[5] > 0 && (τ = -τ)                 # decrease v_ref first

    path = [copy(z)]
    ds = 0.05
    for step in 1:400
        ok = false
        while !ok
            zpred = z + ds .* τ
            zc = copy(zpred)
            for _ in 1:30
                G = vcat(R4(zc), [dot(τ, zc - zpred)])
                norm(G[1:4]) < 1e-7 && (ok = true; break)
                JG = Matrix{Float64}(undef, 5, 5)
                for i in 1:5
                    e = zeros(5); e[i] = 1e-7
                    @views JG[:, i] .= (vcat(R4(zc .+ e), [dot(τ, zc .+ e - zpred)]) .-
                                        vcat(R4(zc .- e), [dot(τ, zc .- e - zpred)])) ./ (2e-7)
                end
                dz = JG \ (-G); all(isfinite, dz) || break; zc .+= dz
            end
            ok && break
            ds /= 2
            ds < 5e-4 && (println("step too small at step ", step); break)
        end
        ok || (println("corrector failed at step ", step); break)
        z = zc; τ = tangent(z, τ); push!(path, z)
        ds = min(ds * 1.3, 0.1)
        @printf("step %3d  v_ref=%.5f  T=%.4f  tf=%.3f s  (ds=%.3f)\n", step, z[5], z[4], z[4] * T0(z[5]), ds)
        flush(stdout)
        z[4] * T0(z[5]) > 36 && break
    end

    vrt = 2.47
    ifold = argmin([zz[5] for zz in path])
    leg1 = 1:ifold; leg2 = ifold:length(path)
    zshort = path[leg1[argmin([abs(path[i][5] - vrt) for i in leg1])]]
    zlong = path[leg2[argmin([abs(path[i][5] - vrt) for i in leg2])]]
    sshort, rshort = solve_fixed(zshort, vrt)
    slong, rlong = solve_fixed(zlong, vrt)
    @printf("\nv_ref = %.3f\n  short: tf = %.4f s  res=%.1e  |mu| = %s\n", vrt,
        sshort[4] * T0(vrt), rshort, join([@sprintf("%.4f", m) for m in multipliers(sshort, vrt)], "  "))
    @printf("  long : tf = %.4f s  res=%.1e  |mu| = %s\n",
        slong[4] * T0(vrt), rlong, join([@sprintf("%.4f", m) for m in multipliers(slong, vrt)], "  "))

    t1, U1, X1 = traj(sshort, vrt)
    t2, U2, X2 = traj(slong, vrt)

    p = plot(layout=(2, 2), size=(900, 650))
    plot!(p[1], t1, U1[1, :]; label="short (tf=$(round(t1[end],digits=2))s)", xlabel="t (s)", ylabel="α (rad)", title="α")
    plot!(p[1], t2, U2[1, :]; label="long", lw=2)
    plot!(p[2], t1, U1[3, :]; label="short", xlabel="t (s)", ylabel="dα (rad/s)", title="dα")
    plot!(p[2], t2, U2[3, :]; label="long", lw=2)
    plot!(p[3], t1, U1[2, :]; label="short", xlabel="t (s)", ylabel="τ (rad)", title="τ")
    plot!(p[3], t2, U2[2, :]; label="long", lw=2)
    plot!(p[4], X1[1, :], X1[2, :], X1[3, :]; label="short", xlabel="x", ylabel="y", zlabel="z", title="kite path")
    plot!(p[4], X2[1, :], X2[2, :], X2[3, :]; label="long", lw=2)
    savefig(p, joinpath(@__DIR__, "BK_tests_0910_two_cycles.png"))
    println("saved BK_tests_0910_two_cycles.png")
end

main()
