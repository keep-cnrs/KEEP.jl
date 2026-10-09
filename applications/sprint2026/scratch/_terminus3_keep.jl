# Pseudo-arclength continuation of the shooting map to pass the lower limit
# point in v_ref, and watch a nontrivial Floquet multiplier -> +1 (fold =
# saddle-node of cycles) while the period stays finite (not homoclinic).
#
# z = [α0, dα0, dτ0, T, v_ref]  (5 unknowns)
# R(z)[1:4] = endpoint residuals ; augmented constraint  τ·(z - z_pred) = 0.
using Pkg
Pkg.activate(@__DIR__)

using KEEP
using KEEP.PointMassPara: build_para, build_vbpara
using KEEP.PointMass4: integrate
using KEEP.LimitCycle:
    compute_limit_cycle, build_shooting, unpack_shooting, endpoint_residuals
using StaticArrays, LinearAlgebra, Printf, Logging
Logging.disable_logging(Logging.Warn)   # keep the ODE divergence warnings out of the log

const NT_P0 = NamedTuple(build_para(build_vbpara()))
vbp_at(vr) = build_vbpara(merge(NT_P0, (v_ref=vr,)))

function R4(z)
    u0, T = unpack_shooting(SA[z[1], z[2], z[3], z[4]])
    sol = integrate(u0, T, vbp_at(z[5]); tol=1e-11)
    return collect(endpoint_residuals(sol; sense=(+)))
end

function fd_jac!(J, f, x, fx; h=1e-7)
    n = length(x)
    e = zeros(n)
    for i in 1:n
        fill!(e, 0.0)
        e[i] = h
        @views J[:, i] .= (f(x .+ e) .- f(x .- e)) ./ (2h)
    end
    return J
end

"Right singular vector of the 4×5 Jacobian for the smallest singular value."
function tangent(z, τprev)
    z0 = copy(z)
    r0 = R4(z)
    J = Matrix{Float64}(undef, 4, 5)
    fd_jac!(J, R4, z0, r0; h=1e-7)
    τ = svd(J).V[:, end]
    dot(τ, τprev) < 0 && (τ = -τ)
    return τ ./ norm(τ)
end

function multipliers(z)
    u0, T = unpack_shooting(SA[z[1], z[2], z[3], z[4]])
    vbp = vbp_at(z[5])
    v0 = collect(Array(u0)[1:4])
    pm =
        v -> copy(
            Array(integrate(SA[v[1], v[2], v[3], v[4], 0.0], T, vbp; tol=1e-11).u[end])[1:4],
        )
    J = fd_jac!(Matrix{Float64}(undef, 4, 4), pm, v0, pm(v0))
    return sort(abs.(eigvals(J)); rev=true)
end

function main()
    s = [-0.7771, 1.3899, 1.8811, 2.82632]
    for vr in (6.0, 4.0, 3.0, 2.75, 2.6, 2.55, 2.52, 2.5)
        s = build_shooting(
            compute_limit_cycle(unpack_shooting(s)[1], vbp_at(vr); tol=1e-11)
        )
    end
    z = [s[1], s[2], s[3], s[4], 2.5]
    τ = tangent(z, ones(5))
    # orient so v_ref initially decreases
    τ[5] > 0 && (τ = -τ)

    println("\n  v_ref      T        tf (s)     |mu| (sorted)                    res")
    ds = 0.02
    for step in 1:400
        zpred = z + ds .* τ
        zc = copy(zpred)
        ok = false
        for _ in 1:25
            G = vcat(R4(zc), [dot(τ, zc - zpred)])
            norm(G[1:4]) < 1e-10 && (ok=true; break)
            JG = fd_jac!(
                Matrix{Float64}(undef, 5, 5),
                w -> vcat(R4(w), [dot(τ, w - zpred)]),
                zc,
                G;
                h=1e-7,
            )
            dz = try
                JG \ (-G)
            catch
                break
            end
            all(isfinite, dz) || break
            zc .+= dz
        end
        res = norm(R4(zc))
        mu = multipliers(zc)
        @printf(
            "%8.4f  %9.5f  %8.4f  %s  %.1e%s\n",
            zc[5],
            zc[4],
            zc[4] * (NT_P0.l / zc[5]),
            join([@sprintf("%.5f", m) for m in mu], "  "),
            res,
            ok ? "" : "  <-- no conv"
        )
        flush(stdout)
        ok || break
        z = zc
        τ = tangent(z, τ)
    end
end

main()
