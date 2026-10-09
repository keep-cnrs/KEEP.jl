# Landscape summary, KEEP-only: branch geometry (from the BK runs, reloaded),
# ground-truth Floquet stability, coexisting cycles, and equilibria/Hopf.
using Pkg
Pkg.activate(@__DIR__)

using KEEP
using KEEP.PointMassPara: build_para, build_vbpara, lmt
using KEEP.PointMass4: integrate, dynamics
using KEEP.LimitCycle:
    compute_limit_cycle, build_shooting, unpack_shooting, all_limit_cycles
using KEEP.SteadyState: all_steady_states_halton
using StaticArrays, LinearAlgebra, Printf, Serialization
using Plots

const NT_P0 = NamedTuple(build_para(build_vbpara()))
T0(vr) = NT_P0.l / vr

function fd_jac(f, x; h=1e-7)
    return (
        n=length(x);
        J=Matrix{Float64}(undef, n, n);
        for i in 1:n
            e = zeros(n)
            e[i] = h
            @views J[:, i] .= (f(x .+ e) .- f(x .- e)) ./ (2h)
        end;
        J
    )
end

function multipliers(s, vbp; tol=1e-11)
    u0, T = unpack_shooting(s)
    pm =
        v -> copy(
            Array(integrate(SA[v[1], v[2], v[3], v[4], 0.0], T, vbp; tol=tol).u[end])[1:4],
        )
    return eigvals(fd_jac(pm, collect(Array(u0)[1:4])))
end

function sweep(s0, vrefs)
    rows = []
    s = s0
    for vr in vrefs
        vbp = build_vbpara(merge(NT_P0, (v_ref=vr,)))
        lc = compute_limit_cycle(unpack_shooting(s)[1], vbp; tol=1e-11)
        s = build_shooting(lc)
        mu = sort(abs.(multipliers(s, vbp)); rev=true)
        push!(rows, (vr=vr, T=s[4], tf=s[4] * T0(vr), mu=mu))
    end
    return rows
end

s0 = [-0.7771, 1.3899, 1.8811, 2.82632]
rows_dn = sweep(s0, collect(9.0:-0.5:2.5))
rows_up = sweep(s0, collect(9.5:0.5:17.0))
allm = vcat(rows_dn, rows_up)

## ---- branch geometry (reloaded from the BK runs) ----
bA = deserialize(joinpath(@__DIR__, "brA_ds002.jls"))   # collocation, ds=0.02
bS = deserialize(joinpath(@__DIR__, "brS_shoot.jls"))   # shooting

pa = plot(
    bA.p,
    bA.tf;
    lw=1,
    label="collocation (brA)",
    xlabel="v_ref",
    ylabel="period tf (s)",
    title="limit-cycle branch, default params",
)
plot!(pa, bS.p, bS.tf; lw=1, alpha=0.5, label="shooting (brS)")

pb = plot(
    [r.vr for r in rows_dn],
    [r.mu[2] for r in rows_dn];
    label="from v_ref=9 down",
    yscale=:log10,
    xlabel="v_ref",
    ylabel="max nontrivial |mu|",
    title="ground-truth Floquet stability",
    ylims=(1e-4, 10),
    lw=2,
)
plot!(
    pb, [r.vr for r in rows_up], [r.mu[2] for r in rows_up]; label="from v_ref=9 up", lw=2
)
hline!(pb, [1.0]; label="|mu| = 1", color=:gray, ls=:dash)

plt = plot(pa, pb; layout=(2, 1), size=(680, 760))
savefig(plt, joinpath(@__DIR__, "BK_tests_0910_landscape.png"))
println("saved landscape figure")

## ---- coexisting cycles (Poincaré sampling) ----
println("\n=== coexisting limit cycles ===")
for vr in (9.0, 5.0, 3.0)
    vbp = build_vbpara(merge(NT_P0, (v_ref=vr,)))
    lcs = all_limit_cycles(vbp; αmin=(-π), αmax=π, vmax=12, N=40)
    for lc in lcs
        sh = build_shooting(lc)
        @printf(
            "v_ref=%5.2f  α0=%+7.4f  dα0=%+7.4f  dτ0=%+7.4f  T=%8.5f  P=%9.2f\n",
            vr,
            sh[1],
            sh[2],
            sh[3],
            sh[4],
            lc.u[end][5] / lc.t[end]
        )
    end
end

## ---- equilibria and their stability (Hopf hunt) ----
println("\n=== equilibria (α, τ) and Re(eig) of the 4x4 jacobian ===")
f4(x, p) = collect(dynamics(SA[x[1], x[2], x[3], x[4], 0.0], p)[1:4])
for vr in (2.0, 3.0, 5.0, 7.0, 9.0, 12.0, 15.0)
    vbp = build_vbpara(merge(NT_P0, (v_ref=vr,)))
    ss = all_steady_states_halton(vbp; N=60)
    @printf("v_ref=%5.2f : %d equilibria\n", vr, length(ss))
    for q in ss
        xeq = [q[1], q[2], 0.0, 0.0]
        ev = eigvals(fd_jac(x -> f4(x, vbp), xeq))
        @printf(
            "   α=%+7.4f τ=%+7.4f  Re(eig)=%s  %s\n",
            q[1],
            q[2],
            join([@sprintf("%+.4f", real(z)) for z in sort(ev, by=real, rev=true)], " "),
            maximum(real.(ev)) > 1e-6 ? "UNSTABLE" : "stable"
        )
    end
end
