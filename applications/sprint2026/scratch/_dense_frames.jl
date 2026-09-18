# Dense frame states for anim_three_cycles.jl, produced the BifurcationKit-native
# way: the shooting unknown x is re-solved through BK's own ensemble flow
# (`BifurcationKit.evolve(sh.flow, Val(:Full), um, p, sh.ds)` -> one dense
# ODESolution per shooting arc) and sampled at the animation frame phases.
#
# Why per-arc: a single whole-period solve from x diverges (closure drift
# 8 / 1e67 / 1e109 for short/long/long-long — the saddles amplify local error by
# e^{λT}); the M-arc decomposition is what BK itself stores. Why not
# BVP.get_po_solution: it reads the period as getperiod(sh,x)=x[end]=tf
# (physical) while the BVP flow is normalized to tspan (0,1), so its arcs are
# tf× too long and its POInterpolation walks an EnsembleSolution as a matrix.
#
# Serializes scratch/three_cycles_frames.jls:
#   (; short=(;frames, trail, tf), long=..., longlong=..., v_ref, SNAPS, M)
# frames = 2×SNAPS (α, τ) at s = (i-1)/SNAPS (normalized, [0,1));
# trail  = Vector of (α,τ) mesh points for the faint background.
#
# Run: julia --project=applications/sprint2026 scratch/_dense_frames.jl
using Serialization, Printf, LinearAlgebra, StaticArrays
import BifurcationKit
include(joinpath(@__DIR__, "..", "BK_tests_0910.jl"))

const OPT = (
    params_opt=(r=47.44987027243979, I_eq=3541.2653832051565, torque_slope=3862.561181135744),
    shooting=[-1.2304408672910867, 1.3582650527729334, 1.3900409853869808, 3.8474750196009815],
)
const SNAPS = 300

function dense_frames(c, bvp, sh, p)
    N = BVP.state_dimension(bvp)
    M = length(sh.ds)
    um = reshape(c.x[1:(N * M)], N, M)                     # M shooting states
    ens = BifurcationKit.evolve(sh.flow, Val(:Full), um, p, sh.ds)  # M dense arcs
    t0 = cumsum(vcat(0.0, sh.ds))
    F = Matrix{Float64}(undef, 2, SNAPS)
    for i in 1:SNAPS
        s = (i - 1) / SNAPS
        k = clamp(searchsortedlast(t0, s), 1, M)
        q = ens.u[k](s - t0[k])
        F[1, i] = q[1]
        F[2, i] = q[2]
    end
    # x closes under the rebuilt BVP (validates it is the periodic orbit)
    xend = ens.u[M](sh.ds[M]); xs = um[:, 1]
    res = maximum(abs.((xend[1] - xs[1], xend[2] - xs[2] - 2π, xend[3] - xs[3], xend[4] - xs[4])))
    # dense frames agree with the archived mesh (Hermite, same phase)
    mesh = BVP.get_solution_bvp(bvp, c.x[1:(N * M)], p)
    U = Array(mesh.u); sn = collect(mesh.t)
    err = 0.0
    for i in (1, 75, 150, 225, 299)
        s = (i - 1) / SNAPS
        j = clamp(searchsortedlast(sn, s), 1, length(sn) - 1)
        w = (s - sn[j]) / (sn[j + 1] - sn[j]); h = (sn[j + 1] - sn[j]) * c.tf
        h00 = 2w^3 - 3w^2 + 1; h10 = w^3 - 2w^2 + w; h01 = -2w^3 + 3w^2; h11 = w^3 - w^2
        her = h00 .* U[1:2, j] .+ h10 .* h .* U[3:4, j] .+ h01 .* U[1:2, j + 1] .+ h11 .* h .* U[3:4, j + 1]
        err = max(err, norm(F[:, i] .- her))
    end
    @printf("closure res=%.2e   |frames-mesh|=%.2e\n", res, err)
    @assert res < 1e-6 "x does not close: res=$res"
    @assert err < 1e-5 "dense frames inconsistent with the mesh: err=$err"
    (; frames=F, trail=[(U[1, k], U[2, k]) for k in 1:2:size(U, 2)], tf=c.tf)
end

function main()
    cy = deserialize(joinpath(@__DIR__, "three_cycles_shooting.jls"))
    setup = make_setup(opt=OPT)
    method = BVP.Shooting(cy.M, ODE_ALG, true)
    bvp = make_bvp(make_model(method, setup), method)
    sh = bvp.cache
    p = merge(setup.nt_p0, (v_ref=cy.v_ref,))
    out = Dict{Symbol,NamedTuple}()
    for (name, c) in ((:short, cy.short), (:long, cy.long), (:longlong, cy.longlong))
        @printf("  %-8s ", name)
        out[name] = dense_frames(c, bvp, sh, p)
    end
    serialize(joinpath(@__DIR__, "three_cycles_frames.jls"),
        (; short=out[:short], long=out[:long], longlong=out[:longlong],
            v_ref=cy.v_ref, SNAPS=SNAPS, M=cy.M))
    @printf("saved three_cycles_frames.jls (SNAPS=%d, M=%d)\n", SNAPS, cy.M)
end

main()
