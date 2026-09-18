# _fold_tangle_probe.jl  —  TEMPORARY exploration (not part of the deliverable)
#
# QUESTION
#   `BK_tests_0910.jl` / REPORT.md report a "fold tangle" near v_ref ≈ 9.2–9.8 m/s
#   (several :fold special points crowding sheets, plus corrector sheet-jumps).
#   Is anything genuinely new there — real saddle-nodes / extra periodic orbits —
#   or is it a continuation artifact (BK's `detect_fold` firing on a nearly flat
#   tf(v_ref) plus the corrector hopping between symmetry-related cycles)?
#
# TESTS (all independent of the BK continuation that produced the tangle)
#   A. Ground-truth Floquet multipliers of the tracked stable cycle over
#      [8.8, 10.2]. A genuine fold of THIS sheet requires a nontrivial |mu| → 1.
#   B. Orbit census by KEEP Poincaré sampling (`all_limit_cycles`) in the band:
#      how many distinct cycles exist, and are they just the α ↔ −α mirror pair?
#   C. The archived BK branch trace (`scratch/brA_ds002.jls`): does the recorded
#      parameter actually turn around (fold) in the band, or just adapt ds?
#
# Report appended at the bottom once the exploration finished.

using Pkg
Pkg.activate(@__DIR__)

using KEEP
using KEEP.PointMassPara: build_para, build_vbpara
using KEEP.PointMass4: integrate
using KEEP.LimitCycle: compute_limit_cycle, build_shooting, unpack_shooting,
    all_limit_cycles
using StaticArrays, LinearAlgebra, Printf, Serialization

const NT_P0 = NamedTuple(build_para(build_vbpara()))
const T0 = vr -> NT_P0.l / vr            # normalized period -> physical seconds
const S9 = [-0.7771, 1.3899, 1.8811, 2.82632]   # stable operating cycle @ v_ref = 9

fd_jac(f, x; h=1e-7) = (n = length(x); J = Matrix{Float64}(undef, n, n);
    for i in 1:n
        e = zeros(n); e[i] = h
        @views J[:, i] .= (f(x .+ e) .- f(x .- e)) ./ (2h)
    end; J)

"Floquet multipliers of the cycle `s` at parameters `vbp` (4×4 period-map FD)."
function multipliers(s, vbp; tol=1e-11)
    u0, T = unpack_shooting(s)
    pm = v -> begin
        sol = integrate(SA[v[1], v[2], v[3], v[4], 0.0], T, vbp; tol=tol)
        copy(Array(sol.u[end])[1:4])
    end
    return eigvals(fd_jac(pm, collect(Array(u0)[1:4])))
end

"Largest |mu| excluding the trivial multiplier (the one closest to +1)."
function mmax(mu)
    j = argmin(abs.(mu .- 1))
    return maximum(abs(mu[k]) for k in eachindex(mu) if k != j)
end

"Warm-started continuation of the cycle in v_ref, returning rows."
function sweep(s0, vrefs)
    rows = NamedTuple[]
    s = s0
    for vr in vrefs
        vbp = build_vbpara(merge(NT_P0, (v_ref=vr,)))
        s = build_shooting(compute_limit_cycle(unpack_shooting(s)[1], vbp; tol=1e-11))
        mu = sort(abs.(multipliers(s, vbp)), rev=true)
        push!(rows, (vr=vr, tf=s[4] * T0(vr), mmax=mmax(multipliers(s, vbp))))
        flush(stdout)
    end
    return rows
end

## ---------------------------------------------------------------- A ---------
println("=== A. ground-truth |mu| of the tracked stable sheet, v_ref ∈ [8.8, 10.2] ===")
println("     (a genuine fold of this sheet needs max nontrivial |mu| -> 1)")
rows_up = sweep(S9, collect(9.0:0.02:10.2))
rows_dn = sweep(S9, collect(9.0:-0.02:8.8))
@printf("  v_ref     tf (s)     max|mu|!=1\n")
for r in vcat(reverse(rows_dn[2:end]), rows_up)
    @printf("  %6.2f   %9.5f   %.3e\n", r.vr, r.tf, r.mmax)
end
tf_up = [r.tf for r in rows_up]          # v_ref increasing 9 -> 10.2
tf_dn = [r.tf for r in rows_dn]          # v_ref decreasing 9 -> 8.8
println("  tf monotone-decreasing over [8.8, 10.2] ? ",
    issorted(tf_up, rev=true) && issorted(tf_dn))
println("  max nontrivial |mu| over the band = ",
    maximum(r.mmax for r in vcat(rows_dn, rows_up)))

## ---------------------------------------------------------------- B ---------
println("\n=== B. independent orbit census (Poincaré sampling), v_ref ∈ [9.0, 10.0] ===")
for vr in 9.0:0.25:10.0
    vbp = build_vbpara(merge(NT_P0, (v_ref=vr,)))
    lcs = all_limit_cycles(vbp; αmin=-π, αmax=π, vmax=12, N=40)
    println("v_ref = ", vr, " : ", length(lcs), " cycle(s)")
    for lc in lcs
        s = build_shooting(lc)
        mus = multipliers(s, vbp)
        @printf("    α0=%+8.5f  dα0=%+8.5f  dτ0=%+8.5f  tf=%8.4f s  max|mu|!=1=%.3e\n",
            s[1], s[2], s[3], s[4] * T0(vr), mmax(mus))
    end
end

## ---------------------------------------------------------------- C ---------
println("\n=== C. archived BK collocation trace (scratch/brA_ds002.jls) ===")
bA = deserialize(joinpath(@__DIR__, "scratch", "brA_ds002.jls"))
m = (bA.p .> 8.8) .& (bA.p .< 10.4)
pp, tt, dd = bA.p[m], bA.tf[m], bA.ds[m]
println("  points in band: ", length(pp), "   ds range: [",
    round(minimum(dd), sigdigits=3), ", ", round(maximum(dd), sigdigits=3), "]")
# a PALC fold shows a sign change / reversal of consecutive Δp at roughly-constant ds
Δp = diff(pp)
revs = [i for i in 2:length(Δp) if sign(Δp[i]) != sign(Δp[i-1])]
println("  parameter reversals (sign change of Δp) at indices: ", revs)
for i in revs
    @printf("    p=%.5f  tf=%.5f  ds=%.3e\n", pp[i], tt[i], dd[i])
end
println("  min |ds| in band = ", minimum(abs.(dd)),
    "   ds collapsed (< 1e-3) ? ", minimum(abs.(dd)) < 1e-3)
println("  p monotone across band ? ", all(Δp .>= -1e-9) || all(Δp .<= 1e-9))

## ------------------------------------------------------------- REPORT -------
# VERDICT: the "fold tangle" near v_ref ≈ 9.2–9.8 m/s is a CONTINUATION
# ARTIFACT, not a physical saddle-node, and it hides no new orbits.
#
# Evidence (all at default params, the setting the tangle was reported for):
#   A. Ground-truth Floquet multipliers of the tracked stable sheet are
#      ≤ 4.8e-3 across v_ref ∈ [8.8, 10.2] (rising smoothly with v_ref) —
#      nowhere near the +1 a genuine fold requires. tf(v_ref) is smooth and
#      monotone decreasing (0.657 s @ 8.8 → 0.489 s @ 10.2): no turnaround.
#   B. Independent Poincaré census finds exactly 2 cycles at every v_ref in
#      [9.0, 10.0]: the known α ↔ −α mirror pair, both stable. No subharmonics,
#      no extra sheets.
#   C. The archived BK collocation trace in the band shows parameter reversals
#      but NO step-size collapse (|ds| ∈ [0.005, 0.02]); the reversals sit at
#      the `bothside` seam (p = 9.0, ds = ∓0.005) plus a mild wobble near 9.21.
#      With A excluding a multiplier → +1, these are PALC/corrector artifacts on
#      a nearly flat segment, not folds.
#
# What the sheet-jumps actually hop between: the two symmetry-related cycles of
# (B), not distinct physical branches. So the tangle caveat in REPORT.md is a
# continuation diagnostic, not a physical feature.
#
# Run: julia +1.12.7 --project=. _fold_tangle_probe.jl
# (juliaup's `release` channel moved to 1.13.0, which cannot use this project's
#  1.12 images; pin 1.12.7 explicitly.)


