# Wind scaling and periodic-orbit continuation of the point-mass kite

Result of the ECC2026 parametric study (`applications/sprint2026`, optimized
kite of `BK_parametric.jl`; reference wind `v_ref = 9 m/s`). The objective of
`6_fixed_parametrized_optimization.jl` is to probe the system across wind
speeds with every other physical parameter held fixed: the cycles at different
winds are therefore expected *not* to be self-similar — that is the intended
physics, and the fixed implementation realizes it exactly. The deliverable
code is `BK_parametric_fixed.jl` and
`test/publications/ECC2026/6_fixed_parametrized_optimization_rescaling.jl`.

## 1. Verification: the wind enters alone (Vaschy–Buckingham control)

The kite dynamics are self-similar along the fibers of `build_vbpara`: there is
a *unique* parameter transform that multiplies the wind by `k` while leaving
every dimensionless group invariant. Gravity pins the gauge: the Froude group
forces `λ = k²`, i.e. a full Froude rescaling of the hardware,

| quantity | law | | quantity | law |
|---|---|---|---|---|
| `l, r, h_ref` | ∝ k² | | `I_eq` | ∝ k⁴ |
| `Cmax` | ∝ k² | | `Ω`, `d_l`, `S` | ∝ k⁻¹, k⁻², k⁻⁴ |
| `torque_slope` | ∝ k³ | | `ρ_l` | ∝ k⁶ |

with `m` fixed. On the family of periodic orbits this implies **period ∝ k,
mean power ∝ k, line tension invariant** (physical gauge), verified to ~1e-8
relative (`rescale_para` / `rescale_wind` in the test file, Tests A–C: 32/32
pass).

The fiber is the *control* that certifies the wind sweep. A wind-only change
at fixed physical parameters is not on the fiber: normalized gravity, line
mass, `Cmax`, pump rate and torque slope (each ∝ 1/k² or 1/k) all shift with
the wind, so different winds are genuinely different normalized problems and
the cycles should *not* overlay — non-self-similarity is the expected,
correct outcome. The fixed physical-time BVP (`BK_parametric_fixed.jl`)
implements exactly that: `v_ref` changes only the wind profile. The figure
contrasts the two admissible situations at k = 2: the exact `rescale_wind`
fiber overlays the reference orbit (α(τ) to 2.6e-7; period, power and tension
ratios = 1 to ~1e-8), certifying the pipeline, while the old `v_ref` sweep
leaves the fiber — as any wind-only change must. What disqualified the old
script was not its objective but its implementation: the RHS evaluated the
dynamics through `build_vbpara`, which renormalizes with `T(v_ref)`, so
continuing in `v_ref` rescaled the state units *in addition* to the wind — the
sweep implemented neither the wind-only physics nor the fiber.

![VB rescaling](fig_vb_rescale.png)

*Top: normalized phase portrait at k = 2 — the exact `rescale_wind` fiber
overlays the reference; the old `v_ref`-only sweep does not, as any wind-only
change should. Bottom: similarity-law ratios — gray = exact fiber (laws hold),
colored = old sweep, correctly off the fiber (its actual flaw is the unit
leak described above).*

## 2. Family of periodic orbits vs wind speed

With the wind entering only physically (fixed `BK_parametric_fixed.jl`), the
optimized cycle continues as the family below (`tf` = physical period):

![Branch topology](fig_branch.png)

* **Main sheet** (top panel): `tf` falls monotonically from 0.628 s at 9 m/s to
  ≈ 0.12 s at 20 m/s. Collocation, multiple shooting and the bothside runs all
  coincide on it.
* **Fold tangle at v_ref ≈ 9.2–9.8 m/s** (shaded): several periodic-orbit
  sheets crowd together (detected folds at 9.24, 9.30, 9.46, 9.54, 9.84 m/s).
  Near it, Newton's corrector can land on a different sheet than the predictor —
  the red ✕ markers flag recorded steps with |Δp| > 1.5·dsmax, impossible for a
  true PALC step and therefore the signature of a sheet-jump.
* **Low-wind fold** (bottom): the period grows steeply as the wind decreases
  (1/(1/k) law-ish), but it does **not** diverge: the family folds at
  `v* ≈ 2.4586 m/s` where `tf ≈ 49 s` and the branch turns back toward higher
  wind. Both legs overlap in the (v_ref, tf) projection.

## 3. Collocation and multiple shooting agree

Both discretizations trace the same sheets. On sheet-matched, smooth segments
the relative period difference is **≤ 3.0e-5 on the main sheet** and **≤ 8.6e-4
on the low-wind sheet** (where tf is large and sensitive).

![MS vs collocation](fig_ms_vs_coll.png)

The historically reported discrepancy ("multiple-shooting values are much
smaller") had two independent, non-physical causes, both fixed in
`BK_parametric_fixed.jl`:

1. **Record function.** Without an explicit `record_from_solution`,
   `plot(br)` uses `norm(x)` — a norm of different unknown-vector layouts
   (collocation: `n·(1+Ntst)`, shooting: `nM+1`), hence incomparable y-axes.
   Fix: `record_from_solution = (x, p; kwargs...) -> raw_x(x)[5]`, the physical
   period, with `raw_x` unwrapping the `BVPSavedSolutionAndState` wrapper.
2. **Branch ordering.** With `bothside=true`, BifurcationKit merges the legs as
   `_cat!(_reverse(backward), forward)` (`Results.jl:464`): the stored branch
   *begins at the backward leg's endpoint*, so the plotted "second branch" is a
   reversed, sheet-jumped leg — indistinguishable from a second family. The
   fixed script replaces `bothside=true` by two explicit single-direction
   continuations (each with its own budget and a `|Δp| ≤ dsmax` sanity check).

## 4. Practical recipe (what the fixed script does)

1. **Physical-time BVP**: state `(α, τ, α̇, τ̇, tf)` in rad, rad/s, s; the
   characteristic time `T = l/v_ref` is used only to convert to the internally
   normalized dynamics, so `v_ref` changes *only* the wind magnitude.
2. **Deterministic multiple shooting** despite the BifurcationKit 0.8.2
   residual bug (`out[end]` uninitialized): override `bvp_residual` to pin the
   auxiliary unknown (`out[end] = X[end]`, Jacobian row `J[end,end] = 1`).
3. **Explicit record** of the physical period (above).
4. **Two explicit single-direction continuations** (`ds = ±0.01`) instead of
   `bothside=true`, each with its own budget, plus a `sheet_jumps()` sanity
   check (|Δp| > 1.5·dsmax ⇒ Newton left the sheet near the fold tangle).
5. Warm start: damped-Armijo pre-solve Newton (BK's plain `Newton` overshoots),
   multiple shooting warm-started from the converged collocation orbit.
   Integrator: **Tsit5** with abstol = reltol = 1e-10 (an order of magnitude
   faster than the earlier Vern9 setup at equal accuracy).

Runtimes (MacBook Pro, Apple Silicon): FastPara-typed RHS 387 ns/eval
(16.4× faster, 0 allocs), collocation residual 70.7 µs (13.2×), full script
~2.5 min including IPOPT optimization and both continuation pairs.

## 5. Reproducibility

| file | role |
|---|---|
| `applications/sprint2026/BK_parametric_fixed.jl` | deliverable: optimization → physical-time BVP → collocation + multiple-shooting continuation |
| `test/publications/ECC2026/6_fixed_parametrized_optimization_rescaling.jl` | VB derivation, `rescale_para`/`rescale_wind`, Tests A–C (32/32) |
| `applications/sprint2026/_order2_tmp.jl` | topology diagnostics behind figures 1–2 (CSVs in the scratch dir) |
| `fig_branch.png`, `fig_ms_vs_coll.png`, `fig_vb_rescale.png` | figures shown above |
| `applications/sprint2026/REPORT.md` | this report |

---

## How the right code came about

The final form of `BK_parametric_fixed.jl` resulted from a chain of five
distinct discoveries, each forced by an observed symptom:

1. **The original sweep's implementation was not physical.** Its objective —
   vary the wind with every other physical parameter fixed, yielding
   intentionally non-self-similar cycles — was sound, but the RHS evaluated
   the dynamics through `build_vbpara`, which renormalizes with `T(v_ref)`: the
   continuation acted in Vaschy–Buckingham space, silently changing normalized
   gravity, inertia and torque *in addition* to the wind. The fix keeps the
   objective and changes the implementation: a state written in *physical*
   units so that `v_ref` shifts only the wind profile, with the VB fiber
   (gravity forcing `λ = k²`) kept in the ECC test file as the §1 verification
   control.

2. **Multiple shooting gave nondeterministic, and seemingly "smaller",
   results.** Two unrelated bugs: (a) BifurcationKit 0.8.2's `Shooting`
   residual leaves `out[end]` uninitialized (`out = similar(X)`), feeding
   Newton random memory — fixed by overriding `bvp_residual` and pinning the
   auxiliary unknown to zero; (b) `plot(br)` defaults to `norm(x)` as record,
   a different number for each discretization's unknown layout — the
   "MS values are much smaller" was `norm(x)` collocation vs a *different*
   default for the ODEProblem-wrapped shooting, not dynamics.

3. **"Both branches go right from v_ref = 9" was a plotting artifact.**
   `bothside=true` merges its two legs as `[reverse(backward) ++ forward]`
   (`Results.jl:464`), so the *backward* leg appears first, reversed, starting
   at its endpoint — indistinguishable from a second family going right.
   Single-direction runs behave perfectly; the merged result does not.

4. **Sheet jumps near v_ref ≈ 9.2–9.8.** Recorded steps with |Δp| > dsmax
   (impossible for PALC) revealed Newton corrections landing on adjacent sheets
   in the fold tangle, then riding the wrong sheet to the budget. Fresh-mesh
   runs jump; runs with mesh history reached the low-wind sheet — the tangle
   makes the corrector basin history-dependent. This motivated the explicit
   `sheet_jumps()` check and replacing `bothside=true` with two explicit runs.

5. **The "low-wind divergence" is a fold.** Wide-budget runs stalled at
   v_ref = 2.4585 with huge tf — read as divergence — but the 200-step multiple
   shooting run passed v* = 2.4586, peaked at tf ≈ 49 s and turned back:
   a fold, not an endpoint. Faster integration (Tsit5 at 1e-10) made this
   affordable at all.

Methodological choices that stuck: fast typed parameter struct + `v_ref`-keyed
cache (84% of residual cost), container-matched RHS (`SVector`/`Vector`),
damped-Armijo pre-solve Newton with step capping, warm-starting MS from the
converged collocation orbit, and doing all diagnostics on *one* persistent
parameter set (snapshot `_opt_result_tmp.jl`) so any discrepancy is
methodological, not setup noise.
