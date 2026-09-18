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
| `applications/sprint2026/BK_tests_0910.jl` | diagnostics: method switch, manual block shooting Jacobian (`ManualJacFwd`/`ManualJacFD`), ground-truth Floquet (`floquet_true`, segment-wise `monodromy`), fold-vs-homoclinic fit (`fit_fold`/`fit_homoclinic`), coexisting cycles, equilibria/Hopf (Findings) |
| `test/publications/ECC2026/6_fixed_parametrized_optimization_rescaling.jl` | VB derivation, `rescale_para`/`rescale_wind`, Tests A–C (32/32) |
| `applications/sprint2026/scratch/` | archived investigation traces (superseded `_*_tmp.jl` / `_*_keep.jl` diagnostics + their `.jls` data); topology diagnostics behind figures 1–2 in `scratch/_order2_tmp.jl` |
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
parameter set (snapshot `scratch/_opt_result_tmp.jl`) so any discrepancy is
methodological, not setup noise.

---

## 3. The two folds of the low-wind branch, the arc-count sweep, and coexisting cycles

### 3.1 Arc-count sweep (multiple shooting)

`scratch/_arcs_sweep_one_tmp.jl` (fresh collocation warm start per M, one
bothside continuation, loose tol `1e-6`, `max_steps=400`):

| M | passes fold | long-leg tfmax | steps | wall | peak RSS* |
|---|---|---|---|---|---|
| 5 | no (p_min=2.876) | 6.3 | 621 | 16.8 s | — |
| 10 | yes | 39.7 | 802 (budget) | 28.7 s | — |
| 20 | yes | 58.8 | 688 | 22.2 s | — |
| 40 | yes | 71.8 | 691 | 16.2 s | — |
| 80 | yes | 111.4 | 802 (budget) | 29.7 s | — |
| 100 | yes, but stalls on `dsmin` | 59.1 | 655 | 59.8 s | 1412 MiB |

\* cumulative process peak (Julia + BK + Plots). At full tolerance
(`scratch/_branch_fulltol.jl`, `max_steps=400`) both M=40 and M=80 reach
`tfmax ≈ 111.4` and are budget-limited, so the finer method does **not** extend
the reach at equal budget; M=100 is worse (smaller reach, slower) because its
corrector underflows `ds`. **Chosen base M = 40** (smallest clean pass, cheapest
per step); the extra arcs buy conditioning, not reach, here.

### 3.2 The two turns are genuine quadratic saddle–nodes of the prograde branch

From `scratch/brS_shoot.jls` (`i0 = argmin p`, `p_f = 2.458474`, `tf_f = 14.296`):

- `p` has an *interior* minimum — strictly decreasing on the long leg, strictly
  increasing on the short leg. A homoclinic/SNIC would instead have
  `tf → ∞` with `p` monotone and no turn.
- Local log–log exponent `(p − p_min) ∝ |tf − tf_f|^α`: **α = 2.10** (long leg),
  **α = 2.03** (short leg) — quadratic tangency.
- Tightest window (long leg, tf∈[14,15], 5 pts): fold law `p = p_f + γ(tf−tf_f)²`
  gives RMS **4.6e-7 m/s**; the homoclinic log-law **9.6e-5** (200× worse).

So "no prograde cycles below `v_ref*`" is correct **for this branch**.

**Second fold.** The long-period side does not run away monotonically: on the
full-tol branches the branch turns a *second* time, at a `p` **maximum**
(`i2 = argmax p[1:i0]`). Identical at M=40 and M=80, so it is genuine, not a
discretization artifact:

| fold | `tf` [s] | `v_ref` | type |
|---|---|---|---|
| 1 | 14.44 | 2.45848 | `p` min (low-wind saddle–node) |
| 2 | 67.67 | 2.48168 | `p` max (long-sheet saddle–node) |

The branch is therefore an **S-curve**: short leg (tf 0.07→14.44), long sheet
(14.44→67.67, `p` rising), and a **long-long sheet** (67.67→111+, `p` falling,
still budget-truncated at `tf≈111.4`). For every `v_ref ∈ (2.45848, 2.48168)` a
horizontal cut meets it **three times**, i.e. three coexisting prograde periods
(see §3.3).

### 3.3 The Poincaré sampler reproduces the branch (a units trap), and three periods at the second fold

The BVP boundary condition (line 167 of `BK_tests_0910.jl`) enforces
`τ(tf) − τ(0) = +2π`, so a cycle "verifies the BVP" iff it is **prograde**
(`dτ0 > 0`). A retrograde cycle (`Δτ = −2π`) is a genuine flow solution but
**not** a BVP solution, so it is out of scope by construction.
`KEEP.LimitCycle.all_limit_cycles` (`scratch/_bvp_cycles_poincare_tmp.jl`,
`N=40`, `vmax=12`, `α∈[−π,π]`, filtered to `build_shooting(lc)[3] > 0`)
**returns NORMALIZED units** (L=2 m, M=6 kg, time `T0 = l/v_ref = 2.0/v_ref` s,
velocities physical×`T0`). Converting to the branch's PHYSICAL/SI convention,
`tf_SI = T·T0`, `dα_SI = dα0/T0`, `dτ_SI = dτ0/T0`:

| `v_ref` | cycles | `T` (norm.) | `tf = T·T0` (SI) | branch short `tf` (SI) | `α0` | `dα0`(SI) | `dτ0`(SI) | power [W] |
|---|---|---|---|---|---|---|---|---|
| 2.42 | **0** | — | — | (below fold) | — | — | — | — |
| 2.45 | **0** | — | — | (below fold) | — | — | — | — |
| 2.4585 | **0** | — | — | 14.12 (at fold) | — | — | — | — |
| 2.46 | 1 | 15.6537 | 12.727 | 12.858 (@2.4597) | −0.219 | +0.041 | +0.655 | 42 |
| 2.47 | 1 | 14.0367 | **11.366** | **11.3658** | −0.225 | +0.044 | +0.671 | 48 |
| 2.50 | 1 | 12.6874 | 10.150 | 10.025 | −0.233 | +0.049 | +0.695 | 56 |
| 2.60 | 1 | 11.0111 | 8.470 | 8.467 | −0.251 | +0.062 | +0.760 | 75 |
| 2.80 | 1 | 9.4972 | 6.784 | 6.652 | −0.279 | +0.089 | +0.886 | 110 |
| 3.00 | 1 | 8.5795 | 5.720 | 5.790 | −0.303 | +0.121 | +1.017 | 150 |
| 3.50 | 1 | 7.1232 | 4.070 | 4.019 | −0.359 | +0.230 | +1.377 | 287 |
| 4.00 | 1 | 6.1822 | 3.091 | 3.058 | −0.410 | +0.388 | +1.788 | 487 |
| 5.00 | 1 | 4.9596 | 1.984 | 1.965 | −0.503 | +0.892 | +2.758 | 1140 |
| 7.00 | 1 | 3.6018 | 1.029 | 1.016 | −0.659 | +2.843 | +5.261 | 3847 |
| 9.00 | 1 | 2.8263 | 0.6281 | 0.6281 | −0.777 | +6.255 | +8.465 | 9055 |

Data: `scratch/bvp_cycles_poincare.jls`; branch column from
`scratch/brS_shoot_M40_fulltol.jls`. Consequences:

- **The sampler finds the blue short branch, exactly.** At `v_ref = 2.47` the
  converted sampler cycle `(α0,dα0,dτ0,tf) = (−0.225, +0.044, +0.671, 11.366)`
  equals the branch's `(−0.22521, +0.04439, +0.67054, 11.36577)`; seeding the
  *normalized* Poincaré flow from the branch state returns the sampler's cycle
  to 5 digits. An earlier figure overlaid `T` (normalized) against the branch's
  `tf` (SI) and concluded, wrongly, that the sampler traced a **separate family**
  with ≥3 coexisting prograde cycles. There is no separate family — that was a
  unit mismatch. `presentation_figs.jl` no longer overlays it.
- **No prograde cycle below the fold**: the sampler finds none at
  `v_ref = 2.42, 2.45, 2.4585`, exactly as the fold requires.
- **The sampler is attractor-only.** Its callback terminates when successive
  section crossings agree, so it only sees *attracting* cycles. It therefore
  finds the stable short branch and **misses** the long / long-long saddle
  sheets. Enumerating those needs deflation (§3.4).
- **Three periods coexist** for `v_ref ∈ (2.45848, 2.48168)` (the two folds of
  §3.2). A cut at `v_ref = 2.480` meets the S-curve three times — short
  `10.83 s`, long `49.22 s`, long-long `96.35 s` (§3.5).
- A retrograde attractor also exists below the fold (`T ≈ 8.5 s`; diagnostic
  `scratch/attractor_family.jls`, probe `scratch/_attractor_probe_tmp.jl`). It is
  a real limit cycle but fails the BVP boundary condition, so it is not counted.

### 3.4 Deflation: recovering the non-attracting cycles

The Poincaré sampler is attractor-only, so the saddle sheets are invisible to it.
Deflation (Farrell–Birkisson–Funke, `BifurcationKit.DeflationOperator`, penalising
`‖u−u_i‖^{-2p}+α`) is what reaches them. `scratch/_deflate_cycles.jl`
(M=20, physical/SI BVP) at the two anchors:

| `v_ref` | # cycles | `tf` [s] | `α0` [rad] | max non-trivial \|μ\| | power [W] | cycle res |
|---|---|---|---|---|---|---|
| 2.47 | 2 | 11.3658 | −0.22521 | 2.94e−01 | 47.3 | 1.6e−16 |
| 2.47 | | 28.8205 | −0.19298 | 3.47e+16 | 24.4 | 1.9e−12 |
| 9.00 | 1 | 0.6281 | −0.77710 | 2.80e−03 | 8736.6 | 1.8e−15 |

- **Deflation demonstration (v_ref = 2.47):** with the *short* cycle deflated, a
  deflated Newton started from the long state converges back to the **long
  saddle** (`tf = 28.8205 s`) — i.e. deflation does reach the non-attracting
  sheet, which the sampler cannot. The archived short/long periods are
  reproduced exactly (`11.3658` / `28.8205`), and the freshly regenerated short
  from the sampler agrees (`tf = 11.36575`), so the archived orbits are not stale.
- **No extra cycle found.** Searching beyond the seeds (perturbed states) is not
  usable with this machinery: a divergent shooting guess makes the adaptive
  integrator crawl and `with_timeout` cannot interrupt it, so the search was
  abandoned. A complete enumeration must be seeded from continuation-derived
  states (see §3.7).
- **Reference wind `v_ref = 9`** (the optimization wind, `build_vbpara()` default)
  has exactly **one** prograde BVP cycle — the operating cycle. So a deflated
  Newton there returns nothing new to seed an optimization with.

### 3.5 Three coexisting periods at the second fold

Within `v_ref ∈ (2.45848, 2.48168)` the S-curve is cut three times. Crossings
from the M=40 full-tol branch (`tf` interpolated at fixed `v_ref`):

| `v_ref` | cycle 1 (short) | cycle 2 (long) | cycle 3 (long-long) |
|---|---|---|---|
| 2.470 | 11.376 | 28.821 | — (not reached in budget) |
| 2.478 | 10.919 | 42.442 | — |
| **2.480** | **10.833** | **49.216** | **96.348** |
| 2.4810 | 10.790 | 55.085 | 84.181 |
| 2.4816 | 10.764 | 62.887 | 72.702 |

The three merge at the second fold (`tf = 67.67 s`) and cycle 1's value
`10.833 s` matches the independent deflation/sampler result `10.8232 s`. The
long-long sheet is **budget-truncated** at `tf ≈ 111.4 s` (it is still
descending in `v_ref`), not closed by a fold, so a third fold or a homoclinic
terminus beyond `tf = 111` is unresolved. Its stability was not measured (the
state would have to come from a continuation); only the two 2.47 sheets were
deflated.

### 3.6 Two-cycle extraction at M=40 (and a silent non-convergence, fixed)

`scratch/_regen_two_cycles_tmp.jl` (`MSTAR=40`, Newton tol `1e-10`) re-extracts
the `v_ref = 2.47` pair at 40 arcs. First attempt gave short `tf = 11.36577 s`
(reproduces the archive) but long `tf = 17.15 s`, not the archived `28.82 s`.
The `17.15` orbit did **not** close — cycle residual `2.1e-3`
(`|α(tf)−α(0)|` etc.) — i.e. the fixed-parameter Newton from a near-fold
candidate failed to converge and the script reported it anyway. Diagnosis: at
`v_ref = 2.47` the long branch crosses once, at `tf ≈ 29 s` (from the M=40
full-tol branch: `p = 2.47014` at `tf = 28.998`), but `longs[1]` picked the
smallest candidate above `1.3·tf_short` — a near-fold point whose refine did not
close. The archive's `28.82 s` orbit closes to `5.7e-15`, so it was right and the
M=40 run was a tolerance/candidate-selection failure, not a different solution.

Fix: `_regen` now tries candidates ordered by proximity to `v_ref` and accepts
only an orbit that closes (cycle residual < `1e-7`). Rerun at M=40:

| cycle | tf [s] | max non-trivial \|μ\| |
|---|---|---|
| short | 11.36577 | 0.2943 |
| long | 28.82052 | 3.47e16 |

both closing to ~`1e-14` and matching the archived pair (Finding 9's long-sheet
`|μ| ≈ 3.5e16`). So the documented pair is confirmed at M=40, and
`two_cycles_shooting_M40.jls` / `two_cycles_M40_validated.jls` now hold it; the
figures keep using `two_cycles_shooting.jls` (identical values).

### 3.7 What's next

- **Optimization seeding.** The optimization is warm-started from a single
  operating orbit. If a deflated Newton at the reference wind (`v_ref = 9`)
  returned a *distinct, stable, higher-power* cycle, it could seed a second
  `optimize()` run to test for a better operating point. At present deflation
  finds only the known operating cycle there, so there is nothing to seed with —
  but the hook is cheap to add once a richer anchor (e.g. the 2.48 three-cycle
  window) yields a stable alternative.
- **Full deflation enumeration.** The perturbation search was abandoned because
  divergent shooting guesses stall the integrator. A complete sweep should seed
  from *continuation-derived* states on every sheet (short, long, long-long at
  2.48) and deflate them jointly, then probe with continuation-adjacent guesses —
  bounded per attempt (a cooperative cancellation hook, or a per-attempt thread
  kill).
- **Third fold / homoclinic.** The long-long sheet is budget-truncated at
  `tf ≈ 111 s` and still descending in `v_ref`; a longer continuation would show
  whether it folds a third time (a fourth coexisting period) or runs to a
  homoclinic terminus.
