# Optimization over Θ = (ℓ, I, b, Δφ)

Report for `Section_Optimization.jl`. All numbers come from a full run of that
script at the reference wind speed 9 m/s, tolerance `1e-5`. The kite
acceleration is the norm of the Cartesian acceleration of the kite,
`a = ‖d²x/dt²‖` in m s⁻², computed from the 4D model at each state of the
limit cycle (`PointMass4.kite_acceleration`).

## Summary

- Optimizing Θ = (ℓ, I, b, Δφ) at 9 m/s raises the average power from
  **9.64 kW** (kinetic-criterion start) to **12.25 kW** (+27.1 %). All four
  variables end strictly inside their boxes: Δφ = 16.6°, ℓ = 27.0 m,
  I = 352 kg m², b = 737.
- Re-optimizing Θ at each wind speed gains up to **+13.5 % at 5 m/s** and
  **+8.6 % at 15 m/s**, and nothing at the 9 m/s reference by construction.
- Δθ is not a power-relevant variable: over ±10 % it changes the power by at
  most 3.1 %, while it changes the peak acceleration by up to 20 % and the
  tension spread by 14 %. Holding it fixed is a load/trajectory choice, not a
  lost power opportunity.
- The 10×10 shape map shows that the shape choice dominates both power and
  loads: at the optimal (ℓ, I, b), power peaks at 12.46 kW near Δφ = 19° and
  falls to 1.2–6.4 kW at the narrowest shape and to 2–3 kW for Δφ ≥ 36.7°,
  while the peak acceleration ranges from 129 to 4731 m/s². At Δφ = 90° no
  pumping cycle exists at all.

## Setup

The optimization variables are

```
Θ = (ℓ, I, b, Δφ) = (l, I, b, Δφ).
```

Here `ℓ` (`p.l`) is the **tether length**, `I` is the **inertia of the rotating
arm** about its pivot, `b` is the generator braking slope and `Δφ` is the
lemniscate half-width. The arm length itself is fixed at `p.a = 2 m`. The
lemniscate half-height Δθ is held at its nominal 6° (0.10472 rad) and is never
optimized.

The boxes are

| variable | lower | upper |
|---|---:|---:|
| ℓ (m) | 15 | 200 |
| I (kg m²) | 62 | 1e6 |
| b (kg m² s⁻¹ rad⁻¹) | 1 | 1e6 |
| Δφ (°) | 10 | 90 |

The inertia floor is a structural one. The 2 m aluminium arm (`p.a` is the
lever arm) sees the worst line tension in the study, `T ≈ 34 kN` at 15 m/s, so
the root bending moment is `M = T p.a ≈ 68 kN m`. For 6061-T6
(`σy ≈ 276 MPa`) with a safety factor of 2, the required section modulus is
`W = M/σallow ≈ 491 cm³`, met by a Ø250 mm tube with an 11.5 mm wall (47 kg).
Its inertia about the pivot is `I = m p.a²/3 ≈ 62 kg m²` (a solid bar would be
165 kg m²). There is deliberately no structural upper bound: extra mass is
allowed if it buys power.

The objective is the average mechanical power over one limit cycle. The
reference wind is 9 m/s. The starting point is the default parameters with
`b = I ω` (`PMP.optimal_torque_slope`), not the default torque slope.

## Reference optimum at 9 m/s

| | init | Θ optimum |
|---|---:|---:|
| ℓ (m) | 31.0 | 26.97 |
| I (kg m²) | 1300 | 351.6 |
| b (kg m² s⁻¹ rad⁻¹) | 2888.3 | 737.1 |
| Δθ (rad) | 0.10472 | 0.10472 (fixed) |
| Δφ (rad) | 0.872665 (50°) | 0.290 (16.62°) |
| period (s) | 2.828 | 0.949 |
| average power (kW) | 9.639 | 12.253 |
| min line tension (kN) | 3.356 | 5.748 |
| max line tension (kN) | 8.856 | 10.827 |
| max kite acceleration (m/s²) | 517.1 | 938.6 |

No box constraint is active: Δφ sits at 16.6° between 10° and 90°, ℓ at 27 m
between 15 and 200 m, b at 737 between 1 and 1e6, and I at 352 kg m², 5.7 times
the structural floor of 62 kg m². The power optimum therefore prefers an arm
heavier than strength requires. The tether shortens, the inertia and the
period fall by factors of 3.7 and 3.0, and the faster cycle extracts 27 % more
power from the same wind.

## Wind sweep

![wind sweep](figs/optimization_wind_sweep_theta.pdf)

The dotted black curve is the 9 m/s optimum held fixed; the blue curve
re-optimizes Θ inside the same boxes; the grey dashed line is 9 m/s.

| v (m/s) | P fixed (kW) | P opt (kW) | gain | Tmax opt (kN) | ℓ (m) | I (kg m²) | b | Δφ (°) |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 5 | 1.894 | 2.149 | +13.5 % | 3.49 | 30.34 | 451.3 | 497.8 | 17.60 |
| 6 | 3.463 | 3.675 | +6.1 % | 4.93 | 28.78 | 405.4 | 551.3 | 17.22 |
| 7 | 5.672 | 5.801 | +2.3 % | 6.63 | 27.89 | 379.1 | 610.0 | 16.95 |
| 8 | 8.585 | 8.627 | +0.5 % | 8.60 | 27.34 | 362.7 | 672.3 | 16.76 |
| 9 | 12.253 | 12.253 | 0.0 % | 10.83 | 26.97 | 351.6 | 737.1 | 16.62 |
| 10 | 16.717 | 16.780 | +0.4 % | 13.32 | 26.71 | 343.9 | 803.7 | 16.51 |
| 11 | 22.008 | 22.306 | +1.4 % | 16.08 | 26.52 | 338.3 | 871.6 | 16.43 |
| 12 | 28.155 | 28.933 | +2.8 % | 19.09 | 26.37 | 334.0 | 940.5 | 16.37 |
| 13 | 35.180 | 36.759 | +4.5 % | 22.38 | 26.26 | 330.7 | 1010.2 | 16.33 |
| 14 | 43.103 | 45.885 | +6.5 % | 25.92 | 26.17 | 328.1 | 1080.4 | 16.29 |
| 15 | 51.940 | 56.412 | +8.6 % | 29.72 | 26.10 | 326.1 | 1151.2 | 16.26 |

The re-optimization adapts the tether and the inertia to the wind: ℓ shortens
from 30.3 to 26.1 m and I falls from 451 to 326 kg m² as the wind rises, while
the braking coefficient rises from 498 to 1151 and Δφ closes slightly from
17.6° to 16.3°. At 15 m/s the re-optimized cycle carries 10 % less maximum
tension (33.2 → 29.7 kN) and a 33 % higher minimum (12.2 → 16.1 kN); at 5 m/s
the maximum rises by 4.5 % (3.34 → 3.49 kN) and the minimum falls from 1.99 to
1.72 kN. The gain vanishes at 9 m/s because that is where the fixed parameters
were optimized.

With the wide boxes the objective is multimodal, so each point is first solved
from the 9 m/s optimum and retried from the neighbouring wind's optimum if that
fails. Only 7 m/s needs the retry; at 12 m/s the reverse happens, and a
continuation from 11 m/s alone would stall on a 3.3 kW non-pumping branch. The
sweep is a stationary solution path, not a global search.

## Objective landscape

![landscape](figs/optimization_landscape.pdf)

One parameter at a time around the optimum, the others fixed; the power is
normalized by its value at the optimum.

| parameter | at −10 % | at +10 % | box factor (÷, ×) |
|---|---:|---:|---:|
| ℓ | 0.885 | 0.919 | 0.556, 7.416 |
| I | 0.974 | 0.973 | 0.176, 2843.9 |
| b | 0.996 | 0.997 | 0.0014, 1356.7 |
| Δφ | 0.946 | 0.971 | 0.602, 5.417 |

ℓ is the stiffest direction (−11.5 % power at −10 %); I is about four times
flatter and b more than twenty times flatter. The Δφ direction loses 5.4 % at
15.0° and 2.9 % at 18.3°. No box constraint enters the ±10 % window, so none
of them is visible in the panes — the optimum is interior in all four
variables.

## Why Δθ is not optimized

![Δθ sensitivity](figs/optimization_dtheta.pdf)

Δθ swept ±10 % around nominal with Θ fixed at the optimum. The dashed line is
the nominal value.

| Δθ | P (kW) | Tmin (kN) | Tmax (kN) | Tmax − Tmin (kN) | max a (m/s²) |
|---:|---:|---:|---:|---:|---:|
| −10 % | 12.145 | 5.437 | 11.238 | 5.801 | 1035.1 |
| −8 % | 12.215 | 5.549 | 11.141 | 5.592 | 1011.6 |
| −6 % | 12.259 | 5.647 | 11.052 | 5.405 | 990.2 |
| −4 % | 12.279 | 5.675 | 10.970 | 5.296 | 970.9 |
| −2 % | 12.276 | 5.707 | 10.896 | 5.189 | 953.7 |
| 0 % | 12.253 | 5.748 | 10.827 | 5.079 | 938.6 |
| +2 % | 12.211 | 5.799 | 10.764 | 4.965 | 925.8 |
| +4 % | 12.152 | 5.861 | 10.705 | 4.844 | 957.3 |
| +6 % | 12.075 | 5.934 | 10.651 | 4.717 | 1011.7 |
| +8 % | 11.982 | 6.019 | 10.600 | 4.581 | 1068.8 |
| +10 % | 11.875 | 6.115 | 10.552 | 4.438 | 1128.1 |

Checking the expectations:

- **Power is flat — confirmed.** The whole ±10 % window spans 3.4 % of power;
  the maximum is 12.279 kW at −4 %, only 0.2 % above nominal.
- **Smaller Δθ lowers the minimum tension and widens the spread — confirmed,
  but the minimum does not approach zero.** Tmin falls 5.4 % (5.75 → 5.44 kN)
  and the spread grows 14 % at −10 %, yet the minimum stays at about half the
  maximum. The near-zero minima (down to 0.13 kN) appear at **large Δφ**, not
  at small Δθ — see the shape map below.
- **Smaller Δθ does not reduce the maximum tension — confirmed.** Tmax is
  3.8 % higher at −10 % (11.24 vs 10.83 kN).
- **Smaller Δθ gives a higher peak acceleration — only on one side.** At the
  optimum shape the peak acceleration is U-shaped in Δθ: 1035 m/s² at −10 %,
  a minimum of 925.8 m/s² at +2 %, and 1128 m/s² at +10 %. Reducing Δθ does
  raise the peak (the tight-turn intuition), but increasing it raises it too.
  The shape map shows the sign is set by the aspect ratio: at Δφ = 10° the
  peak is 3530 m/s² at Δθ = −20 % and 4731 m/s² at +20 %, while at
  Δφ = 18.9° it rises monotonically from 592 to 1343 m/s² over the same Δθ
  range.

**Conclusion.** Δθ is not a power variable (3.1 % at most), but it is a load
and trajectory variable (3.8 % in Tmax, 14 % in the tension spread, 20 % in
the peak acceleration). A power-only objective must not optimize it; the
nominal 6° is within 1.4 % of the acceleration minimum and within 0.2 % of the
power maximum. Choosing Δθ for loads or handling qualities would require a
multi-objective criterion.

## Shape map at the optimal ℓ, I, b

![shape map](figs/optimization_shape_grid.pdf)

10×10 limit cycles over Δθ ∈ [−20 %, +20 %] of nominal and Δφ ∈ [10°, 90°].
ℓ, I and b are held at the optimum; nothing is re-optimized. The white dashed
lines mark the nominal Δθ and the optimal Δφ = 16.6°. Of the 100 cells, 90
converge; the whole Δφ = 90° row has no dτ > 0 pumping cycle at these fixed
parameters.

| quantity | minimum | maximum |
|---|---:|---:|
| average power (kW) | 1.21 (Δθ −20 %, Δφ 10°) | 12.46 (Δθ −20 %, Δφ 18.9°) |
| max line tension (kN) | 10.06 (Δθ +15.6 %, Δφ 18.9°) | 14.96 (Δθ −20 %, Δφ 10°) |
| min line tension (kN) | 0.13 (Δθ +20 %, Δφ 81.1°) | 10.51 (Δθ −20 %, Δφ 10°) |
| max acceleration (m/s²) | 129 (Δθ +20 %, Δφ 81.1°) | 4731 (Δθ +20 %, Δφ 10°) |

The map separates four regimes:

- **Power is a Δφ effect with a sharp band near 19°.** It peaks at 12.4 kW
  there and collapses away from it: 1.2–6.4 kW at Δφ = 10°, 2–3 kW for
  Δφ ≥ 36.7°. Within the band, Δθ matters: at Δφ = 18.9° the power falls from
  12.46 to 9.90 kW as Δθ goes from −20 % to +20 %, and at Δφ = 10° it rises
  from 1.21 to 6.44 kW over the same sweep.
- **Maximum tension has two ridges.** It peaks at 14.96 kN at the narrowest
  shape (Δφ = 10°, Δθ = −20 %), has a secondary ridge of about 12.6 kN near
  Δφ = 36.7°, and is mildest along Δφ = 18.9° (10.1–10.6 kN).
- **The acceleration is largest at narrow shapes.** It reaches 3.5–4.7 km s⁻²
  at Δφ = 10°, forms a ridge of 1.1–1.8 km s⁻² at intermediate Δφ, and falls
  to 0.13–0.3 km s⁻² at Δφ = 81°. The power-optimal band sits on the mild side
  of that ridge.
- **The minimum tension is a large-Δφ effect.** It is 6–10.5 kN near the
  power optimum, falls through 1–2 kN at Δφ ≈ 72°, and reaches 0.13 kN at
  Δφ = 81° with Δθ = +20 %, where the kite slows at the ends of the wide
  sweep. This is where the "min tension goes to zero" expectation is realized
  — through Δφ, not through reducing Δθ.

At Δφ = 90° the fixed optimal (ℓ, I, b) admit no pumping cycle at any Δθ, so
the last row of the map is empty. That is a feasibility boundary of the device
parameters, not a plotting choice; a wider shape would require re-optimizing
the tether and the inertia.

## Discussion and limitations

- The optimum is interior in all four variables, so it is a local optimum of a
  multimodal objective rather than a box corner. The wide boxes make the
  solver's path matter: the wind sweep needs a retry at 7 m/s, and a pure
  continuation at 12 m/s would report a 3.3 kW branch. The reported gains are
  stationary-solution gains, not global optima.
- The inertia floor (62 kg m²) comes from a strength-sized 2 m aluminium arm
  and is inactive at the optimum (I = 352 kg m²). The power objective wants a
  heavier arm than strength requires; a mass-limited design would have to add
  a mass or stress constraint, not just a box.
- The shape map holds the device parameters fixed, which is exactly what the
  question "what does the shape choice cost?" asks, but it is a local
  sensitivity: a different shape could be paired with a different optimal
  tether and inertia.
- Only the dτ > 0 pumping cycle is considered. A dτ < 0 cycle also exists at
  the same parameters; it is not analyzed here.
- The 4D model is quasi-steady: no line dynamics, no structural flexibility,
  no unsteady aerodynamics. The acceleration is a model diagnostic, not a
  measured load; it includes the sharp turns but not jerk or tension rates.
- The power objective ignores loads. The shape map shows why a power-only
  optimum is not enough for a shape choice: the same power can be reached with
  very different tension and acceleration levels.

## Figures

| file | content |
|---|---|
| `figs/optimization_wind_sweep_theta.pdf` | wind sweep over Θ, three panels |
| `figs/optimization_power_theta.pdf` | wind sweep, power panel |
| `figs/optimization_tension_theta.pdf` | wind sweep, tension panel |
| `figs/optimization_params_theta.pdf` | wind sweep, parameter panel |
| `figs/optimization_landscape.pdf` | objective landscape around the Θ optimum |
| `figs/optimization_dtheta.pdf` | Δθ sensitivity: power, tension, acceleration |
| `figs/optimization_shape_grid.pdf` | 10×10 shape map over (Δθ, Δφ) |
| `figs/optimization_timeseries_init.pdf` | time series at the initial parameters |
| `figs/optimization_timeseries_opt.pdf` | time series at the Θ optimum |
| `figs/optimization_timeseries_all.pdf` | both time series side by side |
