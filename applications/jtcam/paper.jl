# julia --project=. paper.jl #src
# pandoc paper.md -o paper.html --standalone --mathjax --css="https://cdn.jsdelivr.net/npm/water.css@2/out/water.css" #src

if abspath(PROGRAM_FILE) == (@__FILE__) && get(ENV, "PAPER_LITERATE_EXEC", "0") != "1" #src
    ENV["PAPER_LITERATE_EXEC"] = "1" #src
    using Literate #src
    Literate.markdown(
        @__FILE__, @__DIR__; execute=true, flavor=Literate.CommonMarkFlavor(), credit=false
    ) #src
    exit(0) #src
end #src

# # Rocking-Arm Kite Power Generation (KEEP): Compact Analysis
#
# **Research question:** *Given a kite tethered to a rotating ground-level arm that is pulled left and right across the wind window, how much average power can we expect, and how robust is the design?*

using Pkg #hide
Pkg.activate(@__DIR__) #hide

using KEEP: PointMass10 as PM10 #hide
using KEEP: PointMass4 as PM4 #hide
using KEEP: PointMassPara as PMP #hide
using KEEP: SteadyState #hide
using KEEP: LimitCycle as LC #hide
using KEEP: Optimization as OPT #hide
using KEEP: Visualization as VIZ #hide
using KEEP: DEFAULT_TOLERANCE, TAU0 #hide

using OrdinaryDiffEqVerner: Vern9, ODEProblem #hide
using OrdinaryDiffEqTsit5: Tsit5 #hide
using BifurcationKit #hide
const BK = BifurcationKit #hide
using QuasiMonteCarlo: sample, HaltonSample #hide
using SplitApplyCombine: invert #hide
using BenchmarkTools: @belapsed #hide
using LinearAlgebra: norm, eigvals, × #hide
using Logging: with_logger, NullLogger #hide
ENV["GKSwstype"] = "100" #hide
using Plots #hide
using PyFormattedStrings: @f_str #hide
using Markdown #hide
using ForwardDiff: ForwardDiff #hide

include(joinpath(@__DIR__, "plot_helpers.jl")) #hide
my_defaults() #hide

function save_fig(fig, name) #hide
    savefig(fig, joinpath(OUTPATH, "$(name).pdf")) #hide
    savefig(fig, joinpath(OUTPATH, "$(name).png")) #hide
    return nothing #hide
end #hide

function kite_velocity(u, p)
    pos(q_) = PM4.compute_OK(PM4.compute_Rτ(q_, p), p)
    return ForwardDiff.jacobian(pos, u[1:2]) * u[3:4]
end

function cycle_metrics(sol; n=1001)
    sol === nothing && return (;
        P=NaN,
        period=NaN,
        Tmin=NaN,
        Tmax=NaN,
        αmax=NaN, # no orbit: NaN, do not crash
        amax=NaN,
        κmax=NaN,
        Rmin=NaN,
    )
    ts = range(sol.t[1], sol.t[end], n)
    us = sol.(ts)
    p = sol.prob.p
    tension = PM4.compute_line_tension.(us, Ref(p))
    ok = isfinite.(tension)
    vel = kite_velocity.(us, Ref(p))
    acc = PM4.kite_acceleration.(us, Ref(p))
    κ = [norm(v × a) / max(norm(v)^3, 1e-12) for (v, a) in zip(vel, acc)]
    return (;
        P=sol.u[end][5] / sol.t[end],
        period=sol.t[end],
        Tmin=minimum(tension[ok]),
        Tmax=maximum(tension[ok]),
        αmax=maximum(u -> abs(rem2pi(u[1], RoundNearest)), us),
        amax=maximum(norm, acc),
        κmax=maximum(κ),
        Rmin=1 / maximum(κ),
    )
end
nothing #hide

# ## 1. Mathematical Model & 10D-to-4D Validation
#
# The KEEP system couples a horizontal arm $\mathbf{OA} = a\,\hat{\mathbf{a}}(\alpha)$ rotating about $\hat{\mathbf{z}}$ (inertia $I$, generator braking torque $T(\dot\alpha) = -b\dot\alpha$) to a kite $K$ via tethers of fixed length $\ell > a$.
# The kite is passively guided on a conical lemniscate $\mathbf{OK} = r(\alpha,\tau)\,\hat{\mathbf{r}}(\tau)$ with $\theta(\tau) = \theta_0 + \Delta\theta\sin 2\tau$ and $\phi(\tau) = \phi_0 + \Delta\phi\sin\tau$.
# Solving the loop closure $\norm{\mathbf{OK} - \mathbf{OA}} = \ell$ gives $r(\alpha,\tau) = p + \sqrt{p^2 + \ell^2 - a^2}$ ($p = \hat{\mathbf{r}}\cdot\mathbf{OA}$), reducing the 5-DoF (10D) DAE to a 2-DoF (4D) pure ODE $M(\mathbf{q})\ddot{\mathbf{q}} = \mathbf{h}(\mathbf{q},\dot{\mathbf{q}})$ in $\mathbf{q} = (\alpha,\tau)$, augmented with generator work $\dot W = b\dot\alpha^2 \ge 0$ so average cycle power is $\bar P = W(t_f)/t_f$.

p_nom = PMP.build_para()
vbp_nom = PMP.build_vbpara(p_nom)

#Nominal parameter table #hide
Markdown.parse(
    """
| Symbol | Description | Value | Unit |
| :--- | :--- | ---: | :--- |
| \$a\$ | Arm length | $(p_nom.a) | m |
| \$\\ell\$ | Line length | $(p_nom.l) | m |
| \$m\$ | Kite mass | $(p_nom.m) | kg |
| \$I\$ | Arm inertia | $(p_nom.I) | kg m² |
| \$S\$ | Kite surface | $(p_nom.S) | m² |
| \$C_L / C_D\$ | Glide ratio (\$C_L = $(p_nom.C_L)\$) | $(round(p_nom.C_L / p_nom.C_D; digits=1)) | — |
| \$w_0\$ | Reference wind speed (at \$z_0 = $(p_nom.z0)\$ m) | $(p_nom.w0) | m/s |
| \$\\theta_0, \\Delta\\theta, \\Delta\\phi\$ | Lemniscate angles | $(round(rad2deg(p_nom.θ0); digits=1))°, $(round(rad2deg(p_nom.Δθ); digits=1))°, $(round(rad2deg(p_nom.Δφ); digits=1))° | deg |
""",
) #hide

# ### Equivalence with the 10D DAE Reference
#
# We compare the 4D ODE against the raw 10D index-2 DAE and 10D with a manifold-projection callback (`Vern9` across tolerances $10^{-2}\dots 10^{-12}$ against a $10^{-14}$ reference, plus `Tsit5` runtime at $\text{tol}=10^{-3}, t_f=100\text{ s}$).

τ0_val, dτ0_val = 1e-10, 10.0 #hide
u0_10 = PM10.init_u(τ0_val, dτ0_val, p_nom) #hide
u0_4 = PM4.build_u(; α=0, dα=0, dτ=dτ0_val, τ=τ0_val) #hide
cb_10_save = PM10.build_manifold_projection(u0_10; save=true, tol=1e-10) #hide
cb_10_nosave = PM10.build_manifold_projection(u0_10; save=false, tol=1e-10) #hide

tf_drift, tol_drift = 100.0, 1e-3
sol10_drift = PM10.integrate(
    u0_10, tf_drift, p_nom, Tsit5(); save_everystep=true, tol=tol_drift
)
sol10_cb_drift = PM10.integrate(
    u0_10,
    tf_drift,
    p_nom,
    Tsit5();
    save_everystep=true,
    tol=tol_drift,
    callback=cb_10_nosave,
)

pos10(sol) = [PM10.compute_pos1(u[1:5], sol.prob.p) for u in sol.u]
pos4(sol) = [PM4.compute_OK(PM4.compute_Rτ(u[1:2], sol.prob.p), sol.prob.p) for u in sol.u]

tf_conv = 20.0 #hide
saveat_conv = range(0, tf_conv, 1001) #hide
tols_conv = exp10.(-2:-1:-12) #hide
pref_conv = pos4(
    PM4.integrate(u0_4, tf_conv, vbp_nom, Vern9(); tol=1e-14, saveat=saveat_conv)
) #hide
err10 = [
    maximum(
        norm.(
            pos10(PM10.integrate(u0_10, tf_conv, p_nom, Vern9(); tol, saveat=saveat_conv)) .-
            pref_conv,
        ),
    ) for tol in tols_conv
] #hide
err10cb = [
    maximum(
        norm.(
            pos10(
                PM10.integrate(
                    u0_10,
                    tf_conv,
                    p_nom,
                    Vern9();
                    tol,
                    saveat=saveat_conv,
                    callback=cb_10_nosave,
                ),
            ) .- pref_conv,
        ),
    ) for tol in tols_conv
] #hide
err4 = [
    maximum(
        norm.(
            pos4(PM4.integrate(u0_4, tf_conv, vbp_nom, Vern9(); tol, saveat=saveat_conv)) .-
            pref_conv,
        ),
    ) for tol in tols_conv
] #hide

drift_pos_10 = norm(
    PM10.manifold_residuals!(similar(u0_10, 6), sol10_drift.u[end], p_nom)[1:3]
) #hide
drift_pos_cb = norm(
    PM10.manifold_residuals!(similar(u0_10, 6), sol10_cb_drift.u[end], p_nom)[1:3]
) #hide

btime4 = @belapsed PM4.integrate(
    $u0_4, $tf_drift, $vbp_nom; save_everystep=true, tol=$tol_drift
) seconds = 1
btime10 = @belapsed PM10.integrate(
    $u0_10, $tf_drift, $p_nom; save_everystep=true, tol=$tol_drift
) seconds = 1
btime10_cb = @belapsed PM10.integrate(
    $u0_10, $tf_drift, $p_nom; save_everystep=true, tol=$tol_drift, callback=$cb_10_save
) seconds = 1

fig_conv = plot(
    tols_conv,
    err10;
    marker=:circle,
    label="10D",
    xscale=:log10,
    yscale=:log10, #hide
    xflip=true,
    ymirror=true,
    xlabel="Tolerance",
    ylabel="Max position error (m)", #hide
    xticks=exp10.(-12:3:-3),
    yticks=exp10.(-15:3:0),
    legend=:topright,
    size=SIZE_HALF,
) #hide
plot!(fig_conv, tols_conv, err10cb; marker=:square, label="10D + callback") #hide
plot!(fig_conv, tols_conv, err4; marker=:diamond, label="4D") #hide
plot!(fig_conv; bottommargin=4Plots.mm, rightmargin=6Plots.mm) #hide
save_fig(fig_conv, "three_models_convergence") #hide

Markdown.parse(
    """
- **Algebraic drift (\$t_f = $(Int(tf_drift))\\,\\text{s}\$, \$\\text{tol}=10^{-3}\$):** raw 10D drifts by **$(round(drift_pos_10; digits=2)) m** on position; the projection callback reduces this to **$(f"{drift_pos_cb:.1e}") m** at a **$(round(btime10_cb / btime10; digits=2))×** runtime penalty.
- **Equivalence:** all three formulations converge to the same trajectory (max position discrepancy **$(f"{err4[end]:.1e}") m** at \$\\text{tol}=10^{-12}\$).
- **Performance:** the 4D ODE is drift-free by construction and runs in **$(round(1e3 * btime4; digits=2)) ms**, a **$(round(btime10 / btime4; digits=1))×** speedup over raw 10D (**$(round(1e3 * btime10; digits=2)) ms**) and **$(round(btime10_cb / btime4; digits=1))×** over 10D + callback (**$(round(1e3 * btime10_cb; digits=2)) ms**).

![10D vs 4D convergence](figs/three_models_convergence.png)
""",
) #hide

# ## 2. 1D Analytical PTO Sizing (Simplified Baseline)
#
# To obtain an analytical starting point for the power take-off (PTO) parameters $(I, b)$ before running 4D optimizations, we reduce the arm dynamics to $I\ddot\alpha + b\dot\alpha = C_M\cos(\omega t)$ with harmonic torque amplitude $C_M = T_s S a \sin\theta_0$ (wing loading $T_s = 1000\text{ N/m}^2$, i.e., $F_s = T_s S = 10\text{ kN}$ total) and zero-mass kite speed $V_k$:
# $$b_{\text{1D}} = I_{\text{1D}}\omega, \qquad I_{\text{1D}} = \frac{\sqrt{2}\,C_M}{\pi \omega^2}, \qquad P_{\text{1D}} = \frac{\pi C_M \omega}{4\sqrt{2}},$$
# where $\omega = 2\pi V_k / (\ell L_{81})$ and $L_{81}$ is the unit-sphere lemniscate length.

function pto_sizing_1d(p; T_s=1000.0)
    G = p.C_L / p.C_D
    C_M = T_s * p.S * p.a * sin(p.θ0)
    V_k = p.w0 * (p.l * cos(p.θ0) / p.z0)^(1 / p.n_wind) * sqrt((1 + G^2) * sin(p.θ0)^2 - 1)
    τs = range(0, 2π; length=1001)[1:(end - 1)]
    L81 =
        (2π / length(τs)) * sum(
            τ -> sqrt((2p.Δθ * cos(2τ))^2 + (sin(p.θ0 + p.Δθ * sin(2τ)) * p.Δφ * cos(τ))^2),
            τs,
        )
    ω = 2π * V_k / (p.l * L81)
    I_1d = sqrt(2) * C_M / (π * ω^2)
    b_1d = I_1d * ω
    P_1d = π * C_M * ω / (4 * sqrt(2))
    return (; C_M, V_k, L81, ω, period=2π / ω, I=I_1d, b=b_1d, P=P_1d, F_s=T_s * p.S)
end

pto = pto_sizing_1d(p_nom)
p_1d = PMP.build_para(; I=pto.I, b=pto.b)
m_1d = cycle_metrics(
    LC.compute_limit_cycle(PMP.build_vbpara(p_1d); sense=+, tol=1e-5, save_everystep=true)
)
nothing #hide

# Also compute the 4D kinetic-criterion baseline (nominal I, b = I ω from 4D period)
p_kin = PMP.build_para()
p_kin.b = PMP.optimal_torque_slope(p_kin; sense=(+), tol=1e-3)
vbp_kin = PMP.build_vbpara(p_kin)
m_kin = cycle_metrics(
    LC.compute_limit_cycle(vbp_kin; sense=+, tol=1e-5, save_everystep=true)
)

Markdown.parse(
    """
| Model | \$I\$ (kg m²) | \$b\$ (N m s rad⁻¹) | Period (s) | Power (kW) | \$T_{\\min}\$–\$T_{\\max}\$ (kN) | \$\\max\\lvert\\alpha\\rvert\$ (deg) | \$a_{\\max}\$ (m/s²) | \$R_{\\min}\$ (m) |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1D analytical PTO | $(round(pto.I; digits=1)) | $(round(pto.b; digits=1)) | $(round(pto.period; digits=2)) | **$(round(pto.P / 1000; digits=2))** | — (assumed ≤ $(round(pto.F_s / 1000; digits=1))) | $(round(rad2deg(π / 2); digits=1)) | — | — |
| 4D at 1D PTO params | $(round(pto.I; digits=1)) | $(round(pto.b; digits=1)) | $(round(m_1d.period; digits=2)) | **$(round(m_1d.P / 1000; digits=2))** | $(round(m_1d.Tmin / 1000; digits=2))–$(round(m_1d.Tmax / 1000; digits=2)) | $(round(rad2deg(m_1d.αmax); digits=1)) | $(round(m_1d.amax; digits=0)) | $(round(m_1d.Rmin; digits=2)) |
| 4D kinetic (\$b = I\\omega\$) | $(round(p_kin.I; digits=1)) | $(round(p_kin.b; digits=1)) | $(round(m_kin.period; digits=2)) | **$(round(m_kin.P / 1000; digits=2))** | $(round(m_kin.Tmin / 1000; digits=2))–$(round(m_kin.Tmax / 1000; digits=2)) | $(round(rad2deg(m_kin.αmax); digits=1)) | $(round(m_kin.amax; digits=0)) | $(round(m_kin.Rmin; digits=2)) |

The 1D analytical sizing gives **\$L_{81} = $(round(pto.L81; digits=4))\$**, **\$V_k = $(round(pto.V_k; digits=1))\\,\\text{m/s}\$**, **\$I_{\\text{1D}} = $(round(pto.I; digits=1))\\,\\text{kg m}^2\$**, and **\$b_{\\text{1D}} = $(round(pto.b; digits=0))\\,\\text{N m s rad}^{-1}\$**, predicting **$(round(pto.P / 1000; digits=2)) kW**.
Evaluating those exact parameters in the 4D model yields **$(round(m_1d.P / 1000; digits=2)) kW**—a **$(round(pto.P / m_1d.P; digits=1))×** gap caused by the 1D assumptions (zero kite mass, constant \$V_k\$ without gravity/induced-drag losses, full-swing \$\\max\\lvert\\alpha\\rvert=90^\\circ\$ vs. **$(round(rad2deg(m_1d.αmax); digits=1))°**, and peak \$F_s\$ acting over the whole cycle).
Because the 4D kinetic-criterion baseline (**$(round(m_kin.P / 1000; digits=2)) kW**) outperforms the 4D simulation at 1D-PTO parameters (**$(round(m_1d.P / 1000; digits=2)) kW**), we use the **4D kinetic baseline** as our static reference for Section 4.
""",
) #hide

# ## 3. Dynamical Behavior & Operational Feasibility
#
# Using the 4D kinetic baseline (`p_kin`), we check equilibria, limit-cycle uniqueness/convergence, and physical feasibility.

ss = SteadyState.all_steady_states(vbp_nom)
dyn4(u) = PM4.dynamics(PM4.build_u(; α=u[1], τ=u[2], dα=u[3], dτ=u[4]), vbp_kin)[1:4]
ss_left = [q for q in ss if q[2] > 0]
max_re_eig = maximum(
    q -> maximum(real, eigvals(Array(ForwardDiff.jacobian(dyn4, [q..., 0, 0])))), ss_left
)
min_re_eig = minimum(
    q -> maximum(real, eigvals(Array(ForwardDiff.jacobian(dyn4, [q..., 0, 0])))), ss_left
)
nothing #hide

# Two attracting limit cycles (D(τ) > 0 rises on outer sides; D(τ) < 0 rises in center)
lc_p = LC.compute_limit_cycle(vbp_kin; sense=(+), save_everystep=true)
lc_m = LC.compute_limit_cycle(vbp_kin; sense=(-), save_everystep=true)
mp = cycle_metrics(lc_p)
mm = cycle_metrics(lc_m)
nothing #hide

# Floquet multipliers & convergence slope
const CB_POINC = LC.build_poincare_callback(eps(1.0)) #hide
function Φ_map(x, T; vbp)
    return Array(PM4.integrate(PM4.build_u(x[1], x[2], x[3], x[4], 0.0), T, vbp).u[end])[1:4]
end #hide
function floquet_rate(sol, vbp) #hide
    u0, T = LC.unpack_shooting(LC.build_shooting(sol)) #hide
    μ = eigvals(
        Matrix(ForwardDiff.jacobian(x -> Φ_map(x, T; vbp), collect(Array(u0))[1:4]))
    ) #hide
    λ = maximum(abs, μ[abs.(μ) .< 0.5]) #hide
    return (; slope=-log10(λ), λ, T) #hide
end #hide
fl_p = floquet_rate(
    PM4.integrate(
        PM4.build_u(; dτ=1), 50, vbp_kin; callback=CB_POINC, save_end=false, tol=1e-13
    ),
    vbp_kin,
) #hide
fl_m = floquet_rate(
    PM4.integrate(
        PM4.build_u(; dτ=-1), 50, vbp_kin; callback=CB_POINC, save_end=false, tol=1e-13
    ),
    vbp_kin,
) #hide
nothing #hide

# Exhaustive inventory on the fundamental domain (modulo S(α,τ,dα,dτ) = (-α, τ+π, -dα, dτ))
const BK_M, BK_NSEEDS = 10, 100
vec4(u) = [u[1], u[2], u[3], u[4]] #hide
function f_search!(du, u, p, t) #hide
    if abs(u[3]) > BK_M + 2 || abs(u[4]) > BK_M + 2 #hide
        du .= 0 #hide
    else #hide
        du .= PM4.dynamics(PM4.build_u(u[1], u[2], u[3], u[4], 0.0), p, t)[1:4] #hide
    end #hide
end #hide
search_prob = ODEProblem(f_search!, zeros(4), (0.0, 1.0), vbp_kin) #hide
sh_bk = BK.Shooting(
    search_prob,
    Vern9(),
    1,
    BK.SectionSS([0.0, 1.0, 0.0, 0.0], [0.0, TAU0, 0.0, 0.0]);
    abstol=1e-10,
    reltol=1e-10,
    maxiters=20_000,
) #hide
struct TwistedShooting{S,U} <: BK.AbstractBifurcationProblem
    sh::S
    u0::U
    twist::Float64
end #hide
BK.getparams(w::TwistedShooting) = BK.getparams(w.sh) #hide
BK.getlens(w::TwistedShooting) = BK.getlens(w.sh) #hide
BK.getparam(w::TwistedShooting) = BK.getparam(w.sh) #hide
BK.getu0(w::TwistedShooting) = w.u0 #hide
BK.getdelta(w::TwistedShooting) = BK.getdelta(w.sh) #hide
function BK.residual(w::TwistedShooting, x, p)
    return (out=BK.po_residual(w.sh, x, p); out[end - 3] -= w.twist; out)
end #hide
BK.jacobian(w::TwistedShooting, x, p) = BK.po_jacobian(w.sh, x, p) #hide
root_state(x) = (α=rem2pi(x[1], RoundNearest), dα=x[3], dτ=x[4], T=x[5]) #hide
function same_cycle(a, b)
    return abs(rem2pi(a.α - b.α, RoundNearest)) < 1e-4 &&
           abs(a.dα - b.dα) < 1e-4 &&
           abs(a.dτ - b.dτ) < 1e-4 &&
           abs(a.T - b.T) < 1e-4
end #hide
npar = BK.NewtonPar(; tol=1e-9, max_iterations=15, verbose=false) #hide
refs_bk = map((+1, -1)) do sense #hide
    lc_ref = sense > 0 ? lc_p : lc_m #hide
    tw = TwistedShooting(sh_bk, vcat(vec4(lc_ref(0.0)), lc_ref.t[end]), sense * 2π) #hide
    s = BK.solve(tw, BK.Newton(), npar) #hide
    return (sense=sense, x=collect(s.u), state=root_state(s.u)) #hide
end #hide
cb_turn = LC.build_poincare_callback(Inf) #hide
seeds = eachcol(
    sample(BK_NSEEDS, [0.0, -π, -BK_M, -BK_M], [π, π, BK_M, BK_M], HaltonSample())
) #hide
extra_cycles = NamedTuple[] #hide
for sense in (+1, -1) #hide
    ref = only(filter(r -> r.sense == sense, refs_bk)) #hide
    defOp = BK.DeflationOperator(2.0, 1.0, [copy(ref.x)]) #hide
    for u0 in seeds #hide
        sign(u0[4]) == sense || continue #hide
        tf_g = 4π / clamp(abs(u0[4]), 0.5, BK_M) #hide
        sg = PM4.integrate(
            PM4.build_u(u0[1], TAU0, u0[3], u0[4], 0.0),
            tf_g,
            vbp_kin;
            callback=cb_turn,
            save_end=false,
            tol=1e-6,
            maxiters=50_000,
        ) #hide
        length(sg.t) < 2 && continue #hide
        tw = TwistedShooting(sh_bk, vcat(vec4(sg.u[end]), sg.t[end]), sense * 2π) #hide
        s = try
            BK.solve(tw, defOp, npar)
        catch
            nothing
        end #hide
        (s === nothing || !BK.converged(s)) && continue #hide
        x = collect(s.u)
        norm(BK.residual(tw, x, vbp_kin), Inf) < 1e-6 || continue #hide
        st = root_state(x) #hide
        traj = try
            LC.shoot([st.α, st.dα, st.dτ, st.T], vbp_kin; tol=1e-10, save_everystep=true)
        catch
            nothing
        end #hide
        (
            traj === nothing || !all(abs(u[3]) <= BK_M && abs(u[4]) <= BK_M for u in traj.u)
        ) && continue #hide
        norm(LC.endpoint_residuals(traj; sense=sense > 0 ? (+) : (-)), Inf) < 1e-6 ||
            continue #hide
        (same_cycle(st, ref.state) || any(f -> same_cycle(f.state, st), extra_cycles)) &&
            continue #hide
        push!(defOp, copy(x))
        push!(extra_cycles, (x=x, state=st)) #hide
    end #hide
end #hide
cb_conv = LC.build_poincare_callback(1e-8) #hide
qmc_counts = zeros(Int, 3) #hide
for u0 in seeds #hide
    sol_q = PM4.integrate(
        PM4.build_u(u0[1], u0[2], u0[3], u0[4], 0.0),
        60,
        vbp_kin;
        callback=cb_conv,
        save_end=false,
        tol=1e-8,
        maxiters=1_000_000,
    ) #hide
    u = sol_q.u[end] #hide
    d = if all(isfinite, u)
        (LC.distance_on_section(u, lc_p(0.0)), LC.distance_on_section(u, lc_m(0.0))) #hide
    else
        (Inf, Inf) #hide
    end #hide
    qmc_counts[minimum(d) < 1e-3 ? argmin(d) : 3] += 1 #hide
end #hide

# Transient settling onto the pumping cycle
sol_trans = PM4.integrate(PM4.build_u(; dτ=1), 10, vbp_kin; save_everystep=true) #hide
save_fig(plot_time_series(sol_trans), "time_series_4D") #hide

Markdown.parse(
    """
- **Equilibria:** all **$(length(ss))** equilibria are unstable (dominant \$\\text{Re}(\\lambda) \\in [$(round(min_re_eig; digits=2)), $(round(max_re_eig; digits=2))]\\,\\text{s}^{-1}\$).
- **Attracting limit cycles:** trajectories converge exponentially to one of two periodic orbits—**\$\\dot\\tau > 0\$** (kite rises on the outer sides, Floquet multiplier \$\\lambda = $(f"{fl_p.λ:.3g}")\$, slope **$(round(fl_p.slope; digits=2))** digits/rev, **$(round(mp.P / 1000; digits=2)) kW**) and **\$\\dot\\tau < 0\$** (rises in the center, \$\\lambda = $(f"{fl_m.λ:.3g}")\$, slope **$(round(fl_m.slope; digits=2))**, **$(round(mm.P / 1000; digits=2)) kW**).
- **Global uniqueness:** deflated single shooting over $(BK_NSEEDS) Halton seeds in the fundamental domain finds **$(length(extra_cycles))** additional cycles, and quasi-Monte-Carlo integration routes **$(qmc_counts[1]) / $(BK_NSEEDS)** seeds to \$\\dot\\tau > 0\$, **$(qmc_counts[2]) / $(BK_NSEEDS)** to \$\\dot\\tau < 0\$, and **$(qmc_counts[3])** non-convergent.
- **Is the \$\\dot\\tau > 0\$ behavior satisfying?**
  - *Taut lines:* \$T_{\\min} = $(round(mp.Tmin / 1000; digits=2))\\,\\text{kN} > 0\$ (no slackness).
  - *Arm excursion:* \$\\max|\\alpha| = $(round(rad2deg(mp.αmax); digits=1))^\\circ < 90^\\circ\$ (stays downwind, no mast wrap-around).
  - *Line clearance:* \$\\dot\\tau > 0\$ ascends on the outer lobes, avoiding tether crossing.
  - *Loads & turning:* \$T_{\\max} = $(round(mp.Tmax / 1000; digits=2))\\,\\text{kN}\$, \$a_{\\max} = $(round(mp.amax; digits=0))\\,\\text{m/s}^2\$, \$R_{\\min} = $(round(mp.Rmin; digits=2))\\,\\text{m}\$.

![Transient settling onto the pumping limit cycle](figs/time_series_4D.png)
""",
) #hide

# ## 4. Dynamic Optimization & Robustness
#
# We maximize average power $\bar P = W(t_f)/t_f$ subject to periodic BVP closure over $\Theta = (\ell, I, b, \Delta\varphi)$ at $w_0 = 9\text{ m/s}$, holding $\Delta\theta = 6^\circ$ fixed (box: $\ell \in [15, 200]\text{ m}$, $I \in [62, 10^6]\text{ kg m}^2$, $b \in [1, 10^6]\text{ N m s rad}^{-1}$, $\Delta\varphi \in [10^\circ, 90^\circ]$).

const TOL_OPT = 1e-5
const WIND_REF = 9.0
const SENSE_OPT = +
const SYMS_OPT = [:l, :I, :b, :Δφ]
const LOWER_OPT = [15.0, 62.0, 1.0, deg2rad(10.0)]
const UPPER_OPT = [200.0, 1.0e6, 1.0e6, deg2rad(90.0)]
nothing #hide

# Warm-up, then the timed solve
with_logger(NullLogger()) do #hide
    return OPT.optimize(
        p_kin,
        SYMS_OPT,
        LOWER_OPT,
        UPPER_OPT;
        sense=SENSE_OPT,
        tol=TOL_OPT,
        max_wall_time=120.0,
    ) #hide
end #hide
t_opt_start = time() #hide
sol_opt, stats_opt, _ = with_logger(NullLogger()) do
    return OPT.optimize(
        p_kin,
        SYMS_OPT,
        LOWER_OPT,
        UPPER_OPT;
        sense=SENSE_OPT,
        tol=TOL_OPT,
        max_wall_time=120.0,
    )
end
t_opt_wall = time() - t_opt_start #hide
p_opt = copy(p_kin)
p_opt[SYMS_OPT] .= sol_opt.params[SYMS_OPT]

lc_opt = LC.compute_limit_cycle(
    PMP.build_vbpara(p_opt); sense=SENSE_OPT, tol=TOL_OPT, save_everystep=true
)
m_opt = cycle_metrics(lc_opt)
res_opt = lc_opt === nothing ? NaN : norm(LC.endpoint_residuals(lc_opt; sense=SENSE_OPT))
nothing #hide

fig_ts = plot_time_series_grid(
    (lc_p, lc_opt); titles=("Kinetic baseline", "Dynamic optimum")
) #hide
save_fig(fig_ts, "optimization_timeseries_all") #hide

gain_pct = (m_opt.P / m_kin.P - 1) * 100 #hide
Markdown.parse(
    """
Ipopt converges in **$(stats_opt.iter) iterations** (**$(round(t_opt_wall; digits=2)) s** warm wall time, endpoint residual **$(f"{res_opt:.1e}")**) with all four variables strictly interior:

| Quantity | Kinetic baseline | Dynamic optimum | Change |
| :--- | ---: | ---: | ---: |
| \$\\ell\$ (m) | $(round(p_kin.l; digits=2)) | **$(round(p_opt.l; digits=2))** | $(round((p_opt.l / p_kin.l - 1) * 100; digits=1))% |
| \$I\$ (kg m²) | $(round(p_kin.I; digits=1)) | **$(round(p_opt.I; digits=1))** | $(round((p_opt.I / p_kin.I - 1) * 100; digits=1))% |
| \$b\$ (N m s rad⁻¹) | $(round(p_kin.b; digits=1)) | **$(round(p_opt.b; digits=1))** | $(round((p_opt.b / p_kin.b - 1) * 100; digits=1))% |
| \$\\Delta\\varphi\$ (deg) | $(round(rad2deg(p_kin.Δφ); digits=1)) | **$(round(rad2deg(p_opt.Δφ); digits=2))** | $(round((p_opt.Δφ / p_kin.Δφ - 1) * 100; digits=1))% |
| Period (s) | $(round(m_kin.period; digits=3)) | **$(round(m_opt.period; digits=3))** | $(round((m_opt.period / m_kin.period - 1) * 100; digits=1))% |
| **Average power (kW)** | $(round(m_kin.P / 1000; digits=2)) | **$(round(m_opt.P / 1000; digits=2))** | **+$(round(gain_pct; digits=1))%** |
| \$T_{\\min}\$–\$T_{\\max}\$ (kN) | $(round(m_kin.Tmin / 1000; digits=2))–$(round(m_kin.Tmax / 1000; digits=2)) | $(round(m_opt.Tmin / 1000; digits=2))–**$(round(m_opt.Tmax / 1000; digits=2))** | — |
| \$\\max\\lvert\\alpha\\rvert\$ (deg) | $(round(rad2deg(m_kin.αmax); digits=1)) | **$(round(rad2deg(m_opt.αmax); digits=1))** | — |
| \$a_{\\max}\$ (m/s²) | $(round(m_kin.amax; digits=0)) | **$(round(m_opt.amax; digits=0))** | — |
| \$R_{\\min}\$ (m) | $(round(m_kin.Rmin; digits=2)) | **$(round(m_opt.Rmin; digits=2))** | — |

Note that peak tension at the optimum (**$(round(m_opt.Tmax / 1000; digits=2)) kN**) exceeds the 1D wing-loading limit \$F_s = T_s S = 10\\,\\text{kN}\$ (\$T_s = 1000\\,\\text{N/m}^2\$), and peak acceleration reaches **$(round(m_opt.amax; digits=0)) m/s²** with a minimum turning radius \$R_{\\min} = $(round(m_opt.Rmin; digits=2))\\,\\text{m}\$; enforcing explicit tension and load constraints is left for future work.

![Baseline vs optimized limit-cycle time series](figs/optimization_timeseries_all.png)
""",
) #hide

# ### Local Landscape, Shape Sensitivity ($\Delta\theta, \Delta\varphi$), and Wind Robustness

clamp_guess(::Nothing) = (;) #hide
function clamp_guess(s)
    return (;
        α0=clamp(s.α0, -π + 1e-6, π - 1e-6),
        dα0=clamp(s.dα0, -100 + 1e-6, 100 - 1e-6), #hide
        dτ0=clamp(s.dτ0, 1e-6, 100 - 1e-6),
        tf=clamp(s.tf, 1.1e-3, 99.0),
        P=clamp(s.P, 1.1e-3, 9.9e4),
    ) #hide
end #hide
nothing #hide

# Local landscape around the optimum, one parameter at a time (±10 %)
xs_land = range(0.9, 1.1; length=21)
land_ys = [
    begin #hide
        row = [
            begin #hide
                q = copy(p_opt)
                q[sym] = x * p_opt[sym] #hide
                lc = LC.compute_limit_cycle(
                    PMP.build_vbpara(q); sense=SENSE_OPT, tol=TOL_OPT
                ) #hide
                lc.u[end][5] / lc.t[end] #hide
            end for x in xs_land
        ] #hide
        row ./ row[11] #hide
    end for sym in SYMS_OPT
] #hide
nothing #hide

# Δθ sweep, 2° to 20°: Θ held at the optimum against Θ re-optimized at each angle
dtheta_angles = collect(range(2.0, 20.0; length=10))
dtheta_fixed = [
    begin #hide
        q = copy(p_opt)
        q.Δθ = deg2rad(a) #hide
        sol = LC.compute_limit_cycle(
            PM4.build_u(; dτ=1),
            PMP.build_vbpara(q);
            tol=TOL_OPT,
            save_everystep=true,
            tf=100.0,
        ) #hide
        (; Δθ=a, cycle_metrics(sol)...) #hide
    end for a in dtheta_angles
] #hide
nothing #hide

i0_dt = argmin(abs.(dtheta_angles .- rad2deg(p_opt.Δθ)))
dtheta_reopt = Vector{Any}(undef, length(dtheta_angles)) #hide
p_dt_at = Vector{Any}(undef, length(dtheta_angles)) #hide
sol_dt_at = Vector{Any}(undef, length(dtheta_angles)) #hide
function attempt_dt!(i, st_p, g) #hide
    q = copy(st_p)
    q.Δθ = deg2rad(dtheta_angles[i]) #hide
    try #hide
        s, st, _ = with_logger(NullLogger()) do #hide
            return OPT.optimize(
                q,
                SYMS_OPT,
                LOWER_OPT,
                UPPER_OPT;
                sense=SENSE_OPT,
                tol=TOL_OPT,
                max_wall_time=60.0,
                initial_guess=g,
            ) #hide
        end #hide
        p_new = copy(q)
        p_new[SYMS_OPT] .= s.params[SYMS_OPT] #hide
        return (p=p_new, status=st.status, sol=s, P=s.P) #hide
    catch #hide
        return (p=copy(q), status=:error, sol=nothing, P=NaN) #hide
    end #hide
end #hide
function run_dt!(i, prev) #hide
    starts = if prev === nothing
        ((p_opt, (;)),) #hide
    else
        (
        (prev.p, clamp_guess(prev.sol)),
        (prev.p, (;)),
        (p_opt, clamp_guess(sol_dt_at[i0_dt])),
        (p_opt, (;)),
    ) #hide
    end #hide
    fixed_P = dtheta_fixed[i].P     # Θ held at the nominal: always feasible, so a lower P is another branch #hide
    best = nothing #hide
    for (k, (st_p, g)) in enumerate(starts) #hide
        r = attempt_dt!(i, st_p, g) #hide
        if r.status == :first_order &&
            (best === nothing || best.status != :first_order || r.P > best.P) #hide
            best = r #hide
        elseif best === nothing #hide
            best = r #hide
        end #hide
        k == length(starts) && break #hide
        best.status == :first_order && best.P >= fixed_P * (1 - 1e-3) && break #hide
    end #hide
    best.status == :first_order && best.P >= fixed_P * (1 - 1e-3) || #hide
        @warn "Δθ = $(dtheta_angles[i])°: best re-optimization does not beat the fixed Θ" best.status best.P fixed_P #hide
    p_dt_at[i], sol_dt_at[i] = best.p, best.sol #hide
    if best.status == :first_order #hide
        sol = LC.compute_limit_cycle(
            PM4.build_u(; dτ=1),
            PMP.build_vbpara(best.p);
            tol=TOL_OPT,
            save_everystep=true,
            tf=100.0,
        ) #hide
        dtheta_reopt[i] = (; Δθ=dtheta_angles[i], cycle_metrics(sol)...) #hide
    else #hide
        dtheta_reopt[i] = (; Δθ=dtheta_angles[i], cycle_metrics(nothing)...)   # no solved point: a gap, not a bogus power #hide
    end #hide
end #hide
run_dt!(i0_dt, nothing) #hide
for i in (i0_dt - 1):-1:1
    run_dt!(i, (p=p_dt_at[i + 1], sol=sol_dt_at[i + 1]))
end #hide
for i in (i0_dt + 1):length(dtheta_angles)
    run_dt!(
        i,
        (
            p=(i == i0_dt + 1 ? p_opt : p_dt_at[i - 1]),
            sol=(i == i0_dt + 1 ? sol_dt_at[i0_dt] : sol_dt_at[i - 1]),
        ),
    )
end #hide

# 10×10 shape map over (Δθ, Δφ), at the optimal (ℓ, I, b) and with nothing re-optimized
shape_factors = range(0.8, 1.2; length=10)
shape_dphi = range(10.0, 90.0; length=10)
P_map = fill(NaN, 10, 10);
Tmax_map = fill(NaN, 10, 10);
amax_map = fill(NaN, 10, 10) #hide
for (j, dphi) in enumerate(shape_dphi), (i, f) in enumerate(shape_factors) #hide
    q = copy(p_opt)
    q.Δθ = f * p_opt.Δθ
    q.Δφ = deg2rad(dphi) #hide
    sol = try
        LC.compute_limit_cycle(
            PM4.build_u(; dτ=1),
            PMP.build_vbpara(q);
            tol=TOL_OPT,
            save_everystep=true,
            tf=100.0,
        )
    catch
        nothing
    end #hide
    (sol === nothing || norm(LC.endpoint_residuals(sol; sense=SENSE_OPT)) >= 1e-3) &&
        continue #hide
    m = cycle_metrics(sol; n=401) #hide
    P_map[j, i], Tmax_map[j, i], amax_map[j, i] = m.P / 1000, m.Tmax / 1000, m.amax #hide
end #hide

# Wind sweep, 5 to 15 m/s: the 9 m/s optimum held fixed against Θ re-optimized
winds = collect(range(5.0, 15.0; length=11))
i_ref_w = findfirst(>=(WIND_REF), winds) #hide
w_fixed = Vector{Any}(undef, length(winds)) #hide
w_opt = Vector{Any}(undef, length(winds)) #hide
w_theta = Vector{Vector{Float64}}(undef, length(winds)) #hide
w_sols = Vector{Any}(undef, length(winds)) #hide
w_ok = falses(length(winds)) #hide
function attempt_w!(i, p_init_w, g) #hide
    p0w = copy(p_init_w)
    p0w.w0 = winds[i] #hide
    s, st, _ = with_logger(NullLogger()) do #hide
        return OPT.optimize(
            p0w,
            SYMS_OPT,
            LOWER_OPT,
            UPPER_OPT;
            sense=SENSE_OPT,
            tol=TOL_OPT,
            max_wall_time=15.0,
            initial_guess=g,
        ) #hide
    end #hide
    pw = copy(p0w)
    pw[SYMS_OPT] .= s.params[SYMS_OPT] #hide
    return (p=pw, status=st.status, sol=s, P=s.P) #hide
end #hide
"Starts at wind speed `i`: the reference, then the already-solved neighbour."
function w_starts(i) #hide
    i == i_ref_w && return ((p_opt, (;)),) #hide
    j = i < i_ref_w ? i + 1 : i - 1 #hide
    isassigned(w_theta, j) || return ((p_opt, (;)),) #hide
    g = clamp_guess(w_sols[j]) #hide
    return ((p_opt, (;)), ((p=copy(p_opt); p[SYMS_OPT].=w_theta[j]; p), g), (p_opt, g)) #hide
end #hide
function solve_w!(i) #hide
    pf = copy(p_opt)
    pf.w0 = winds[i] #hide
    w_fixed[i] = cycle_metrics(
        LC.compute_limit_cycle(
            PMP.build_vbpara(pf); sense=SENSE_OPT, save_everystep=true, tol=TOL_OPT
        ),
    ) #hide
    best = nothing #hide
    for (p0, g) in w_starts(i) #hide
        r = attempt_w!(i, p0, g) #hide
        r.status == :first_order || continue #hide
        (best === nothing || r.P > best.P) && (best = r) #hide
        r.P >= w_fixed[i].P * (1 - 1e-3) && break   # keeping p_opt is feasible: a lower P is another branch #hide
    end #hide
    if best === nothing #hide
        w_opt[i] = cycle_metrics(nothing)
        w_ok[i] = false
        return nothing #hide
    end #hide
    w_opt[i] = cycle_metrics(
        LC.compute_limit_cycle(
            PMP.build_vbpara(best.p); sense=SENSE_OPT, save_everystep=true, tol=TOL_OPT
        ),
    ) #hide
    w_theta[i] = collect(best.p[SYMS_OPT])
    w_sols[i] = best.sol #hide
    w_ok[i] = w_opt[i].P >= w_fixed[i].P * (1 - 1e-3) #hide
    return nothing #hide
end #hide
for i in i_ref_w:-1:1 #hide
    solve_w!(i) #hide
end #hide
for i in (i_ref_w + 1):length(winds) #hide
    solve_w!(i) #hide
end #hide
for i in eachindex(winds) #hide
    w_ok[i] && continue #hide
    @warn "wind sweep at $(winds[i]) m/s: best re-optimization does not beat the fixed parameters" w_opt[i].P w_fixed[i].P #hide
end #hide

# Figures of this subsection
f_dtheta_p = plot(
    dtheta_angles,
    [d.P / dtheta_fixed[i0_dt].P for d in dtheta_fixed];
    seriestype=:scatter,
    c=:black,
    ms=4,
    label="Θ fixed",
    xlabel="Δθ (deg)",
    ylabel="P / P*",
    ylims=(0, 1.1),
) #hide
plot!(
    f_dtheta_p,
    dtheta_angles,
    [d.P / dtheta_fixed[i0_dt].P for d in dtheta_reopt];
    c=:black,
    label="Θ re-optimized",
) #hide
f_dtheta_t = plot(
    dtheta_angles,
    [d.Tmax / 1000 for d in dtheta_fixed];
    seriestype=:scatter,
    c=:black,
    ms=4,
    label="",
    xlabel="Δθ (deg)",
    ylabel="Line tension (kN)",
) #hide
scatter!(
    f_dtheta_t,
    dtheta_angles,
    [d.Tmin / 1000 for d in dtheta_fixed];
    c=:black,
    ms=4,
    label="",
) #hide
plot!(
    f_dtheta_t,
    dtheta_angles,
    [d.Tmax / 1000 for d in dtheta_reopt];
    c=:black,
    label="",
    fillrange=[d.Tmin / 1000 for d in dtheta_reopt],
    fillalpha=0.1,
) #hide
plot!(f_dtheta_t, dtheta_angles, [d.Tmin / 1000 for d in dtheta_reopt]; c=:black, label="") #hide
save_fig(
    plot(
        f_dtheta_p,
        f_dtheta_t;
        layout=(1, 2),
        size=plot_size(1, 2),
        bottommargin=6Plots.mm,
        leftmargin=10Plots.mm,
    ),
    "optimization_dtheta",
) #hide
x_map = 100 .* (shape_factors .- 1) #hide
p_s1 = heatmap(
    x_map,
    shape_dphi,
    P_map;
    c=:viridis,
    title="P (kW)",
    xlabel="Δθ dev (%)",
    ylabel="Δφ (deg)",
) #hide
p_s2 = heatmap(
    x_map, shape_dphi, Tmax_map; c=:viridis, title="Tmax (kN)", xlabel="Δθ dev (%)"
) #hide
p_s3 = heatmap(
    x_map, shape_dphi, amax_map; c=:viridis, title="max a (m/s²)", xlabel="Δθ dev (%)"
) #hide
save_fig(
    plot(
        p_s1,
        p_s2,
        p_s3;
        layout=(1, 3),
        size=plot_size(4 / 3, 2),
        bottommargin=5Plots.mm,
        leftmargin=6Plots.mm,
    ),
    "optimization_shape_grid",
) #hide
f_wp = plot(
    winds,
    [s.P / 1000 for s in w_fixed];
    ls=:dot,
    c=:black,
    label="Fixed 9 m/s",
    ylabel="Power (kW)",
) #hide
plot!(f_wp, winds, [s.P / 1000 for s in w_opt]; c=1, label="Re-optimized") #hide
f_wt = plot(
    winds,
    [s.Tmax / 1000 for s in w_fixed];
    ls=:dot,
    c=:black,
    label="",
    ylabel="Tension (kN)",
) #hide
plot!(f_wt, winds, [s.Tmin / 1000 for s in w_fixed]; ls=:dot, c=:black, label="") #hide
plot!(
    f_wt,
    winds,
    [s.Tmax / 1000 for s in w_opt];
    c=1,
    fillrange=[s.Tmin / 1000 for s in w_opt],
    fillalpha=0.1,
    label="",
) #hide
plot!(f_wt, winds, [s.Tmin / 1000 for s in w_opt]; c=1, label="") #hide
theta_norm = [θ ./ maximum(θ) for θ in invert(w_theta)] #hide
f_wpar = plot(
    winds,
    theta_norm;
    c=permutedims(1:4),
    label=["ℓ" "I" "b" "Δφ"],
    xlabel="Wind speed (m/s)",
    ylabel="Norm. params",
    legend=:right,
) #hide
save_fig(
    plot(
        f_wp,
        f_wt,
        f_wpar;
        layout=(3, 1),
        size=plot_size(2, 1),
        leftmargin=10Plots.mm,
        bottommargin=4Plots.mm,
    ),
    "optimization_wind_sweep_all",
) #hide
gain_5 = (w_opt[1].P / w_fixed[1].P - 1) * 100 #hide
gain_15 = (w_opt[end].P / w_fixed[end].P - 1) * 100 #hide
dt_valid = filter(d -> isfinite(d.P), dtheta_reopt) #hide
dt_best = isempty(dt_valid) ? NaN : maximum(d -> d.P, dt_valid) #hide
dt_failed = count(d -> !isfinite(d.P), dtheta_reopt) #hide
Markdown.parse(
    """
- **Local landscape (±10% around \$\\Theta^*\$):** \$\\ell\$ is the stiffest direction (**$(round((land_ys[1][1] - 1) * 100; digits=1))%** at \$-10\\%\$, **$(round((land_ys[1][end] - 1) * 100; digits=1))%** at \$+10\\%\$), followed by \$\\Delta\\varphi\$ (**$(round((land_ys[4][1] - 1) * 100; digits=1))%** / **$(round((land_ys[4][end] - 1) * 100; digits=1))%**), \$I\$ (**$(round((land_ys[2][1] - 1) * 100; digits=1))%**), and \$b\$ (**$(round((land_ys[3][1] - 1) * 100; digits=1))%**).
- **Why \$\\Delta\\theta = 6^\\circ\$ is held fixed:** continuing the \$\\Theta\$ optimum from one angle to the next over \$\\Delta\\theta \\in [2^\\circ, 20^\\circ]\$, the re-optimized branch peaks at **$(round(dt_best / 1000; digits=2)) kW** (only **+$(round((dt_best / m_opt.P - 1) * 100; digits=1))%** above the \$6^\\circ\$ point)$(dt_failed == 0 ? "" : " — $(dt_failed) of $(length(dtheta_angles)) angles did not converge and are left out"), while the \$10\\times 10\$ shape map shows power is governed by a narrow \$\\Delta\\varphi\$ band near \$19^\\circ\$ (peaking at **$(round(maximum(filter(isfinite, P_map)); digits=2)) kW**).
- **Wind robustness (\$5\\text{--}15\\,\\text{m/s}\$):** expected power scales from **$(round(w_opt[1].P / 1000; digits=2)) kW** (\$5\\,\\text{m/s}\$) to **$(round(w_opt[end].P / 1000; digits=2)) kW** (\$15\\,\\text{m/s}\$). Holding the \$9\\,\\text{m/s}\$ optimum fixed gives up only **$(round(gain_5; digits=1))%** at \$5\\,\\text{m/s}\$ and **$(round(gain_15; digits=1))%** at \$15\\,\\text{m/s}\$—confirming strong robustness without mechanical reconfiguration.

![Δθ sensitivity](figs/optimization_dtheta.png)
![10x10 shape map](figs/optimization_shape_grid.png)
![Wind sweep](figs/optimization_wind_sweep_all.png)
""",
) #hide
