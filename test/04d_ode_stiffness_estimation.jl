# cf. DiffEqDevTools.analyticless_test_convergence

using LinearAlgebra: norm
using StaticArrays
using Plots
using Test
using StrFormat

using KEEP: TAU0
using KEEP.PointMass4: dynamics, integrate
import KEEP.PointMass4 as PM4
using KEEP.PointMassPara
using KEEP.LimitCycle

using KEEP.TorqueFunction

## Calcul d'une solution avec un grand tf pour analyser la convergence vers le cycle limite
vbp = normalize_vbpara(build_vbpara())
tf = 14  # Approximate duration of a limite cycle in normalized time

α0, dα0, dτ0, W0 = 0., 0., 0., 0.
u0 = SA[α0, TAU0, dα0, dτ0, W0]
sol = PM4.integrate(u0, tf, vbp)
Nmax = 30_000
Nmin = 100
n = 10
Ns = exp10.(log10(Nmin):1/n:log10(Nmax))
dts = tf ./ Ns
vals = [integrate(u0, tf, vbp; dt=dt, adaptive=false).u[end] for dt in dts];


skip = 2
Ns_error = Ns[1+skip:end-skip]
errors = norm.(vals[1+2skip:end] .- vals[1:end-2skip])

Ns_order = Ns_error[1+skip:end-skip]
orders = -log.(errors[1+2skip:end] ./ errors[1:end-2skip]) ./ log.(Ns_error[1+2skip:end] ./ Ns_error[1:end-2skip])

logticks = 10. .^ (-16:16)
plot(Ns_error, errors, title="Error vs dt", label="", xaxis=:log, yaxis=:log, marker=:circle, xticks=logticks, yticks=logticks, minorgrid=true, xlabel="dt")
plot(Ns_order, orders, label="Order", marker=:circle, xaxis=:log, minorgrid=true, xlabel="dt")
plot!(Ns_order[1:1], [0], alpha=1, label="")

# Find first dt such that the order
i = findnext(orders .< 1, 1)

@info "Observation: the system is NOT stiff since there is no regime for large enough timesteps where the order of the method is not obtained; thus, an explicit method is enough."
@info f"First N such that the order is less than 1: \%(Ns_order[i]).
dt = tf / N = \%.1e(tf / Ns_order[i]).
We can safely choose 1e-3 (slightly less than the above value) as our dt, where numerical errors just begin to appear.
Note that it is important to normalize the quantities, else the data points are shifted and we would need to adjust our dt accordingly."
