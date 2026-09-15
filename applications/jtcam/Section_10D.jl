using KEEP: PointMass10 as PM10
using KEEP: Visualization as VIZ

using Plots
using PyFormattedStrings
using Logging

using LinearAlgebra: norm

p = PM10.build_para()
τ0, dτ0 = 1e-10, 10
u0 = PM10.init_u(τ0, dτ0, p)
tf, tol = 80, 1e-3
sol = PM10.integrate(u0, tf, p; save_everystep=true, tol)
callback = PM10.build_manifold_projection(u0; save=true, tol)
sol_cb = PM10.integrate(u0, tf, p; save_everystep=false, tol, callback)

VIZ.plot_trajectory_10D(sol)
VIZ.plot_trajectory_10D(sol_cb)


plot(title="Residuals without vs. with callback")
plot!(sol.t, [norm(PM10.manifold_residuals!(similar(u0, 6), u, p)) for u in sol.u], label="Without")
plot!(sol_cb.t, [norm(PM10.manifold_residuals!(similar(u0, 6), u, p)) for u in sol_cb.u], label="With")
hline!([tol], label="Tolerance", c=:black)
display(plot!(yscale=:log, yticks=exp10.(-15:3:3), legend=:bottomright))

btime = @belapsed PM10.integrate(u0, tf, p; save_everystep=true, tol) seconds=0.1
btime_cb = @belapsed PM10.integrate(u0, tf, p; save_everystep=true, tol, callback) seconds=0.1

ratio = btime_cb / btime
@info f"Callback-induced slowdown
without = {btime:.2g} ms
with = {btime_cb:.2g} ms
ratio = {ratio:.2f}x"

@info "Conclusion
The overdetermined system, which requires carefull initialization (cf. algebraic constraints) experiences drifting when numerically integrated, which is not suitable for finer analysis.
By adding a manifold projection callback, we effectively nullify the drift -- at the cost of performance (about 2.5 times slower), and we would still be working with an inherently inconsistent state (residuals are not exactly 0)"