using KEEP: PointMass10 as PM10
using KEEP: PointMass4 as PM4
using KEEP: Visualization as VIZ

using OrdinaryDiffEqVerner: Vern9

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
display(plot!())
VIZ.plot_trajectory_10D(sol_cb)
display(plot!())

plot(title="Configuration residuals without vs. with callback")
plot!(sol.t, [norm(PM10.manifold_residuals!(similar(u0, 6), u, p)) for u in sol.u], label="No callback")
plot!(sol_cb.t, [norm(PM10.manifold_residuals!(similar(u0, 6), u, p)) for u in sol_cb.u], label="Callback")
hline!([tol], label="Tolerance", c=:black)
display(plot!(xlabel="Time (s)", yscale=:log, yticks=exp10.(-15:3:3), legend=:bottomright))



vbp = PM4.build_vbpara()
u0 = PM4.init_u(τ=τ0)
sol = PM4.integrate(u0, tf, vbp; save_everystep=true, tol)

# Run a low tolerance (Alg = Vern9, abstol = reltol 1e-10) integration and compare error between that and 10D, 10D+callback and 4D in configuration space (x, y, z)
throw("implement comment above and add line to conclusion")

btime4 = @belapsed PM4.integrate(u0, tf, vbp; save_everystep=true, tol) seconds = 0.1

btime = @belapsed PM10.integrate(u0, tf, p; save_everystep=true, tol) seconds=0.1
btime_cb = @belapsed PM10.integrate(u0, tf, p; save_everystep=true, tol, callback) seconds=0.1

throw("add PM4 timing and update conclusion (it can be called 'proposed model')")

ratio_cb = btime_cb / btime
@info f"Timing
reference (10D, no callback) = {btime:.2g} ms
reference + callback = {btime_cb:.2g} ms

reference to reference+callback = {ratio:.2f}x slowdown
reference to better model = {something:.2f}x speedup"

@info "Conclusion
The overdetermined system, which requires carefull initialization (cf. algebraic constraints) experiences drifting when numerically integrated, which is not suitable for finer analysis.
By adding a manifold projection callback, we effectively nullify the drift -- at the cost of performance (about 2.5 times slower), and we would still be working with an inherently inconsistent state (residuals are not exactly 0)
The proposed model"
