using Pkg
Pkg.activate("CT_keep")

using OptimalControl
using ComponentArrays: ComponentArray as CA
using Plots

using KEEP.PointMassPara: build_vbpara
using KEEP.PointMass4: dynamics
using KEEP.TorqueFunction: torque_function
using KEEP.LimitCycle: compute_limit_cycle

const vbp0 = build_vbpara()
lc = compute_limit_cycle(vbp0; sense=+, save_everystep=true)

const syms = (:r, :I_eq, :torque_slope)

f(x, p) = begin
    vbp = CA(vbp0; (syms .=> p)...)
    x_dyn = vcat(x, 0)
    return dynamics(x_dyn, vbp)[1:4]
end

generated_power(dα, p) = begin
    vbp = CA(vbp0; (syms .=> p)...)
    return -dα(t) * torque_function(dα, vbp)
end

f(rand(4), [15, 50, 40])
torque(rand(), [15, 50, 40])

ocp = @def begin
    # tf ∈ R, variable
    # p ∈ R³, variable
    vars ∈ R⁴, variable
    tf = vars[1]
    p = vars[2:4]
    t ∈ [0, tf], time
    X = [α, τ, dα, dτ] ∈ R⁴, state

    X(0) - X(tf) == zeros(4)

    Ẋ(t) == f(X(t), p)

    ∫(generated_power(dα, p)) → max
end

init = @init ocp begin
    X(t) := lc(t)[1:4]
end
init = nothing

solve(ocp; initial_guess=init)
solve(ocp; initial_guess=init, backend=:manual)