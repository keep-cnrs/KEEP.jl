
using StaticArrays
using Plots
using SplitApplyCombine
using LaTeXStrings
using PyFormattedStrings
using Roots
using LinearAlgebra: norm, ⋅, ×

import KEEP.PointMass4 as PM4
import KEEP.PointMassPara as PMP
import KEEP.SteadyState as SS
import KEEP.TorqueFunction as TF

include("state_funcs.jl")
include("plots_default.jl")

function main()
    p = PMP.build_vbpara()

    α0, τ0, dα0, dτ0 = 0, 1, 0, 0
    tf = 100
    sol_plus = PM4.integrate(SA[α0, 1, dα0, dτ0, 0], tf, p, save_everystep=true)
    sol_minus = PM4.integrate(SA[α0, -1, dα0, dτ0, 0], tf, p, save_everystep=true)
    # TOP FRONT LEFT view
    # using code in @Visualisation.jl , write a function here that will plot the trajectory in XY, XZ and YZ planes, aswell as from the current angle

    t = range(tf - 5, tf, step=0.01)
    plot(aspect_ratio=:equal, xlabel="\$y\$ (m)", ylabel="\$z\$ (m)", xtick=-20:5:20, ytick=5:2.5:10, size=plot_size(1.8))
    plot!(invert(r_vec.(t, Ref(sol_plus)))[2:3]..., label="\$\\dot{\\tau} > 0\$")
    plot!(invert(r_vec.(t, Ref(sol_minus)))[2:3]..., label="\$\\dot{\\tau} < 0\$")
    lens!([-11, -10], [8, 9], inset=(1, bbox(0.25, 0, 0.3, 0.3)))
    l = plot!().inset_subplots[1]
    xticks!(l, -11:0.5:-10)
    yticks!(l, 8:0.5:9)
    plot!([0], [13])
    display(plot!())

    # DEPRECATED, use the one found in limit_cycles.jl
    # savefig("test/publications/ECC2026/figs/two_trajectories.pdf")

    ##
    plot(title="Two different limit cycles", aspect_ratio=:equal, axis=false, xlabel=L"\alpha", ylabel=L"\tau")
    plot!(α.(t, Ref(sol_plus)), τ.(t, Ref(sol_plus)) .% 2π, label="+")
    plot!(α.(t, Ref(sol_minus)), τ.(t, Ref(sol_minus)) .% 2π, label="-")
    display(plot!())


    ## This is very cool
    xyz_torus = (t, sol, sym=false) -> begin
        α_ = α(t, sol)
        τ_ = τ(t, sol)
        if sym
            α_ = -α_
            τ_ = -τ_
        end
        x = (2 + cos(τ_)) * cos(α_)
        y = (2 + cos(τ_)) * sin(α_)
        z = 1.5 * sin(τ_)
        return SA[x, y, z]
    end

    lw = 3
    plot_sym = false
    α_ = τ_ = range(-π, π, 50)
    X = [(2 + cos(v)) * cos(u) for u in α_, v in τ_]
    Y = [(2 + cos(v)) * sin(u) for u in α_, v in τ_]
    Z = [1.5 * sin(v) for u in α_, v in τ_]
    surface(X, Y, Z, lims=(-3, 3), size=(600, 600), cbar=:none, legend=false, alpha=0.3, color=:diff, title="A trajectory and its symmetry (opposite phase)", axis=false)
    plot!(invert(xyz_torus.(t, Ref(sol_plus)))..., lw=lw, label="+", c=1)
    plot!(invert(xyz_torus.(t, Ref(sol_minus)))..., lw=lw, label="-", c=2)
    if plot_sym
        plot!(invert(xyz_torus.(t, Ref(sol_plus), true))..., lw=lw, label="sym(+)", c=1, alpha=0.7)
        plot!(invert(xyz_torus.(t, Ref(sol_minus), true))..., lw=lw, label="sym(-)", c=2, alpha=0.7)
    end
    main = deepcopy(plot!(title=""))
    f1 = deepcopy(plot!(camera=(0, 90), title=""))
    f2 = deepcopy(plot!(camera=(90, 0), title=""))
    f3 = deepcopy(plot!(camera=(90, 90), title=""))
    plot(f1, f2, f3, main; layout=4)
    display(plot!())
    gr()
    # Check: the two limit cycles are invariant under the symmetry of the system (opposing the phase)
    # The two limit cycles are clearly distinct in phase space


    ## Finding a point with zero τ velocity at any t > 0 (works for 5 seconds)
    u0(τ0) = SA[π/4, τ0, 0, 0, 0]
    f = (t) -> ((τ0) -> last(PM4.integrate(u0(τ0), t, p).u)[4])

    τ_guess = 1e-100
    tf = 10
    for tf in 1:tf
        tol = 1e-10
        τ_guess = find_zero(f(tf), τ_guess, xatol=tol, xrtol=tol, atol=tol, rtol=tol)
        @show τ_guess, f(tf)(τ_guess)
    end
    τ_guess

    s = PM4.integrate(u0(τ_guess), tf, p, tol=1e-9, save_everystep=true)
    plot(s, idxs=1:4)
end

main()