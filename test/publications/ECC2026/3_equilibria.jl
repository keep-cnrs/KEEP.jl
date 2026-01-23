
import CairoMakie as cm
import Makie
using Colors
using LaTeXStrings
using SplitApplyCombine
using LinearAlgebra
using ForwardDiff: jacobian
using StaticArrays


import KEEP.PointMass4 as PM4
import KEEP.PointMassPara as PMP
using KEEP.SteadyState: ddq_partial, all_steady_states
using KEEP.Visualisation: TICKS_PI, TICKS_HALF_PI

include("plots_default.jl")

function main()
    PT = 4 / 3  # Makie point size ???

    p = PMP.build_vbpara()
    f = q -> ddq_partial(q .+ 1e-100, p)

    # For contour
    N = 20
    αs = range(-π, π, N + 1)
    τs = range(-π, π, 2N + 1)
    ddαs, ddτs = invert([f((α, τ)) for α in αs, τ in τs])

    @time ss = all_steady_states(p)
    function max_eigval(q)
        jacobian(u -> PM4.dynamics(u, p, 0), SA[q..., 0, 0, 0]) |> collect |> eigvals .|> real |> maximum
    end
    @time max_eigvals = max_eigval.(ss)

    # ticks = (-π:π/2:π, [j"-π", L"-\frac{π}{2}", L"0", L"\frac{π}{2}", L"π"])
    # ticks = (-π:π:π, [L"-π", L"0", L"π"])
    Makie.set_theme!(Makie.Theme(
        font="CMU Serif",
        fontsize=2.2 * TICKFONTSIZE,
    ))

    fig = cm.Figure(size=plot_size(4 / 3) .* PT, figure_padding=0)
    ax = Makie.Axis(fig[1, 1],
        xticks=([-π, 0, π], [L"-π", L"α = 0", L"π"]),
        yticks=([-π, 0, π], [L"-π", L"τ = 0", L"π"]),
        yticklabelrotation=π / 2,
        aspect=Makie.DataAspect(),
    )
    cm.hidedecorations!(ax, label=false, ticklabels=false, grid=false)
    cm.hidespines!(ax)
    cm.poly!([Makie.Rect(-π + i * π / 2, -π, π / 2, 2π) for i in 0:3], color=[:red, :green, :green, :red], alpha=0.1)
    col_α = PALETTE[4] # :steelblue2
    col_τ = PALETTE[3] # :chartreuse2

    cm.streamplot!(ax, cm.Point2 ∘ f, (-π, π), (-π, π); colormap=:binary, density=1, color=p -> norm(p), arrow_size=5, linewidth=0.4, gridsize=(20, 20), stepsize=0.1, maxsteps=10)
    # cm.streamplot!(ax, cm.Point2 ∘ f, (-π, π), (-π, π); arrow_size=5, linewidth=.4, density=1, gridsize=(32, 32), colorrange=(0,0),lowclip=colorant"lightgrey", highclip=colorant"lightgrey", alpha=.5)

    cm.contour!(ax, αs, τs, ddτs; levels=[0], color=col_α)
    cm.contour!(ax, αs, τs, ddαs; levels=[0], color=col_τ)
    cm.lines!(ax, 1, 1, color=col_α)
    cm.lines!(ax, 1, 1, color=col_τ)

    s = cm.scatter!(ax, invert(ss)...; markersize=10, color=max_eigvals, marker=:circle, strokewidth=0.5, colormap=Makie.Reverse(:imola))

    line_α = cm.LineElement(color=col_α)
    line_τ = cm.LineElement(color=col_τ)
    L"\max_i \mathrm{Re}(\lambda_i)"
    cb = cm.Colorbar(fig[1, 2], s, label=L"\max_i \ \mathrm{Re}(\lambda_i)", ticklabelfont="CMU Serif", height=cm.Relative(0.9))
    # cm.Legend(fig[1, 2], [line_α, line_τ, marker_ss], [L"\ddot{α}= 0", L"\ddot{τ}= 0", "equilibria"])

    cm.rowsize!(fig.layout, 1, cm.Aspect(1, 1.0))
    cm.resize_to_layout!(fig)


    # line width, border off, legend
    display(fig)

    # cm.save("test/publications/ECC2026/figs/equilibria.png", fig, px_per_unit=PT_PER_INCH)
    cm.save("test/publications/ECC2026/figs/equilibria.pdf", fig, px_per_unit=PT_PER_INCH * PT)

    #=
    We plot in grey the stream lines of the $\ddq(\vq)$, in purple and green the null-clines of the $\ddot{α}(\vq)$ and $\ddot{τ}(\vq)$ respectively, and mark the equilibria at their intersections with a color based on the maximum eigenvalue of the jacobian of $\vb x -> f(x) = \dot{\vb x}$ \alter{vérifier f(x)}. Background is green when the arm points downwind and red when the arm points upwind, which is to be avoided for ease of access to the device by a technician.
    =#

    # sym(q) = (-q[1], π + q[2])
    # f_sym(q) = (-1, 1) .* f(sym(q))
    ##
    τhats = [PM4.compute_τhat(q[2], p)[1:2] for q in ss]
    αhats = [PM4.compute_αhat(q[1])[1:2] for q in ss]
    dot.(τhats, αhats)  # Not aligned



    maximum.(eigvals.(jacobian.(f, ss)))
    u0 = SA[ss[4]..., 0, 0, 0]
    Plots.plot(PM4.integrate(u0, 15, p, save_everystep=true))

    inch = 96
    pt = 4 / 3

    _PT_PER_INCH = inch * pt
end

main()