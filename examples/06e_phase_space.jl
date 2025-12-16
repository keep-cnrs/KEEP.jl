#=
mêmes symétries que (α, τ) -> (cos(α), sin(α), sin(τ), sin(2τ))
equality :
    ddα(-α, π+τ) = -ddα(α, τ)
    ddτ(-α, π+τ) = ddτ(α, τ)

approximate but not equal :
    ddα(-α, -τ) ≈ -ddα(α, τ)
    ddτ(α + π) ≈ ddτ(α, τ)


Domain Ω = (α, τ) ∈ [0, 2π] × [0, 2π], symmetry around α = 0 and α = π and τ = 0 and τ = π
Ω = Ω11 ∪ Ω12 ∪ Ω21 ∪ Ω22 = ([0, π] × [0, π]) ∪ ([0, π] × [π, 2π]) ∪ ([π, 2π] × [0, π]) ∪ ([π, 2π] × [π, 2π])
Il suffit de travailler dans Ω1 = Ω11 ∪ Ω12 par exemple, grâce aux symmetries.

sampler Ω1 = [0, π] × [0, 2π] (Halton)

=#

import CairoMakie as cm
import Makie
using LaTeXStrings
using SplitApplyCombine
using LinearAlgebra: norm

using KEEP.SteadyState: ddq_partial
using KEEP.PointMassPara
using KEEP.Visualisation: TICKS_PI, TICKS_HALF_PI

vbp = build_vbpara()
f = q -> ddq_partial(q .+ 1e-100, vbp)

# For contour
N = 100
αs = range(-π, π, N + 1)
τs = range(-π, π, 2N + 1)
ddαs, ddτs = invert([f((α, τ)) for α in αs, τ in τs])

# ticks = (-π:π/2:π, [j"-π", L"-\frac{π}{2}", L"0", L"\frac{π}{2}", L"π"])
# ticks = (-π:π:π, [L"-π", L"0", L"π"])

fig = cm.Figure(fontsize=32)
ax = Makie.Axis(fig[1,1],
    title = "Dynamics with zero velocity",
    xlabel=L"α", xticks=TICKS_PI,
    ylabel=L"τ", yticks=TICKS_PI
)
cm.streamplot!(ax, cm.Point2 ∘ f, (-π, π), (-π, π); colormap=:binary, density=1, color=p -> norm(p)^2)
cm.contour!(ax, αs, τs, ddτs; levels=[0], color=:cyan, linewidth=5, label="ddα = 0")
cm.contour!(ax, αs, τs, ddαs; levels=[0], color=:orange, linewidth=5, label="ddτ = 0")
cm.lines!(ax, 1, 1, color=:orange, linewidth=5, label="ddot{α}= 0")
cm.lines!(ax, 1, 1, color=:cyan, linewidth=5, label="ddot{τ}= 0")
line_α = cm.LineElement(color=:orange, linewidth=5)
line_τ = cm.LineElement(color=:cyan, linewidth=5)
cm.Legend(fig[1, 2], [line_α, line_τ], [L"\ddot{α}= 0", L"\ddot{τ}= 0"])
display(fig)
# save("tmp.png", fig)

# sym(q) = (-q[1], π + q[2])
# f_sym(q) = (-1, 1) .* f(sym(q))
