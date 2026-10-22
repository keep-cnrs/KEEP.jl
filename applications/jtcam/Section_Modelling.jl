# Spindle torus: set of possible kite positions
#
#   arm tip  A(α) = l * (cos α, sin α, 0)
#   kite     K(α, β) = A(α) + r * (cos β * (cos α, sin α, 0) + sin β * ẑ)
#
# The line of length r points anywhere in the radial-vertical plane. Revolving
# the sphere of radius r about the arm axis gives a torus; r > l => spindle.
using Pkg
Pkg.activate(@__DIR__)

using Plots
using LinearAlgebra: norm

l, r = 2, 3

# sanity check: a sampled point is exactly r from its arm tip
β0 = deg2rad(30)  # line elevation
kite_pos = [l + r * cos(β0), 0, r * sin(β0)]
@assert isapprox(norm(kite_pos - [l, 0, 0]), r; atol=1e-8)

# mesh as flat line segments (GR's wireframe chokes on matrix grids); a mesh
# rather than a surface so the self-intersecting interior stays visible
lg = range(0, 2π, length=120)  # points along each grid line
lk = range(0, 2π, length=40)   # grid lines
xs = Float64[];
ys = Float64[];
zs = Float64[]
for βj in lk
    for αi in lg
        push!(xs, (l + r * cos(βj)) * cos(αi))
        push!(ys, (l + r * cos(βj)) * sin(αi))
        push!(zs, r * sin(βj))
    end
    push!(xs, NaN)
    push!(ys, NaN)
    push!(zs, NaN)
end
for αj in lk
    for βi in lg
        push!(xs, (l + r * cos(βi)) * cos(αj))
        push!(ys, (l + r * cos(βi)) * sin(αj))
        push!(zs, r * sin(βi))
    end
    push!(xs, NaN)
    push!(ys, NaN)
    push!(zs, NaN)
end

fig = plot(xs, ys, zs; c=:black, alpha=0.4, lw=0.5, aspect_ratio=:equal,
    xlabel="x (m)", ylabel="y (m)", zlabel="z (m)",
    title="Spindle torus, not to scale", label="",
    size=SIZE_HALF .* (1.2, 1), margin=1Plots.mm, legend=:outerright)

# set of possible arm-tip positions: circle of radius l in the xy-plane
θc = range(0, 2π, length=200)
arm_tip_pars = (; c=:blue, lw=3, label="Arm-tip locus")
plot!(fig, l * cos.(θc), l * sin.(θc), zero(θc); arm_tip_pars..., label="")

# example arm and kite at α = 0, β = 30°
arm_and_lines_pars = (; c=:black, lw=3, label="Arm and lines")
plot!(fig, [0, l], [0, 0], [0, 0]; arm_and_lines_pars..., label="")
plot!(fig, [l, kite_pos[1]], [0, 0], [0, kite_pos[3]]; arm_and_lines_pars..., label="")

# COMMENT/UNCOMMENT to identify easily those lines
longitude_pars = (; c=:orange, lw=2, label="Kite longitude")
latitude_pars = (; c=:green, lw=2, label="Kite latitude")
plot!(fig, l .+ r .* cos.(θc), zero(θc), r .* sin.(θc); longitude_pars..., label="")
plot!(fig, (l + r * cos(β0)) .* cos.(θc), (l + r * cos(β0)) .* sin.(θc), fill(r * sin(β0), length(θc)); latitude_pars..., label="")

kite_pars = (; c=:red, ms=7, label="Kite")
scatter!(fig, [kite_pos[1]], [0], [kite_pos[3]]; kite_pars..., label="")

scatter!([]; kite_pars...)
plot!([]; arm_and_lines_pars...)
plot!([]; arm_tip_pars...)
plot!([]; longitude_pars...)
plot!([]; latitude_pars...)

savefig(fig, joinpath(@__DIR__, "figs", "spindle_torus.pdf"))
display(fig)

## Parametrisation in tau
# TODO [use the plot from test/publications/ECC2026]

τ_φ = 90.0
τ_θ = 225.0
dx, dy = 1 / 4, 1 / 2
up_col, down_col = 1, 2

τ = range(0, 360, length=501)
θ = τ -> sind(2 * τ)
φ = τ -> sind(τ)

plot(xmirror=true)
plot!(xticks=([-1, 0, 1], (raw"$-Δφ$", raw"$φ = 0$", raw"$Δφ$")),
    yticks=([-1, 0, 1], (raw"$θ_0 + Δθ$", raw"$θ = θ_0$", raw"$θ_0 - Δθ$"))
)

plot!(φ.(τ), θ.(τ), c=:grey)

plot!(φ.([0, τ_φ]), θ.([τ_φ, τ_φ]), color=:black)
annotate!(sum(φ.([0, τ_φ])) / 2, sum(θ.([τ_φ, τ_φ])) / 2 - 0.1, (raw"$Δφ$", ANNOTATIONFONTSIZE, :black, :top))

plot!(φ.([τ_θ, τ_θ]), θ.([τ_θ, 0]), color=:black)
annotate!(sum(φ.([τ_θ, τ_θ])) / 2 + 0.05, sum(θ.([τ_θ, 0])) / 2, (raw"$Δθ$", ANNOTATIONFONTSIZE, :black, :left))
quiver!([0], [0], quiver=([dx], [-dy]), c=up_col)
quiver!([0], [0], quiver=([-dx], [dy]), c=down_col)
annotate!(dx - 0.05, -dy - 0.25, (raw"$\dot{\tau} > 0$", ANNOTATIONFONTSIZE, :center, PALETTE[up_col]))
annotate!(-dx + 0.05, dy + 0.25, (raw"$\dot{\tau} < 0$", ANNOTATIONFONTSIZE, :center, PALETTE[down_col]))
plot!(size=SIZE_HALF .* (1, 0.5), topmargin=1Plots.mm)
display(current())
savefig(current(), joinpath(OUTPATH, "eight_plot.pdf"))