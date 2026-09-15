# Spindle torus: set of possible kite positions
#
#   arm tip  A(α) = l * (cos α, sin α, 0)
#   kite     K(α, β) = A(α) + r * (cos β * (cos α, sin α, 0) + sin β * ẑ)
#
# The line of length r points anywhere in the radial-vertical plane. Revolving
# the sphere of radius r about the arm axis gives a torus; r > l => spindle.

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
    xlabel="x (m)", ylabel="y (m)", zlabel="z (m)", legend=:topright,
    title="Spindle torus, not to scale", label="",
    size=(800, 800), margin=1Plots.mm)

# set of possible arm-tip positions: circle of radius l in the xy-plane
θc = range(0, 2π, length=200)
plot!(fig, l * cos.(θc), l * sin.(θc), zero(θc); c=:blue, lw=3, label="arm-tip circle (radius l)")

# example arm and kite at α = 0, β = 30°
plot!(fig, [0, l], [0, 0], [0, 0]; c=:black, lw=3, label="arm (length l)")
plot!(fig, [l, kite_pos[1]], [0, 0], [0, kite_pos[3]]; c=:red, lw=2, ls=:dash, label="kite line")

# COMMENT/UNCOMMENT to identify easily those lines
plot!(fig, (l + r * cos(β0)) .* cos.(θc), (l + r * cos(β0)) .* sin.(θc), fill(r * sin(β0), length(θc)); c=:green, lw=2, label="kite latitude")
plot!(fig, l .+ r .* cos.(θc), zero(θc), r .* sin.(θc); c=:orange, lw=2, label="kite longitude")

scatter!(fig, [kite_pos[1]], [0], [kite_pos[3]]; c=:red, ms=7, label="kite")

savefig(fig, joinpath(@__DIR__, "spindle_torus.png"))
display(fig)
