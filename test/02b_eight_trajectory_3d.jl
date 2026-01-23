using StaticArrays
using Plots

θ0 = deg2rad(60)
φ0 = 0
Δθ = deg2rad(5)
Δφ = deg2rad(45)

θ(τ) = θ0 + Δθ * sin(2 * τ)
φ(τ) = φ0 + Δφ * sin(τ)

function rθφ2xyz(r, θ, φ)
    x = r * sin(θ) * cos(φ)
    y = r * sin(θ) * sin(φ)
    z = r * cos(θ)
    return SA[x, y, z]
end

τ = range(0, 2π, length=2001)
eight_xyz = rθφ2xyz.(1, θ.(τ), φ.(τ))
segments = stack(reduce.(hcat, (zero(eight_xyz), eight_xyz)); dims=2)

camera_0 = [60, 50]
camera_Δ = [10, 5]

plot(Tuple.(eight_xyz), lw=3)
plot!(eachslice(segments; dims=1)..., c=:black, alpha=20 / length(τ), label="")
plot!(xlabel="x", ylabel="y", zlabel="z", xlims=(0, 1), ylims=(-1, 1), zlims=(0, 1), camera=camera_0, size=(600, 600))

if false
    anim = @animate for α in range(0, 2π, length=60)
        camera = camera_0 + camera_Δ .* [sin(α), cos(α)]
        plot!(camera=camera)
    end
    # gif(anim, "docs/media/eight_anim.gif"; fps=30)
end
