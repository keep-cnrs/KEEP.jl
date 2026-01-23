using Plots

τ = range(0, 360, length=2001)
θ = τ -> sind(2 * τ)
φ = τ -> sind(τ)

fig = plot(φ.(τ), θ.(τ), label="")
plot!(xlabel="\$ϕ(τ)\$", ylabel="\$θ(τ)\$",
    xticks=([-1, 0, 1], ("\$ϕ_0 - Δϕ\$", "\$ϕ_0\$", "\$ϕ_0 + Δϕ\$")),
    yticks=([-1, 0, 1], ("\$θ_0 - Δθ\$", "\$θ_0\$", "\$θ_0 + Δθ\$")),
    title="Base of the cone")
plot!(size=(600, 300))

τ_φ = 90.0
plot!(φ.([0, τ_φ]), θ.([τ_φ, τ_φ]), arrow=true, color=:black, label="")
plot!(φ.([τ_φ, 0]), θ.([τ_φ, τ_φ]), arrow=true, color=:black, label="")
annotate!(sum(φ.([0, τ_φ])) / 2, sum(θ.([τ_φ, τ_φ])) / 2, "\$Δϕ\$", :bottom)

τ_θ = 225.0
plot!(φ.([τ_θ, τ_θ]), θ.([0, τ_θ]), arrow=true, color=:black, label="")
plot!(φ.([τ_θ, τ_θ]), θ.([τ_θ, 0]), arrow=true, color=:black, label="")
annotate!(sum(φ.([τ_θ, τ_θ])) / 2, sum(θ.([τ_θ, 0])) / 2, "\$Δθ\$", :left)

display(fig)

## For presentation
fig = plot(φ.(τ), θ.(τ), lw=3, label="")
plot!(xlabel="\$ϕ(τ)\$", ylabel="\$θ(τ)\$",
    xticks=([-1, 0, 1], ("\$ϕ_0 - Δϕ\$", "\$ϕ_0\$", "\$ϕ_0 + Δϕ\$")),
    yticks=([-1, 0, 1], ("\$θ_0 - Δθ\$", "\$θ_0\$", "\$θ_0 + Δθ\$")),
    title="Base of the cone")
plot!(size=(400, 200), right_margin=3Plots.mm)

τ_φ = 90.0
plot!(φ.([0, τ_φ]), θ.([τ_φ, τ_φ]), arrow=true, color=:black, label="")
plot!(φ.([τ_φ, 0]), θ.([τ_φ, τ_φ]), arrow=true, color=:black, label="")
annotate!(sum(φ.([0, τ_φ])) / 2, sum(θ.([τ_φ, τ_φ])) / 2, "\$Δϕ\$", :bottom)

τ_θ = 225.0
plot!(φ.([τ_θ, τ_θ]), θ.([0, τ_θ]), arrow=true, color=:black, label="")
plot!(φ.([τ_θ, τ_θ]), θ.([τ_θ, 0]), arrow=true, color=:black, label="")
annotate!(sum(φ.([τ_θ, τ_θ])) / 2, sum(θ.([τ_θ, 0])) / 2, "\$Δθ\$", :left)

display(fig)
# savefig(fig, "docs/media/eight_plot.pdf")
