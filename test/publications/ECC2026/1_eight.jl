## θ φ plot with flipped y axis + arrows for orientation of tau

using Plots
using LaTeXStrings
using SplitApplyCombine

include("plots_default.jl")

function main()
    τ = range(0, 360, length=100)
    θ = τ -> sind(2 * τ)
    φ = τ -> sind(τ)

    plot(xmirror=true, size=plot_size(2))
    # plot!(xlabel="\$ϕ(τ)\$", ylabel="\$θ(τ)\$")
    plot!(xticks=([-1, 0, 1], ("\$Δφ\$", "\$φ = 0\$", "\$Δφ\$")),
        yticks=([-1, 0, 1], ("\$θ_0 + Δθ\$", "\$θ = θ_0\$", "\$θ_0 - Δθ\$"))
    )

    # plot!(
    #     xticks=([-1, 0, 1], ("", "\$0\$", "")),
    #     yticks=([-1, 0, 1], ("", "\$θ_0\$", ""))
    # )

    plot!(φ.(τ), θ.(τ), c=:grey)

    delta_col = :black
    τ_φ = 90.0
    plot!(φ.([0, τ_φ]), θ.([τ_φ, τ_φ]), arrow=false, color=delta_col, label="")
    # plot!(φ.([τ_φ, 0]), θ.([τ_φ, τ_φ]), arrow=false, color=delta_col, label="")
    annotate!(sum(φ.([0, τ_φ])) / 2, sum(θ.([τ_φ, τ_φ])) / 2 - 0.1, ("\$Δφ\$", ANNOTATIONFONTSIZE, delta_col, :top))

    τ_θ = 225.0
    # plot!(φ.([τ_θ, τ_θ]), θ.([0, τ_θ]), arrow=false, color=delta_col, label="")
    plot!(φ.([τ_θ, τ_θ]), θ.([τ_θ, 0]), arrow=false, color=delta_col, label="")
    annotate!(sum(φ.([τ_θ, τ_θ])) / 2 + 0.05, sum(θ.([τ_θ, 0])) / 2, ("\$Δθ\$", ANNOTATIONFONTSIZE, delta_col, :left))

    dx, dy = 1 / 4, 1 / 2
    up_col, down_col = 1, 2
    quiver!([0], [0], quiver=([dx], [-dy]), c=up_col)
    quiver!([0], [0], quiver=([-dx], [dy]), c=down_col)
    # annotate!(-dx, -dy-0.1, ("\$τ\$↗", ANNOTATIONFONTSIZE, :center, up_col))
    # annotate!(dx, dy+0.1, ("\$τ\$↘", ANNOTATIONFONTSIZE, :center, down_col))
    annotate!(dx - 0.05, -dy - 0.25, ("\$\\dot{\\tau} > 0\$", ANNOTATIONFONTSIZE, :center, palette(:auto)[up_col]))
    annotate!(-dx + 0.05, dy + 0.25, ("\$\\dot{\\tau} < 0\$", ANNOTATIONFONTSIZE, :center, palette(:auto)[down_col]))

    display(plot!())

    savefig(plot!(), "test/publications/ECC2026/figs/eight_plot.pdf")
end

main()