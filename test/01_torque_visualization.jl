using Plots

using KEEP.PointMassPara
using KEEP.TorqueFunction

## Qualitative representation of the torque function
# adjust xmax to see beyond Ωmax
xmax = 1.2
Xs = [-xmax, -1, -1, -1 / 3, 1 / 3, 1, 1, xmax]
Ys = [0, 0, -1, 0, 0, 1, 0, 0]

xtick_positions = [-1, -1 / 3, 0, 1 / 3, 1]
xtick_labels = [raw"$-Ω_\mathrm{max}$", raw"$-Ω_\mathrm{min}$", raw"$0$", raw"$Ω_\mathrm{min}$", raw"$Ω_\mathrm{max}$"]

ytick_positions = [-1, 0, 1]
ytick_labels = [raw"$-C_\mathrm{max}$", raw"$0$", raw"$C_\mathrm{max}$"]

fig = plot(Xs, Ys, label="")
plot!(xlabel=raw"$ω$", ylabel="torque",
    xticks=(xtick_positions, xtick_labels),
    yticks=(ytick_positions, ytick_labels),
    title="Torque from arm to generator"
)
display(fig)

## Quantitative test of the torque function
p = build_para()
scaled_Xs = p.Ωmax .* LinRange(-2, 2, 1000)
fig = plot(scaled_Xs, ω -> torque_function(ω, p),
    c=:red, lw=8,
    label="Torque function")
plot!(p.Ωmax .* Xs, p.Cmax .* Ys,
    ls=:dash, c=:black, lw=2,
    label="The function should follow this")
plot!(xlabel="ω (rad/s)", ylabel="torque (Nm)", legend=:bottomright,
    title="Torque from arm to generator")
display(fig)


## All torque functions
plot(xlabel="ω (rad/s)", ylabel="torque (Nm)", legend=:bottomright,
title="Torque from arm to generator")
plot!(scaled_Xs, ω -> torque_function(ω, p))

## For presentation
fig = plot(Xs, Ys, lw=3, label="")
plot!(xlabel=raw"$\dot α$", ylabel="torque",
    xticks=(xtick_positions, xtick_labels),
    yticks=(ytick_positions, ytick_labels),
    title="Torque from arm to generator",
    size=(400, 300)
)
# savefig(fig, "docs/media/torque_function.pdf")
