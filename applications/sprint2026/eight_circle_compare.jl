# Compare the two "eight + circle" layouts on one limit cycle, so the layout can
# be chosen before rendering the videos:
#   left  :src       — true scale (faithful to KEEP.Visualization)
#   right :rescaled  — each shape normalized to its own extent (legible)
#
# Run: julia --project=applications/sprint2026 applications/sprint2026/eight_circle_compare.jl
using Plots
using Serialization
import BifurcationKit                                  # resolve RecursiveArrayTools on deserialize
import KEEP.PointMassPara: build_para, build_vbpara

const OPT = (r=47.44987027243979, I_eq=3541.2653832051565, torque_slope=3862.561181135744)
const HERE = @__DIR__

include(joinpath(HERE, "eight_circle.jl"))
using .EightCircle

vbp = build_vbpara(build_para(r=OPT.r, I_eq=OPT.I_eq, torque_slope=OPT.torque_slope))
cyc = deserialize(joinpath(HERE, "scratch", "three_cycles_shooting.jls"))
U = Array(cyc.short.u)
Qs = [(U[1, i], U[2, i]) for i in 1:size(U, 2)]
nh = round(Int, 0.45 * length(Qs))

plt = plot(layout=(1, 2), size=(1500, 780), plot_title="eight-circle layout — short cycle @ v_ref=2.48")
for (k, mode) in enumerate((:src, :rescaled))
    g = eight_circle_geometry(vbp; mode=mode)
    sp = plt[k]
    plot!(sp; aspect_ratio=:equal, axis=false, grid=false,
        title=(mode === :src ? "src proportions (true scale, circle tiny)" :
                              "rescaled (per-shape, both legible)"))
    draw_background!(sp, g)
    add_trail!(sp, Qs[1:3:end], g; color=:royalblue, alpha=0.10, lw=1, label="")   # faint full cycle
    add_trail!(sp, Qs[1:nh], g; color=:royalblue, alpha=0.55, lw=1.5, label="")    # traced so far
    add_state!(sp, Qs[nh], g; color=:crimson, lw=2.5, label="")                    # current spoke
end
savefig(plt, joinpath(HERE, "fig_eight_circle_compare.png"))
println("saved fig_eight_circle_compare.png")
