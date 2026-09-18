ENV["GKSwstype"] = "100"        # headless GR — must precede `using Plots`
using Plots
using Serialization
import KEEP.PointMassPara: build_para, build_vbpara

const HERE = @__DIR__
const OPT = (r=47.44987027243979, I_eq=3541.2653832051565, torque_slope=3862.561181135744)
const EC_MODE = Symbol(get(ENV, "EC_MODE", "src"))            # :src (true scale) | :rescaled
const ARM_SCALING = parse(Float64, get(ENV, "ARM_SCALING", "2"))
const FPS = 30
const DURATION = 10.0
const NFRAMES = round(Int, FPS * DURATION)                    # one full period in 10 s
const SIZE = (2400, 800)
const SPOKE_LW = 5.0

include(joinpath(HERE, "eight_circle.jl"))
using .EightCircle

vbp = build_vbpara(build_para(r=OPT.r, I_eq=OPT.I_eq, torque_slope=OPT.torque_slope))
geom = eight_circle_geometry(vbp; mode=EC_MODE, arm_scaling=ARM_SCALING)
cyc = deserialize(joinpath(HERE, "scratch", "three_cycles_frames.jls"))
@assert cyc.SNAPS == NFRAMES "frames file has $(cyc.SNAPS) phases, need $NFRAMES"

cycles = [
    (name="short",     c=cyc.short,    col=:royalblue),
    (name="long",      c=cyc.long,     col=:crimson),
    (name="long-long", c=cyc.longlong, col=:darkorange),
]

function build_frame(i)
    plt = plot(layout=(1, 3), size=SIZE,
        plot_title="Three prograde limit cycles — each time-scaled to 10 s (v_ref = 2.48)")
    for k in 1:3
        cy = cycles[k].c
        sp = plt[k]
        plot!(sp; aspect_ratio=:equal, axis=false, grid=false,
            title="$(cycles[k].name)   tf = $(round(cy.tf, digits=1)) s")
        draw_background!(sp, geom)
        add_trail!(sp, cy.trail, geom; color=:gray70, alpha=0.25, lw=1, label="")   # faint full cycle
        add_state!(sp, (cy.frames[1, i], cy.frames[2, i]), geom;
            color=cycles[k].col, lw=SPOKE_LW, label="")                             # moving spoke
    end
    return plt
end

if haskey(ENV, "EC_TEST")
    i = round(Int, parse(Float64, get(ENV, "EC_TEST_S", "0.45")) * (NFRAMES - 1)) + 1
    savefig(build_frame(i), joinpath(HERE, "scratch", "anim_testframe.png"))
    println("saved scratch/anim_testframe.png (i=$i)")
else
    t0 = time()
    anim = @animate for i in 1:NFRAMES
        build_frame(i)
    end
    mp4(anim, joinpath(HERE, "anim_three_cycles.mp4"); fps=FPS)
    rm(anim.dir; recursive=true, force=true)
    println("saved anim_three_cycles.mp4  ($(round(Int, time() - t0)) s)")
end
