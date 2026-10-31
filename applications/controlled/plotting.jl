
# =========================================================
# 3D Kite visualization: arm + tether + correctly oriented wing panels
# =========================================================

"""
    plot_kite!(plt, q; u=(0.0, 0.0), chord=1.0, halfspan=1.25, front=false,
               labels=true, shade=false, project=false, lims=nothing)

Draw one kite state `q = (α, θ, φ, β)` on `plt`: pivot→arm, arm→CG tether and
the left (+y) / right (−y) wing panels with the same sweep `delta` and body
orientation as the aerodynamic model. `u = (uₗ, uᵣ)` are the aileron
deflections: each panel is rotated about its span axis with the aerodynamic
sign convention (`i = atan(...) + uᵢ`, i.e. positive `uᵢ` increases that
panel's angle of attack). `front=true` adds an arrow along the nose axis
`Ihat`. When `lims` (the tuple returned in `kite_paths`) is given,
`shade`/`project` add a small CG marker shadow on the floor / back & side
walls.
"""
function plot_kite!(
    plt,
    q;
    u=(0.0, 0.0),
    chord=1.0,
    halfspan=1.25,
    front=false,
    labels=true,
    shade=false,
    project=false,
    lims=nothing,
)
    Ihat, Jhat, Khat = body_axes(q)
    cg = pos_CG(q)
    w2w(p) = cg + p[1] * Ihat + p[2] * Jhat + p[3] * Khat   # body → world

    # arm (pivot → arm tip) and tether (arm tip → CG)
    arm = SVector(pars.larm * cos(q[1]), pars.larm * sin(q[1]), 0.0)
    plot!(
        plt,
        [0.0, arm[1]],
        [0.0, arm[2]],
        [0.0, arm[3]];
        color=:black,
        lw=4,
        label=(labels ? "arm" : false),
    )
    plot!(
        plt,
        [arm[1], cg[1]],
        [arm[2], cg[2]],
        [arm[3], cg[3]];
        color=:orange,
        lw=2,
        label=(labels ? "tether" : false),
    )

    # panels: span axis w and chord axis l, swept by delta; the aileron uᵢ rotates
    # panel i about its span axis (same convention as aerodynamics(), where
    # i = atan(...) + uᵢ: positive uᵢ increases the panel's angle of attack)
    sd, cd = sincos(pars.aero.delta)
    panels = (
        (
            pars.aero.panels[:, 1],
            SVector(-sd, cd, 0.0),
            SVector(cd, sd, 0.0),
            :royalblue,
            "left wing",
        ),
        (
            pars.aero.panels[:, 2],
            SVector(-sd, -cd, 0.0),
            SVector(cd, -sd, 0.0),
            :tomato,
            "right wing",
        ),
    )
    for (i, (R, w, l, col, lbl)) in enumerate(panels)
        l = l * cos(u[i]) + cross3(w, l) * sin(u[i])
        corners = [
            w2w(R + a * halfspan * w + b * (chord / 2) * l) for
            a in (-1.0, 1.0), b in (-1.0, 1.0)
        ]
        quad = corners[[1, 3, 4, 2, 1]]     # closed rectangle loop (column-major)
        plot!(
            plt,
            [p[1] for p in quad],
            [p[2] for p in quad],
            [p[3] for p in quad];
            color=col,
            lw=1.5,
            label=(labels ? lbl : false),
        )
        # ponytail: GR can't png-export filled 3D quads (firstindex(::Surface)) —
        # outline + span/chord cross reads as a plane; swap to surface! if backend allows
        s1, s2 = w2w(R - halfspan * w), w2w(R + halfspan * w)
        c1, c2 = w2w(R - (chord / 2) * l), w2w(R + (chord / 2) * l)
        plot!(
            plt,
            [s1[1], s2[1]],
            [s1[2], s2[2]],
            [s1[3], s2[3]];
            color=col,
            lw=1,
            alpha=0.6,
            label=false,
        )
        plot!(
            plt,
            [c1[1], c2[1]],
            [c1[2], c2[2]],
            [c1[3], c2[3]];
            color=col,
            lw=1,
            alpha=0.6,
            label=false,
        )
    end

    if front
        tip = cg + 2.0 * Ihat
        plot!(
            plt,
            [cg[1], tip[1]],
            [cg[2], tip[2]],
            [cg[3], tip[3]];
            arrow=arrow(:closed, :head, 0.3, 0.3),
            color=:green,
            lw=2,
            label=(labels ? "front (Ihat)" : false),
        )
    end

    # depth cues: CG marker shadowed on the floor / walls
    if lims !== nothing && (shade || project)
        xmin, xmax, ymin, ymax, zmin, zmax = lims
        if shade
            scatter!(
                plt, [cg[1]], [cg[2]], [zmin]; color=:gray40, markersize=3, label=false
            )
        end
        if project
            scatter!(
                plt, [cg[1]], [ymax], [cg[3]]; color=:gray40, markersize=3, label=false
            )
            scatter!(
                plt, [xmin], [cg[2]], [cg[3]]; color=:gray40, markersize=3, label=false
            )
        end
    end

    return plt
end

"""
    plot_kite(q; plt=plot(), kw...)   # single state q = (α, θ, φ, β)
    plot_kite(sol, t; kw...)          # state of `sol` at time `t`

Single-frame kite plot; see `plot_kite!` for the options.
"""
plot_kite(q; plt=plot(), kw...) = plot_kite!(plt, q; kw...)
plot_kite(sol, t; kw...) = plot_kite(state(sol)(t)[1:4]; u=control(sol)(t), kw...)

"""
    kite_paths(sol; N_pts=300)

Evaluate `sol` on `N_pts` times; returns `(; t, qs, cgs, arms, lims)` with the
CG and arm-tip positions and axis limits `(xmin, xmax, ymin, ymax, zmin, zmax)`
shared by all plot functions below.
"""
function kite_paths(sol; N_pts=300)
    tf = time_grid(sol)[end]
    t = range(0.0, tf; length=N_pts)
    qs = [state(sol)(t)[1:4] for t in t]
    cgs = [pos_CG(q) for q in qs]
    arms = [SVector(pars.larm * cos(q[1]), pars.larm * sin(q[1]), 0.0) for q in qs]
    all_pts = vcat(cgs, arms)
    pad = 2.0
    xmin, xmax = minimum(v -> v[1], all_pts) - pad, maximum(v -> v[1], all_pts) + pad
    ymin, ymax = minimum(v -> v[2], all_pts) - pad, maximum(v -> v[2], all_pts) + pad
    zmax = maximum(v -> v[3], cgs) + pad
    return (; t, qs, cgs, arms, lims=(xmin, xmax, ymin, ymax, 0.0, zmax))
end

"""
    plot_background!(plt, paths; shade=false, project=false, camera=(40,25), size=(900,750))

Faint CG trajectory + dashed arm path; sets axis limits and camera so every
plot built on it shares the same framing. `shade=true` adds the floor shadow
of the CG path (the arm already lies in the ground plane), `project=true`
adds back-wall (y=ymax) and side-wall (x=xmin) projections.
"""
function plot_background!(
    plt, paths; shade=false, project=false, camera=(40, 25), size=(900, 750)
)
    xmin, xmax, ymin, ymax, zmin, zmax = paths.lims
    n = length(paths.t)
    X = [c[1] for c in paths.cgs]
    Y = [c[2] for c in paths.cgs]
    Z = [c[3] for c in paths.cgs]
    plot!(
        plt;
        xlabel="x (m)",
        ylabel="y (m)",
        zlabel="z (m)",
        xlims=(xmin, xmax),
        ylims=(ymin, ymax),
        zlims=(zmin, zmax),
        camera=camera,
        size=size,
        legend=:topright,
    )
    if shade
        plot!(plt, X, Y, fill(zmin, n); color=:gray50, lw=1.5, label=false)
    end
    if project
        plot!(plt, X, fill(ymax, n), Z; color=:gray50, lw=1.5, label=false)
        plot!(plt, fill(xmin, n), Y, Z; color=:gray50, lw=1.5, label=false)
    end
    plot!(plt, X, Y, Z; color=:crimson, lw=2.5, label="CG")
    plot!(
        plt,
        [a[1] for a in paths.arms],
        [a[2] for a in paths.arms],
        [a[3] for a in paths.arms];
        color=:gray30,
        lw=2,
        linestyle=:dash,
        label="arm path",
    )
    return plt
end

"""
    plot_kite_positions(sol, times; shade=false, project=false, front=false,
                        camera=(40,25), size=(900,750))
    plot_kite_positions(sol, n::Int; kw...)   # n snapshots evenly spread over the period

Static plot of several kite snapshots over the shared background.
"""
function plot_kite_positions(
    sol, times; shade=false, project=false, front=false, camera=(40, 25), size=(900, 750)
)
    paths = kite_paths(sol)
    plt = plot()
    plot_background!(plt, paths; shade=shade, project=project, camera=camera, size=size)
    for (i, t) in enumerate(times)
        plot_kite!(
            plt,
            state(sol)(t)[1:4];
            labels=(i == 1),
            front=front,
            u=control(sol)(t),
            shade=shade,
            project=project,
            lims=paths.lims,
        )
    end
    return plt
end
function plot_kite_positions(sol, n::Int; kw...)
    return plot_kite_positions(
        sol, collect(range(0.0, time_grid(sol)[end]; length=n)); kw...
    )
end

"""
    animate_kite(sol; fps=30, trail_length=1, trail_step=nothing,
                 start_camera=(30,30), end_camera=(40,30), size=(900,750),
                 shade=false, project=false, front=false, filename="kite_animation.gif")

Animated GIF of the kite (arm + tether + wing panels) over the full solution
period, on the shared background, in the style of
`Visualization.animate_trajectory_4D`: fading trail behind the kite and a
camera sweeping from `start_camera` to `end_camera` across the cycle.

`fps` plays three consistent roles (same method as in `Visualization`):
- frame count = `round(tf * fps)` → 1:1 real-time playback (deliberately no
  min-points clamp: it would slow down short cycles),
- playback rate of the gif,
- trail spacing: `trail_step` defaults to `1 / (5 * fps)`, so the trail holds
  ~`5 * trail_length` samples at any fps.
"""
function animate_kite(
    sol;
    fps=30,
    trail_length=1,
    trail_step=nothing,
    start_camera=(30, 30),
    end_camera=(40, 30),
    size=(900, 750),
    shade=false,
    project=false,
    front=false,
    filename="kite_animation.gif",
)
    tf = time_grid(sol)[end]
    paths = kite_paths(sol)
    t_anim = range(0.0, tf; length=max(2, round(Int, tf * fps)))   # real-time: frames = period × fps
    trail_step = trail_step === nothing ? 1 / (5 * fps) : trail_step

    anim = @animate for t_ in t_anim
        τ = t_ / tf
        camera = start_camera .* (1 - τ) .+ end_camera .* τ
        t_trail = [max(0.0, t_ - trail_length):trail_step:t_; t_]
        trail_cgs = [pos_CG(state(sol)(t)[1:4]) for t in t_trail]
        width = range(0.0, 1.0; length=length(t_trail))

        plt = plot()
        plot_background!(plt, paths; shade=shade, project=project, camera=camera, size=size)
        plot!(
            plt,
            [c[1] for c in trail_cgs],
            [c[2] for c in trail_cgs],
            [c[3] for c in trail_cgs];
            msw=0,
            lw=4 .* width,
            alpha=width,
            color=:crimson,
            label=false,
        )
        plot_kite!(
            plt,
            state(sol)(t_)[1:4];
            labels=false,
            front=front,
            u=control(sol)(t_),
            shade=shade,
            project=project,
            lims=paths.lims,
        )
        plot!(plt; title="t = $(round(t_, digits=2)) s / $(round(tf, digits=2)) s")
    end fps = fps
    gif(anim, filename; fps=fps)
    println(
        "Animation saved to $filename ($(length(t_anim)) frames = real-time at $fps fps)"
    )
    return anim
end
