#=
10D
 - Plot trajectory 10D

4D
 - Plot trajectory 4D
 - Animate trajectory 4D
 - Plot states 2D (eight + circle)
 - plot states 2D (phase space)
=#

module Visualisation

using Plots
using SplitApplyCombine
using StrFormat
using LaTeXStrings

import KEEP.PointMass10 as PM10
import KEEP.PointMass4 as PM4
using KEEP.PointMassPara
using KEEP.TorqueFunction: torque_function

export plot_trajectory_10D
export plot_trajectory_4D
export plot_avg_power_10D
export plot_avg_power_4D
export plot_phase_space
export plot_eight_circle
export make_links
export TICKS_PI
export TICKS_HALF_PI

const DEFAULT_PPS = 30

# For plotting phase space etc...
const TICKS_PI = (-2π:π:2π, [L"-2π", L"-π", L"0", L"π", L"2π"])
const TICKS_HALF_PI = (-2π:π/2:2π, [L"-2π", L"-3π/2", L"-π", L"-π/2", L"0", L"\frac{π}{2}", L"π", L"\frac{3π}{2}", L"2π"])
# const TICKS_HALF_PI = (-2π:π/2:2π, [L"-2π", L"-\frac{3π}{2}", L"-π", L"-\frac{π}{2}", L"0", L"\frac{π}{2}", L"π", L"\frac{3π}{2}", L"2π"])

"""
Accepts either (t0, tf) or tf only
Return a vector t and a tspan = (t0, tf)"""
function build_t(tspan, points_per_second)
    len = round(Int, (tspan[2] - tspan[1]) * points_per_second)
    len = max(300, len)  # For when plotting only one period
    return (range(tspan..., len), tspan)
end
build_t(tspan::Number, points_per_second) = build_t((0, tspan), points_per_second)

function plot_trajectory_10D(sol; tspan=extrema(sol.t), points_per_second=DEFAULT_PPS)
    p = sol.prob.p
    t, _ = build_t(tspan, points_per_second)
    q = sol.(t, idxs=1:5)
    A = [PM10.compute_posA(q, p) for q in q]  # *A* like tip of the *A*rm
    K1 = [PM10.compute_pos1(q, p) for q in q]  # *K* like tip of the *K*ite, computed from R and τ
    K2 = [PM10.compute_pos2(q, p) for q in q]  # *K* like tip of the *K*ite, computed from α, θ2 and φ2

    arm_and_lines = [[0 * A, A, K1, NaN * A] for (A, K1) in zip(A, K1)]

    plot()
    plot!(invert(flatten(arm_and_lines))..., c=:black, alpha=20 / length(t), lw=1, label="")
    plot!(invert(A)..., label="A(α)", c=1)
    plot!(invert(K1)..., label="K1(R, τ)", c=2)
    plot!(invert(K2)..., label="K2(α, θ2, φ2)", c=3)
end

function plot_trajectory_4D(sol; tspan=extrema(sol.t), points_per_second=DEFAULT_PPS)
    vbp = sol.prob.p
    t, _ = build_t(tspan, points_per_second)
    q = sol.(t, idxs=1:2)
    A = [vbp.l * PM4.compute_αhat(q[1]) for q in q]
    K = [PM4.compute_OK(PM4.compute_Rτ(q[1:2], vbp), vbp) for q in q]

    arm_and_lines = [[0 * A, A, K, NaN * A] for (A, K) in zip(A, K)]

    plot()
    plot!(invert(flatten(arm_and_lines))..., c=:black, alpha=20 / length(t), lw=1, label="")
    plot!(invert(A)..., label="A(α)", c=1)
    plot!(invert(K)..., label="K(α, τ)", c=2)
end

function _plot_avg_power(t0, tf, t, power, avg_power)
    plot(t, power, c=:black, lw=2, label="")
    plot!([t0, tf], avg_power .* [1, 1], c=:red, lw=2, label=f"Average power = \%.0f(avg_power) W")
    plot!(legend=:outerbottom, xlabel="t (s)", ylabel="Power (W)")
end

function plot_avg_power_10D(sol; tspan=extrema(sol.t), points_per_second=DEFAULT_PPS)
    dα_ind, P_ind = 8, 11
    t, (t0, tf) = build_t(tspan, points_per_second)
    p = sol.prob.p
    dα = sol.(t, idxs=dα_ind)
    power = [dα * torque_function(dα, p) for dα in dα]
    avg_power = sol(tf, idxs=P_ind) / (tf - t0)
    _plot_avg_power(t0, tf, t, power, avg_power)
end

function plot_avg_power_4D(sol; tspan=extrema(sol.t), points_per_second=DEFAULT_PPS)
    dα_ind, P_ind = 3, 5
    t, (t0, tf) = build_t(tspan, points_per_second)
    p = sol.prob.p
    L, M, T = lmt(p)
    dα = sol.(t, idxs=dα_ind)
    power = dα .* torque_function.(dα .* T, (p,)) .* (M * L^2 * T^-2)
    avg_power = sol(tf, idxs=P_ind) / (tf - t0)
    _plot_avg_power(t0, tf, t, power, avg_power)
end

function animate_trajectory_4D(sol; tspan=extrema(sol.t), fps=DEFAULT_PPS, trail_length=1, trail_step=1 / 5fps, start_camera=(30, 30), end_camera=(40, 30), figsize=(800, 800))
    vbp = sol.prob.p
    t, (t0, tf) = build_t(tspan, fps)
    global K = [PM4.compute_OK(PM4.compute_Rτ(q, vbp), vbp) for q in sol(t, idxs=1:2)]
    xmax, ymax, zmax = 1.1 * maximum.(invert(K))
    @show xmax, ymax, zmax

    @show length(t)
    @time anim = @animate for t_ in t
        t_trail = [max(0, t_ - trail_length):trail_step:t_; t_]
        q = sol(t_trail, idxs=1:2)
        α = first(q[end])
        A = vbp.l * PM4.compute_αhat(α)
        trail_K = [PM4.compute_OK(PM4.compute_Rτ(q, vbp), vbp) for q in q]
        segment = [0 * A, A, trail_K[end]]

        width = range(0.0, 1.0, length=length(t_trail))
        τ = (t_ - t0) / (tf - t0)
        camera = start_camera .* (1 - τ) .+ end_camera .* τ
        plot()
        plot(invert(segment)..., lw=1, c=:black, label="")
        plot!(invert(trail_K)..., msw=0, lw=4 * width, alpha=width, label="")
        scatter!(invert([A])..., msw=0, ms=4, label="")
        plot!(xlims=(0, xmax), ylims=(-ymax, ymax), zlims=(0, zmax),
            xlabel="x", ylabel="y", zlabel="z",
            camera=camera,
            titlefont="monospace", title=f"t = \%.2f(t_) s / \%.2f(tf) s",
            size=figsize
        )
    end fps = fps
    return anim
end

function init_plot_eight_circle(vbp)
    # Eight
    θ0 = vbp.l + vbp.r
    Δφ = θ0 * vbp.Δφ
    Δθ = vbp.Δθ * Δφ / vbp.Δφ
    plot_vbp = build_para(θ0=θ0, Δθ=Δθ, φ0=0, Δφ=Δφ)
    τ = range(-π, π, length=100)
    θ, φ = invert(PM4.τ_to_θφ.(τ, Ref(plot_vbp)))
    plot(φ, θ, aspect_ratio=:equal, color=:black, grid=false, axis=false, label="")

    # Circle
    α = range(-π, π, length=100)
    x, y = vbp.l .* invert(sincos.(α))
    plot!(y, x; color=:black, label="")

    l = plot_vbp.l
    offset = l / 4

    # order: y, x, φ, θ
    x_axes = [-l, -l - offset, plot_vbp.φ0 - plot_vbp.Δφ, plot_vbp.φ0 - plot_vbp.Δφ - offset]
    y_axes = [-l - offset, -l, plot_vbp.θ0 - plot_vbp.Δθ - offset, plot_vbp.θ0 - plot_vbp.Δθ]

    u_axes = [2l, 0, 2l, 0]
    v_axes = [0, 2l, 0, 2l]
    xy_axes = hcat(x_axes, y_axes)
    uv_axes = hcat(u_axes, v_axes)
    quiver!(eachcol(xy_axes)..., quiver=Tuple(eachcol(uv_axes)), c=:black)

    x_lbl = @. x_axes + u_axes / 2 - v_axes / 4
    y_lbl = @. y_axes + v_axes / 2 - u_axes / 4
    labels = text.(["y", "x", "φ", "θ"], :black, :center, 10)
    annotate!(x_lbl, y_lbl, labels)
    plot!([-1, -2vbp.l], alpha=0, label="")  # Else "y" is not shown

    return (plot!(), plot_vbp)
end

function add_point_plot_eight_circle!(q, plot_vbp, plot_args...; plot_kwargs...)
    α, τ = q
    θ, φ = PM4.τ_to_θφ(τ, plot_vbp)
    x, y, _ = plot_vbp.l * PM4.compute_αhat(α)
    return plot!([0, y, φ], [0, x, θ]; plot_args..., plot_kwargs...)
end

#=
TODO:
 - make param plot
 - add point(q, plot_vbp) and add points(qs, plot_vbp)
 - draw eight_circle
=#

function plot_eight_circle(qs, vbp; plot_kwargs=())
    plot_kwargs_arr = ifelse(
        isa(plot_kwargs, AbstractVector),
        plot_kwargs,
        fill(plot_kwargs, length(qs))
    )

    fig, plot_vbp = init_plot_eight_circle(vbp)
    for (q, plot_kwargs) in zip(qs, plot_kwargs_arr)
        add_point_plot_eight_circle!(q, plot_vbp; plot_kwargs...)
    end

    return fig
end


"""
Make the links given by q1s and q2s

If only qs is given, then q1 and q2 alternate in qs"""
function make_links(q1s, q2s)
    return flatten([[q1, q2, NaN * q1] for (q1, q2) in zip(q1s, q2s)])
end
make_links(qs) = make_links(qs[1:2:end], qs[2:2:end])

"""
Plot states on 2d phase space with indications on arm and kite positions
Use make_links to link the symmetric states"""
function plot_phase_space(qs; links=[[]])
    rect(x1, x2, y1, y2) = Shape(π * [x1, x1, x2, x2], π * [y1, y2, y2, y1])

    ticks = TICKS_PI

    plot(xlabel="α", ylabel="τ", title="Equilibriums linked to their symmetries", legend=:outerright, xticks=ticks, yticks=ticks)

    alpha_col = 0.1
    alpha_pattern = 0.1
    fill_up = :/
    fill_down = nothing
    alpha_up = isnothing(fill_up) ? 0 : alpha_pattern
    alpha_down = isnothing(fill_down) ? 0 : alpha_pattern
    plot!(rect(-0.5, 0.5, -1, 1), c=:green, alpha=alpha_col, label="Arm towards the wind")
    plot!(rect(-1, -0.5, -1, 1), c=:red, alpha=alpha_col, label="Arm against the wind")
    plot!(rect(0.5, 1, -1, 1), c=:red, alpha=alpha_col, label="")

    scatter!(invert(qs)..., c=:black, label="Equlibrium")
    plot!(rect(-1, 1, -1, -0.5), c=:black, alpha=alpha_up, fillstyle=fill_up, label="Kite above center of 8")
    plot!(rect(-1, 1, -0.5, 0), c=:black, alpha=alpha_down, fillstyle=fill_down, label="Kite below center of 8")
    plot!(rect(-1, 1, 0, 0.5), c=:black, alpha=alpha_up, fillstyle=fill_up, label="")
    plot!(rect(-1, 1, 0.5, 1), c=:black, alpha=alpha_down, fillstyle=fill_down, label="")

    ## We could use vspan/hspan but there is padding between the plot area and its axis which is ugly
    # vspan!(π * [-.5, .5], c=:green, alpha=alpha_col, label="Arm towards the wind")
    # vspan!(π * [-1, -.5, .5, 1], c=:red, alpha=alpha_col, label="Arm against the wind")

    # hspan!(π * [-1, -.5, 0, .5], c=:black, alpha=alpha_up, fillstyle=fill_up, label="Kite above center of 8")
    # hspan!(π * [-.5, 0, .5, 1], c=:black, alpha=alpha_down, fillstyle=fill_down, label="Kite below center of 8")

    return plot!(invert(links)..., c=:black, alpha=0.1, lw=5, label="")
end

end  # module