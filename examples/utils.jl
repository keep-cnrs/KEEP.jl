module ExampleUtils

using Plots
using StrFormat
using KEEP.TorqueFunction
import KEEP.PointMass4 as PM4
import KEEP.PointMass10 as PM10
using Logging

export should_verbose, advanced_torque_func, plot_position10, plot_position4, plot_avg_power10, plot_avg_power4, animate_position4
export vecvec2mat

function should_verbose()
    return Logging.min_enabled_level(current_logger()) <= Info
end

"""
Quadratic decay outside [-Ωmax, Ωmax] instead of instantaneous
"""
function advanced_torque_func(dα, p)
    (; Ωmin, Ωmax, Ωlim, Cmax) = p
    abs_out = begin
        if Ωmin < abs(dα) < Ωmax
            (Cmax / (Ωmax - Ωmin)) * (abs(dα) - Ωmin)
        elseif Ωmax < abs(dα) < Ωlim
            Cmax * ((abs(dα) - Ωlim) / (Ωmax - Ωlim))^10
        else
            zero(dα)
        end
    end
    return sign(dα) * abs_out
end

function vecvec2mat(vecvec)
    return reduce(hcat, vecvec)
end

function plot_position10(sol, p; t=sol.t, figure=plot())
    q = sol.(t; idxs=1:5)
    A = PM10.compute_posA.(q, (p,))
    K1 = PM10.compute_pos1.(q, (p,))
    K2 = PM10.compute_pos2.(q, (p,))

    # indices : x/y/z, O/A/K2, indice du segment
    segments = stack(vecvec2mat.((zero(A), A, K2)); dims=2)

    plot!(figure)
    plot!(xlabel=raw"$x$", ylabel=raw"$y$", zlabel=raw"$z$", legend=:topleft, size=(600, 600))
    plot!(eachslice(segments; dims=1)..., c=:black, alpha=20 / length(t), label="")
    plot!(Tuple.(A), lw=1, label=raw"$A(α)$", c=1)
    plot!(Tuple.(K1), lw=1, label=raw"$K_1(R, τ)$", c=2)
    plot!(Tuple.(K2), lw=1, label=raw"$K_2(α, θ2, φ2)$", c=3)
    return plot!()
end


# TODO: Have one function that initializes the plot, and one function that adds a position to it, then loop over all positions
# Reuse these for the animation
function plot_position4(sol, p; t=sol.t, figure=plot())
    q = sol.(t; idxs=1:5)
    α = first.(q)
    A = p.l .* PM4.compute_αhat.(α)
    K = PM4.compute_OK.(PM4.compute_Rτ.(q, (p,)), (p,))

    # indices : x/y/z, O/A/K2, indice du segment
    segments = stack(vecvec2mat.((zero(A), A, K)); dims=2)
    plot!(figure)
    plot!(xlabel=raw"$x$", ylabel=raw"$y$", zlabel=raw"$z$", legend=:topleft, size=(600, 600))
    plot!(eachslice(segments; dims=1)..., c=:black, alpha=20 / length(t), label="")
    plot!(Tuple.(A), label=raw"$A(α)$", c=1)
    plot!(Tuple.(K), label=raw"$K(α, τ)$", c=2)
    return plot!()
end


function plot_avg_power10(sol, p; t=sol.t, figure=plot())
    dα_ind, P_ind = 8, 11
    t0, tf = extrema(t)

    dα = sol.(t, idxs=dα_ind)
    power = dα .* torque_function.(dα, (p,))
    avg_power = sol(tf, idxs=P_ind) / (tf - t0)
    plot!(figure, t, power, c=:black, lw=2, label="")
    plot!([0, tf], avg_power .* [1, 1], c=:red, lw=2, label=f"Average power = \%.0f(avg_power) W")
    plot!(legend=:outerbottom, xlabel="t (s)", ylabel="Power (W)")
end


function plot_avg_power4(sol, p; t=sol.t, figure=plot(), kwargs...)
    dα_ind, P_ind = 3, 5
    L, M, T = lmt(p)
    t0, tf = extrema(t)

    dα = sol.(t, idxs=dα_ind)
    power = dα .* torque_function.(dα .* T, (p,)) .* (M * L^2 * T^-2)
    avg_power = sol(tf, idxs=P_ind) / (tf - t0)
    plot!(figure, t, power, c=:black, lw=2, label="")
    plot!([0, tf], avg_power .* [1, 1], c=:red, lw=2, label=f"Average power = \%.0f(avg_power) W")
    default_kwargs = (legend=:outerbottom, xlabel="t (s)", ylabel="Power (W)")
    return plot!(; merge(default_kwargs, kwargs)...)
end


function animate_position4(sol, p; tspan=extrema(sol.t), fps=30, trail_length=1, trail_step=1 / 5fps, start_camera=(30, 30), end_camera=(40, 30), figsize=(800, 800))
    # Compute box
    K = PM4.compute_OK.(PM4.compute_Rτ.(sol.u, (p,)), (p,))
    xmax, ymax, zmax = 1.1 * maximum(K)

    t0, tf = tspan
    anim = @animate for t_ in t0:(1/fps):tf
        t_trail = [max(0, t_ - trail_length):trail_step:t_; t_]
        q = sol.(t_trail; idxs=1:5)
        α = first(q[end])
        A = p.l * PM4.compute_αhat(α)
        trail_K = PM4.compute_OK.(PM4.compute_Rτ.(q, (p,)), (p,))
        segment = [zero(A), A, trail_K[end]] |> vecvec2mat

        width = range(0.0, 1.0, length=length(t_trail))
        τ = (t_ - t0) / (tf - t0)
        camera = start_camera .* (1 - τ) .+ end_camera .* τ
        plot(eachrow(segment)..., lw=1, c=:black, label="")
        plot!(Tuple.(trail_K), msw=0, lw=4 * width, alpha=width, label="")
        scatter!(Tuple.([A]), msw=0, ms=4, label="")
        plot!(xlims=(0, xmax), ylims=(-ymax, ymax), zlims=(0, zmax),
            xlabel="x", ylabel="y", zlabel="z",
            camera=camera,
            titlefont="monospace", title=f"t = \%.2f(t_) s / \%.2f(tf) s",
            size=figsize
        )
    end fps = fps
    return anim
end

end # module
