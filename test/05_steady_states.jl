using Setfield
using Test
using Plots

import KEEP.PointMass4 as PM4
using KEEP.PointMassPara
using KEEP.SteadyState
using KEEP.Visualisation


default(formatter=:plain, label="")
tol = 1e-12
vbp = build_vbpara()

## Find one steady state
α_init = 1
τ_init = 0.5
ss2 = steady_state(α_init, τ_init, vbp; tol=tol)
ss1 = steady_state(τ_init, vbp; tol=tol)
@test all(ddq_partial(ss2, vbp) .< tol)  # converges with (α, τ)
@test all(ddq_partial(ss1, vbp) .< tol)  # converges with (τ) only

aligned, opposite = aligned_and_opposite_steady_states(τ_init, vbp; tol=tol)
@test -π / 2 < aligned[1] < π / 2  # arm is toward positive x
@test !(-π / 2 < opposite[1] < π / 2)  # arm is toward negative x
@test all(ddq_partial(aligned, vbp) .< tol)  # both are steady states
@test all(ddq_partial(opposite, vbp) .< tol)

## Finding all steady states

ss = all_steady_states(vbp, tol=tol)


println("\n05_steady_states.jl@18: The differentiation here is broken, running on MacOS works, debugging on MacOS raises a 'DualMismatch' error from ForwardDiff, and the github worker has an error with ReverseDiff.

It happens in `all_steady_states`.")
all_steady_states(@set vbp[[:r, :Δφ]] = [85, 1.464945])  # Does not throw
all_steady_states(@set vbp[[:r, :Δφ]] = [85, 1.464944])  # Throws method error

@test all([maximum(abs, ddq_partial(ss_i, vbp)) for ss_i in ss] .< tol)  # All steady states

ss_halton = all_steady_states_halton(vbp, tol=tol)
@test all([maximum(abs, ddq_partial(ss_i, vbp)) for ss_i in ss_halton] .< tol)  # All steady states

@test ss ≈ ss_halton  # Same equlibriums in same order are found with both methods
@test length(ss) == 8  # There *should* be 8 steady states, but it depends on the parameters

## Plots

links = make_links(ss)
plot_phase_space(ss; links=links)

Ax_pos = [-π / 2 < ss_i[1] < π / 2 for ss_i in ss]
plot_eight_circle(ss[Ax_pos], vbp; plot_kwargs=(lw=3, c=:green))
plot_eight_circle(ss[.!Ax_pos], vbp; plot_kwargs=(lw=2, c=:red, ls=:dash))

## Animations

if false
    anim_grouped = @animate for (ss, plot_kwargs, title) in zip([ss[Ax_pos], ss[.!Ax_pos]], [(lw=3, c=:green,), (lw=3, c=:red,)], ["Arm not facing the wind", "Arm facing the wind"])
        plot_eight_circle(ss, vbp; plot_kwargs)
        plot!(title=title)
    end fps = 0.5
    display(gif(anim_grouped, fps=0.5))

    anim_each = @animate for (i, ss) in enumerate(ss)
        Ax_pos = -π / 2 < ss[1] < π / 2
        plot_eight_circle([ss], vbp; plot_kwargs=(lw=3, c=ifelse(Ax_pos, :green, :red)))
        plot!(title="Steady states #$i")
    end fps = 1
    display(gif(anim_each, fps=1))
end