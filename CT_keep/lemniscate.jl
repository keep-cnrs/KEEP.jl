using Pkg
Pkg.activate("CT_keep")

using StaticArrays
using ForwardDiff: derivative
using OptimalControl, NLPModelsIpopt
using ComponentArrays: ComponentArray as CA
using Plots

using KEEP.PointMassPara: build_vbpara, lmt
using KEEP.PointMass4: dynamics
using KEEP.TorqueFunction: torque_function
using KEEP.LimitCycle: compute_limit_cycle
using KEEP: TAU0

const vbp0 = build_vbpara()
const syms = (:r, :I_eq, :torque_slope)

begin  # Functions
    f(x, p) = begin
        T = promote_type(eltype(vbp0), eltype(p))
        vbp = CA(T.(vbp0); (syms .=> p)...)
        x_dyn = vcat(x, 0)
        return dynamics(x_dyn, vbp)[SOneTo(4)]
    end

    generated_power(dα, p) = begin
        T = promote_type(eltype(vbp0), eltype(p))
        vbp = CA(T.(vbp0); (syms .=> p)...)
        L, M, T = lmt(vbp)
        dα_adim = dα * T
        return dα_adim * torque_function(dα_adim, vbp) * M * L^2 * T^-3
    end

    function calc_vars(lc)
        tf = last(lc.t)
        p = lc.prob.p[syms]
        return [tf, p...]
    end

    f(rand(4), [15, 50, 40])
    generated_power(rand(), [15, 50, 40])

    function build_lc_ocp(; tf_init=2.8)
        ocp_lc = @def begin
            vars ∈ R², variable  # tf, dummy
            tf = vars[1]
            t ∈ [0, tf], time
            X = (α, τ, dα, dτ) ∈ R⁴, state

            # Phase fixing
            τ(0) == TAU0

            # Cyclic condition
            α(tf) - α(0) == 0
            τ(tf) - τ(0) == 2π
            dα(tf) - dα(0) == 0
            dτ(tf) - dτ(0) == 0

            Ẋ(t) == f(X(t), collect(vbp0[syms]))

            vars[2]^2 → min
        end

        estimate_τ(t, tf) = 2π * t / tf + TAU0
        estimate_α(t, tf) = sin(estimate_τ(t, tf) - π / 4)
        estimate_dα(t, tf) = derivative(t_ -> estimate_α(t_, tf), t)
        estimate_dτ(t, tf) = derivative(t_ -> estimate_τ(t_, tf), t)

        init_lc = @init ocp_lc begin
            vars := [tf_init, 0]
            α(t) := estimate_α(t, tf_init)
            τ(t) := estimate_τ(t, tf_init)
            dα(t) := estimate_dα(t, tf_init)
            dτ(t) := estimate_dτ(t, tf_init)
        end

        return ocp_lc, init_lc
    end

    function build_optim_ocp(lc; factor)
        ocp = @def begin
            # Initialiser variable, state, control avec des vecteurs vides et ajouter au fur et a mesure ?
            # tf ∈ R, variable
            # p ∈ R³, variable
            vars ∈ R⁴, variable
            tf = vars[1]
            p = vars[2:4]
            t ∈ [0, tf], time
            # X = (α, τ, dα, dτ, W) ∈ R⁵, state
            X = (α, τ, dα, dτ) ∈ R⁴, state

            calc_vars(lc) ./ factor <= vars <= calc_vars(lc) .* factor

            τ(0) == TAU0
            α(tf) - α(0) == 0
            τ(tf) - τ(0) == 2π
            dα(tf) - dα(0) == 0
            dτ(tf) - dτ(0) == 0
            # W(0) == 0

            # Meilleur message d'erreur pour = au lieu de ==: did you mean == ?
            # X \dot et pas \dot X
            Ẋ(t) == f(X(t), p)

            # alias `->` to `→`
            ∫(generated_power(dα(t), p) / tf) → max
            # ∫(generated_power(dα(t), p) / tf) → max
        end

        init = @init ocp begin
            vars := calc_vars(lc)
            X(t) := lc(t)[1:4]
        end
        return ocp, init
    end

    function build_scaled_optim_ocp(lc; factor)
        ocp = @def begin
            vars ∈ R⁴, variable
            tf = vars[1]
            p = vars[2:4]
            s ∈ [0, 1], time                      # Fixed time interval [0, 1]
            X = (α, τ, dα, dτ) ∈ R⁴, state

            calc_vars(lc) ./ factor <= vars <= calc_vars(lc) .* factor

            τ(0) == TAU0
            α(1) - α(0) == 0
            τ(1) - τ(0) == 2π
            dα(1) - dα(0) == 0
            dτ(1) - dτ(0) == 0

            Ẋ(s) == tf * f(X(s), p)               # Explicitly scaled dynamics

            ∫(generated_power(dα(s), p)) → max     # Directly yields P_avg
        end

        # Interpolate the initial guess from [0, tf] to [0, 1]
        init = @init ocp begin
            vars := calc_vars(lc)
            X(s) := lc(s * calc_vars(lc)[1])[1:4]
        end
        #=

        =#
        return ocp, init
    end

    function solve_ocp(ocp, init; kwargs...)
        sol_init = solve(ocp; kwargs..., init=init, max_iter=0)
        sol = solve(ocp; kwargs..., init=init)
        return sol_init, sol
    end

    function make_var_plot(sol, sol_init)
        x_ticks = 1:length(syms)+1
        ratios = variable(sol) ./ variable(sol_init)
        x_labels = (x_ticks, (:tf, syms...))

        plot(
            yscale=:log10,
            xticks=x_labels,
            ylabel="Ratio (Optimized / Initial)",
            title="Parameter Optimization Results",
            framestyle=:box,
            legend=:outertopright
        )

        hline!([1 / factor], fillrange=[factor], color=:grey, alpha=0.15, label="Allowed Range")

        hline!([1], c=:black, ls=:dash, lw=1.5, label="Baseline (1.0)")

        for i in x_ticks
            plot!([i, i], [1.0, ratios[i]], c=:grey, lw=1.0, label="")
        end

        scatter!(x_ticks, ones(length(x_ticks)),
            m=:circle, mc=:white, msc=:black, ms=6, label="Initial")

        scatter!(x_ticks, ratios,
            m=:circle, mc=:crimson, ms=7, label="Optimized")

        ylims!(1 / factor / 1.2, factor * 1.2)

        return plot!()
    end

    function make_plots(sol_init, sol)
        plot(sol_init)
        state_plot = plot!(sol)
        var_plot = make_var_plot(sol, sol_init)
        return state_plot, var_plot
    end
end

oc_kwargs = (grid_size=30, backend=:generic)

lc = compute_limit_cycle(vbp0; sense=+, save_everystep=true);

factor = 5
lc_ocp, lc_init = build_lc_ocp();
optim_ocp, optim_init = build_optim_ocp(lc; factor=factor);
# scaled_optim_ocp, scaled_optim_init = build_scaled_optim_ocp(lc; factor=factor);

lc_init_sol, lc_sol = solve_ocp(lc_ocp, lc_init; oc_kwargs...);
optim_init_sol, optim_sol = solve_ocp(optim_ocp, optim_init; oc_kwargs...);
# scaled_optim_init_sol, scaled_optim_sol = solve_ocp(scaled_optim_ocp, scaled_optim_init; oc_kwargs...);


display.(make_plots(optim_init_sol, optim_sol))
# display.(make_plots(scaled_optim_init_sol, scaled_optim_sol))

using KEEP.PointMassPara: build_para
using KEEP.Optimization: optimize
L, M, T = lmt(vbp0)
p0 = build_para(vbp0)
units = p0[syms] ./ vbp0[syms]
lb = p0[syms] ./ factor
ub = p0[syms] .* factor
solution, stats, model = optimize(p0, syms, lb, ub)
solution.params ./ variable(optim_sol)[2:end] ./ units  # Optimum is OK

using KEEP.LimitCycle: shoot as lc_shoot

vbp = build_vbpara(CA(p0; solution.params...))
shooting = solution[1:4]
solution_sim = lc_shoot(shooting, vbp, save_everystep=true)


###########
# Same values
###########


P_keep = solution.P
P_ct = optim_sol.objective


##########
## Indirect single shooting
##########
using OrdinaryDiffEq

const flow = Flow(optim_ocp);

# https://control-toolbox.org/MagneticResonanceImaging.jl/stable/saturation.html
function shoot!(res, x0, p0, tf, λ)
    xf, pf, pλf = flow(0, x0, p0, tf, [tf, λ...]; augment=true)
    Hf = pf' * f(xf, λ) - generated_power(xf[3], λ) / tf

    res[1] = x0[2] - TAU0
    res[2:5] = xf - x0 - [0, 2π, 0, 0]
    res[6] = Hf
    res[7:9] = (pf-p0)[[1, 3, 4]]
    res[10:12] = pλf[2:4]
    return res
end

@views function shoot!(res, x0, p0, xf, pf, tf, λ, pλf)
    res[1] = x0[2] - TAU0
    res[2:5] = xf - x0 - [0, 2π, 0, 0]
    res[6] = pf' * f(xf, λ) - generated_power(xf[3], λ) / tf # Hf
    res[7:9] = (pf-p0)[[1, 3, 4]]
    res[10:12] = pλf
    return res
end

@views function shoot!(res, x0, p0, tf, λ)
    # λ is previous p
    xf, pf, pvf = flow(0.0, x0, p0, tf, [tf; λ]; augment=true)
    pλf = pvf[2:4]
    return shoot!(res, x0, p0, xf, pf, tf, λ, pλf)
end

@views shoot!(res, y, _) = shoot!(res, y[1:4], y[5:8], y[9], y[10:12])
shoot(y) = shoot!(similar(y), y, nothing)

# Extract solution from direct method for initialization
p_direct = costate(optim_sol)
x_direct = state(optim_sol)
tf_direct, λ_direct... = variable(optim_sol)

# Initial guess
x0_guess = x_direct(0)
p0_guess = p_direct(0)
tf_guess, λ_guess = tf_direct, λ_direct

y_guess = [x0_guess; p0_guess; tf_guess; λ_guess]
shoot(y_guess)
# NLE problem with initial guess (2 unknowns: p0, λ)
if false  # To implement
    prob_indirect = NonlinearProblem(shoot!,)

    # Solve shooting equations
    shooting_sol = solve(prob_indirect; show_trace=Val(true))
    p0_sol, λ_sol = shooting_sol.u

    println("Indirect solution:")
    println("Initial costate: p0 = ", p0_sol)
    println("Parameter: λ = ", λ_sol)
end

using ForwardDiff
using KEEP.PointMass4: integrate
using LinearAlgebra
if false
    begin  # solve_trajectory_with_jacobian
        """
            solve_trajectory_with_jacobian(u0, tf, vbp)

        Integrates the system with ForwardDiff.Dual numbers seeded at `u0`.
        Returns two functions:
        - `f(t)`: returns the N-element state vector at time t.
        - `J(t)`: returns the N×N Jacobian matrix (∂u(t)/∂u0) at time t.
        """
        function solve_trajectory_with_jacobian(u0, tf, vbp)
            N = length(u0)
            T = eltype(u0)
            tag = ForwardDiff.Tag(integrate, T)

            # 1. Seed the initial state with Dual numbers (identity matrix partials)
            u0_dual = SVector{N}(
                ForwardDiff.Dual{typeof(tag)}(
                    u0[i],
                    ForwardDiff.Partials(ntuple(j -> j == i ? one(T) : zero(T), N))
                ) for i in 1:N
            )

            # 2. Integrate using the dual-valued initial state
            ode_sol_dual = integrate(u0_dual, tf, vbp; save_everystep=true)

            # 3. Define the continuous-time lookup functions using the ODE interpolation
            f(t) = ForwardDiff.value.(ode_sol_dual(t))

            J(t) =
                let ut = ode_sol_dual(t)
                    SMatrix{N,N}(ForwardDiff.partials(ut[i], j) for i in 1:N, j in 1:N)
                end

            return f, J
        end


        """
            solve_trajectory_with_jacobian_4d(u0, tf, vbp)

        Integrates the 5D system, but only tracks the sensitivity of the first 4 states 
        with respect to the first 4 initial conditions.

        Returns:
          - `f(t)`: returns the 4-element state vector at time t.
          - `J(t)`: returns the 4×4 Jacobian matrix (∂u_{1:4}(t)/∂u0_{1:4}) at time t.
        """
        function solve_trajectory_with_jacobian_4d(u0, tf, vbp)
            T = eltype(u0)
            tag = ForwardDiff.Tag(integrate, T)

            # We differentiate with respect to the first 4 states
            N_diff = 4

            # Build a 5-element SVector of Duals, but with only 4 derivative partials
            u0_dual = SVector{5}(
                i <= N_diff ?
                ForwardDiff.Dual{typeof(tag)}(u0[i], ForwardDiff.Partials(ntuple(j -> j == i ? one(T) : zero(T), N_diff))) :
                ForwardDiff.Dual{typeof(tag)}(u0[i], ForwardDiff.Partials(ntuple(j -> zero(T), N_diff)))
                for i in 1:5
            )

            # Integrate the 5D system
            ode_sol_dual = integrate(u0_dual, tf, vbp; save_everystep=true)

            # f(t) returns the first 4 nominal states
            f(t) = SVector{4}(ForwardDiff.value(ode_sol_dual(t)[i]) for i in 1:4)

            # J(t) extracts the 4×4 Jacobian for the first 4 states
            J(t) =
                let ut = ode_sol_dual(t)
                    SMatrix{4,4}(ForwardDiff.partials(ut[i], j) for i in 1:4, j in 1:4)
                end

            return f, J
        end
    end

    tf, p... = variable(optim_sol)  # Float64, Vector{Float64}
    u0 = SA[vcat(x_direct(0), 0)...]  # StaticArray{Float64}
    optim_vbp = CA(vbp0; (syms .=> p)...)  # ComponentArray{Float64}

    optim_lc = compute_limit_cycle(u0, optim_vbp)
    u0_bis = first(optim_lc.u)
    tf_bis = last(optim_lc.t)

    _, J = solve_trajectory_with_jacobian_4d(u0, tf, optim_vbp)
    # _, J = solve_trajectory_with_jacobian(u0_bis, tf_bis, optim_vbp);

    # 2. Define a time grid for evaluation
    ts = J.ode_sol_dual.t

    # 3. Compute the singular values at each time point
    # svdvals(J(t)) returns a sorted vector of 5 singular values [σ₁, σ₂, ..., σ₅]
    svals = [svdvals(J(t)) for t in ts]

    # 4. Reshape the data for plotting (matrix of size: length(ts) × 5)
    svals_matrix = reduce(hcat, svals)'

    # 5. Plot the singular values on a logarithmic scale
    plot(
        ts,
        svals_matrix,
        yscale=:log10,
        xlabel="Time (t)",
        ylabel="Singular Values (log scale)",
        label=["σ₁ (Max)" "σ₂" "σ₃" "σ₄" "σ₅ (Min)"],
        title="Singular Values of J(t) over Time",
        lw=2,
        legend=:outerright
    )
    nothing

    ##########
    ## Indirect multiple shooting (General N)
    ##########

    function shoot_multiple_N!(res, y, N)
        # Extract variables
        tf = y[end-3]
        λ = y[end-2:end]

        # Grid of matching times
        ts = range(0, tf, length=N + 1)

        pv_total = zeros(eltype(y), 4)
        x_curr = y[1:4]
        p_curr = y[5:8]

        idx_res = 13 # Match conditions start after the first 12 boundary/transversality conditions

        for i in 1:N
            t_start = ts[i]
            t_end = ts[i+1]

            # Integrate segment
            x_next_flow, p_next_flow, pv_segment = flow(t_start, x_curr, p_curr, t_end, [tf, λ...]; augment=true)
            pv_total .+= pv_segment

            if i < N
                # Extract intermediate state/costate variables
                idx_var = 8 * i + 1
                x_next_var = y[idx_var:idx_var+3]
                p_next_var = y[idx_var+4:idx_var+7]

                # Formulate matching conditions
                res[idx_res:idx_res+3] = x_next_flow - x_next_var
                res[idx_res+4:idx_res+7] = p_next_flow - p_next_var
                idx_res += 8

                # Step forward
                x_curr = x_next_var
                p_curr = p_next_var
            else
                # Enforce boundary conditions on the final segment
                x_final = x_next_flow
                p_final = p_next_flow

                # Boundary & Periodic conditions
                res[1] = y[2] - TAU0                  # τ(0) == TAU0 (using x0[2])
                res[2:5] = x_final - y[1:4] - [0, 2π, 0, 0] # Cyclic conditions on state

                # Hamiltonian condition
                Hf = p_final' * f(x_final, λ) - generated_power(x_final[3], λ) / tf
                res[6] = Hf

                # Periodic conditions on costates
                res[7:9] = (p_final-y[5:8])[[1, 3, 4]]

                # Transversality conditions on parameters
                res[10:12] = pv_total[2:4]
            end
        end
        return res
    end

    # Generate initial guess by evaluating the direct solution at N points
    function generate_initial_guess_N(x_direct, p_direct, tf_direct, λ_direct, N)
        ts = range(0, tf_direct, length=N + 1)
        y_guess = Float64[]
        for i in 1:N
            t = ts[i]
            append!(y_guess, x_direct(t))
            append!(y_guess, p_direct(t))
        end
        push!(y_guess, tf_direct)
        append!(y_guess, λ_direct)
        return y_guess
    end

    # Set segment count
    N_segments = 10

    # Wrappers
    shoot_multN!(res, y, _) = shoot_multiple_N!(res, y, N_segments)
    shoot_multN(y) = shoot_multiple_N!(similar(y), y, N_segments)

    y_guess_multN = generate_initial_guess_N(x_direct, p_direct, tf_direct, λ_direct, N_segments)

    # Verify initial residual
    println("N=$N_segments residual norm: ", sum(abs2, shoot_multN(y_guess_multN)))

    using NonlinearSolve

    # 1. Define the shooting function matching the NonlinearSolve signature: f!(res, y, p)
    # Here, the third argument `N` is passed as the parameter of the NonlinearProblem.
    function shoot_mult_N!(res, y, N)
        return shoot_multiple_N!(res, y, N)
    end

    # 2. Select the number of segments and generate the initial guess
    N_segments = 10

    y_guess_multN = generate_initial_guess_N(
        x_direct,
        p_direct,
        tf_direct,
        λ_direct,
        N_segments
    )

    # 3. Define and solve the Nonlinear Problem
    # We pass N_segments as the third argument to parameterize the problem.
    prob_indirect_N = NonlinearProblem(shoot_mult_N!, y_guess_multN, N_segments)

    # Solve using a Trust Region method (or NewtonRaphson())
    sol_indirect_N = solve(prob_indirect_N, TrustRegion(); show_trace=Val(true))

    # 4. Extract results
    y_sol = sol_indirect_N.u
    tf_sol = y_sol[end-3]
    λ_sol = y_sol[end-2:end]

    println("\n--- Optimization Results ($N_segments segments) ---")
    println("Sol Status: ", sol_indirect_N.retcode)
    println("Final Time (tf): ", tf_sol)
    println("Parameters (λ): ", λ_sol)

    # 5. Reconstruct the continuous trajectory
    # Since multiple shooting only solves for states and costates at the grid nodes,
    # we integrate the flow forward from each node to recover the full continuous path.
    ts_nodes = range(0, tf_sol, length=N_segments + 1)
    points_per_segment = 50

    t_plot = Float64[]
    x_plot = Vector{Float64}[]
    p_plot = Vector{Float64}[]

    for i in 1:N_segments
        t_start = ts_nodes[i]
        t_end = ts_nodes[i+1]

        # Retrieve state and costate values at node i-1
        idx_var = 8 * (i - 1) + 1
        x_node = y_sol[idx_var:idx_var+3]
        p_node = y_sol[idx_var+4:idx_var+7]

        # Generate fine time grid for the current segment
        ts_seg = range(t_start, t_end, length=points_per_segment)

        for t in ts_seg
            # Propagate from t_start to current t
            xt, pt, _ = flow(t_start, x_node, p_node, t, [tf_sol, λ_sol...]; augment=true)
            push!(t_plot, t)
            push!(x_plot, xt)
            push!(p_plot, pt)
        end
    end

    # Convert trajectory lists to matrices for easy plotting
    x_matrix = hcat(x_plot...)'  # Dimensions: (Time steps) x 4
    p_matrix = hcat(p_plot...)'  # Dimensions: (Time steps) x 4

    # Example verification plot
    using Plots
    plot(t_plot, x_matrix, label=["α" "τ" "dα" "dτ"], lw=1.5,
        xlabel="Time [s]", ylabel="States", title="Reconstructed State Trajectories")



















    ## 

    using NonlinearSolve
    using Plots

    """
        shoot_multiple_N!(res, y, N)

    Formulate the multiple shooting equations for N segments.
    y: [ x_node_1, p_node_1, ..., x_node_N, p_node_N, tf, λ... ]
    """
    function shoot_multiple_N!(res, y, N)
        # 1. Unpack decision variables using structured views
        nodes = reshape(view(y, 1:8N), 8, N)
        tf = y[end-3]
        λ = y[end-2:end]

        # Grid of matching times
        ts = range(0, tf, length=N + 1)

        # Views for the matching conditions (8 equations for each of the N-1 internal nodes)
        matching_res = reshape(view(res, 13:8N+4), 8, N - 1)
        pv_total = zeros(eltype(y), 4)

        # Start integration at the first node
        x_curr, p_curr = nodes[1:4, 1], nodes[5:8, 1]

        for i in 1:N
            # Integrate the current segment
            x_flow, p_flow, pv_segment = flow(ts[i], x_curr, p_curr, ts[i+1], [tf, λ...]; augment=true)
            pv_total .+= pv_segment

            if i < N
                # Extract variables for the next node
                x_next, p_next = nodes[1:4, i+1], nodes[5:8, i+1]

                # Enforce matching conditions at the node interface
                matching_res[1:4, i] .= x_flow .- x_next
                matching_res[5:8, i] .= p_flow .- p_next

                # Step forward to the next segment
                x_curr, p_curr = x_next, p_next
            else
                # Final segment: boundary & transversality conditions
                x0, p0 = nodes[1:4, 1], nodes[5:8, 1]

                res[1] = x0[2] - TAU0                       # Initial state constraint τ(0) == TAU0
                res[2:5] .= x_flow .- x0 .- [0, 2π, 0, 0]   # Cyclic state conditions

                # Hamiltonian condition (tf is free)
                Hf = p_flow' * f(x_flow, λ) - generated_power(x_flow[3], λ) / tf
                res[6] = Hf

                # Cyclic costate conditions
                res[7:9] .= (p_flow.-p0)[[1, 3, 4]]

                # Transversality conditions on parameters λ
                res[10:12] .= pv_total[2:4]
            end
        end
        return res
    end

    # Generate initial guess using concise list comprehensions
    function generate_initial_guess_N(x_direct, p_direct, tf_direct, λ_direct, N)
        ts = range(0, tf_direct, length=N + 1)[1:N]
        nodes_guess = [vcat(x_direct(t), p_direct(t)) for t in ts]
        return vcat(nodes_guess..., tf_direct, λ_direct)
    end

    # 1. Define the shooting function with NonlinearSolve's signature: f!(res, y, p)
    shoot_mult_N!(res, y, N) = shoot_multiple_N!(res, y, N)

    # 2. Select segment count and generate the initial guess
    N_segments = 10
    y_guess_multN = generate_initial_guess_N(x_direct, p_direct, tf_direct, λ_direct, N_segments)

    # 3. Define and solve the Nonlinear Problem
    prob_indirect_N = NonlinearProblem(shoot_mult_N!, y_guess_multN, N_segments)
    sol_indirect_N = solve(prob_indirect_N; show_trace=Val(true), maxiters=2)

    # 4. Extract results
    y_sol = sol_indirect_N.u
    tf_sol = y_sol[end-3]
    λ_sol = y_sol[end-2:end]

    println("\n--- Optimization Results ($N_segments segments) ---")
    println("Sol Status: ", sol_indirect_N.retcode)
    println("Final Time (tf): ", tf_sol)
    println("Parameters (λ): ", λ_sol)

    # 5. Reconstruct the continuous trajectory
    ts_nodes = range(0, tf_sol, length=N_segments + 1)
    points_per_segment = 50
    nodes_sol = reshape(view(y_sol, 1:8*N_segments), 8, N_segments)

    t_plot = Float64[]
    x_list = Vector{Float64}[]
    p_list = Vector{Float64}[]

    for i in 1:N_segments
        ts_seg = range(ts_nodes[i], ts_nodes[i+1], length=points_per_segment)
        x_node, p_node = nodes_sol[1:4, i], nodes_sol[5:8, i]

        for t in ts_seg
            xt, pt, _ = flow(ts_nodes[i], x_node, p_node, t, [tf_sol, λ_sol...]; augment=true)
            push!(t_plot, t)
            push!(x_list, xt)
            push!(p_list, pt)
        end
    end

    x_matrix = hcat(x_list...)'
    p_matrix = hcat(p_list...)'

    # Plot the reconstructed states
    plot(t_plot, x_matrix, label=["α" "τ" "dα" "dτ"], lw=1.5,
        xlabel="Time [s]", ylabel="States", title="Reconstructed State Trajectories")
end
##
struct MultipleShooting{F,B,K}
    flow::F
    boundary_conds!::B
    dim_x::Int
    dim_p::Int
    dim_param::Int
    kwargs::K
end

# Outer constructor: Infers dimensions from prototype inputs
function MultipleShooting(flow, boundary_conds!, x_proto, p_proto, λ_proto; kwargs...)
    return MultipleShooting(
        flow,
        boundary_conds!,
        length(x_proto),
        length(p_proto),
        length(λ_proto),
        kwargs
    )
end

# The callable behavior matching the NonlinearSolve signature: f!(res, y, N)
function (ms::MultipleShooting)(res, y, N)
    dim_node = ms.dim_x + ms.dim_p
    dim_bc = dim_node + 1 + ms.dim_param

    # Unpack node matrix and parameters
    nodes = reshape(view(y, 1:dim_node*N), dim_node, N)
    tf = y[end-ms.dim_param]
    λ = view(y, (length(y)-ms.dim_param+1):length(y))

    # Define segment matching times and output views
    ts = range(0, tf, length=N + 1)
    bc_res = view(res, 1:dim_bc)
    matching_res = reshape(view(res, (dim_bc+1):length(res)), dim_node, N - 1)

    pv_total = zeros(eltype(y), ms.dim_p)
    x_curr = nodes[1:ms.dim_x, 1]
    p_curr = nodes[ms.dim_x+1:dim_node, 1]

    for i in 1:N
        # Integrate current segment using the inferred parameters
        flow_params = [tf; λ]
        x_flow, p_flow, pv_segment = ms.flow(ts[i], x_curr, p_curr, ts[i+1], flow_params; ms.kwargs...)
        pv_total .+= pv_segment

        if i < N
            x_next = nodes[1:ms.dim_x, i+1]
            p_next = nodes[ms.dim_x+1:dim_node, i+1]

            # Matching conditions at the node interface
            matching_res[1:ms.dim_x, i] .= x_flow .- x_next
            matching_res[ms.dim_x+1:dim_node, i] .= p_flow .- p_next

            # Step forward to the next segment
            x_curr, p_curr = x_next, p_next
        else
            # Final segment: delegate boundary/transversality conditions to the user callback
            x0 = nodes[1:ms.dim_x, 1]
            p0 = nodes[ms.dim_x+1:dim_node, 1]
            ms.boundary_conds!(bc_res, x0, p0, x_flow, p_flow, tf, λ, pv_total)
        end
    end
    return res
end

# Helper to generate the initial guess vector
function generate_initial_guess(ms::MultipleShooting, x_direct, p_direct, tf_direct, λ_direct, N)
    ts = range(0, tf_direct, length=N + 1)[1:N]
    nodes_guess = [vcat(x_direct(t), p_direct(t)) for t in ts]
    return vcat(nodes_guess..., tf_direct, λ_direct)
end

# Helper to reconstruct the fine continuous trajectory from a solution vector
function reconstruct_trajectory(ms::MultipleShooting, y_sol, N; points_per_segment=50)
    dim_node = ms.dim_x + ms.dim_p
    nodes_sol = reshape(view(y_sol, 1:dim_node*N), dim_node, N)
    tf_sol = y_sol[end-ms.dim_param]
    λ_sol = view(y_sol, (length(y_sol)-ms.dim_param+1):length(y_sol))

    ts_nodes = range(0, tf_sol, length=N + 1)
    t_plot, x_list, p_list = Float64[], Vector{Float64}[], Vector{Float64}[]

    for i in 1:N
        ts_seg = range(ts_nodes[i], ts_nodes[i+1], length=points_per_segment)
        x_node, p_node = nodes_sol[1:ms.dim_x, i], nodes_sol[ms.dim_x+1:dim_node, i]

        for t in ts_seg
            xt, pt, _ = ms.flow(ts_nodes[i], x_node, p_node, t, [tf_sol; λ_sol]; ms.kwargs...)
            push!(t_plot, t)
            push!(x_list, xt)
            push!(p_list, pt)
        end
    end
    return t_plot, hcat(x_list...)', hcat(p_list...)'
end

# shoot!(res, x0, p0, tf, λ)
ms = MultipleShooting(flow, shoot!, x_direct(0), p_direct(0), λ_direct; augment=true);

N_segments = 10
y_guess = generate_initial_guess(ms, x_direct, p_direct, tf_direct, λ_direct, N_segments)

prob = NonlinearProblem(ms, y_guess, N_segments)
sol = solve(prob, TrustRegion(); show_trace=Val(true), maxiters=2)