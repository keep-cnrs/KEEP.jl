using Pkg
Pkg.activate("CT_keep")

using OptimalControl
# import ADNLPModels
import NLPModelsIpopt
using StaticArrays
using ForwardDiff
using LinearAlgebra
using Plots
using ComponentArrays: ComponentArray as CA

# AD Helper functions
grad(f, x) = ForwardDiff.gradient(f, x)
vjp(f, x, dx) = ForwardDiff.derivative(τ -> f(x + τ * dx), 0)

# ================================
# 1. Kinematics
# ================================
function p_c(q, p)
    return SA[q[1], 0]
end

function p_p(q, p)
    return SA[q[1]+p.l*sin(q[2]), p.l*cos(q[2])]
end

# ================================
# 2. Virtual Work (Generalized non-conservative forces)
# ================================
function generalized_forces(q, dq, u, p)
    F = SA[u, 0]
    # vjp computes J(q)*F. Since J is symmetric here, this is equivalent to J'*F.
    return vjp(q_ -> p_c(q_, p), q, F)
end

# ================================
# 3. Lagrangian & Energy
# ================================
function lagrangian(q, dq, u, p)
    # Velocities mapped via Jacobian-Vector Products
    dp_c = vjp(q_ -> p_c(q_, p), q, dq)
    dp_p = vjp(q_ -> p_p(q_, p), q, dq)

    # Kinetic Energy (T) and Potential Energy (V)
    T = 0.5 * p.m_c * sum(abs2, dp_c) + 0.5 * p.m_p * sum(abs2, dp_p)
    V = p.g * p.m_p * p_p(q, p)[2]

    return T - V
end

# ================================
# 4. EOM via Auto-Differentiation (Residual Form)
# ================================
function residuals(q, dq, ddq, u, p)
    # Momentum = ∂L/∂q̇
    momentum(q_, dq_) = grad(v -> lagrangian(q_, v, u, p), dq_)

    # Total time derivative of momentum: d/dt(∂L/∂q̇) 
    A = ForwardDiff.derivative(τ -> momentum(q + τ * dq, dq + τ * ddq), 0)

    # Spatial derivative: ∂L/∂q
    B = grad(q_ -> lagrangian(q_, dq, u, p), q)

    # Generalized non-conservative forces: Q
    Q = generalized_forces(q, dq, u, p)

    return A - B - Q
end

function dyn(x, u, p)
    # Safely extract state components into Static Arrays
    q = x[SA[1, 2]]
    dq = x[SA[3, 4]]

    # Calculate mass matrix (Jacobian of residuals w.r.t ddq) and bias 
    M = ForwardDiff.jacobian(a -> residuals(q, dq, a, u, p), zero(dq))
    b = residuals(q, dq, zero(dq), u, p)

    # Solve for acceleration
    ddq = M \ (-b)

    return vcat(dq, ddq)
end

const p = CA(m_c=5.0, m_p=1.0, l=2.0, g=9.81)

# ================================
# 5. Interface to OptimalControl.jl
# ================================
f(X, u) = dyn(X, u, p)

const tf = 2.0
const X₀ = [0.0, 0.0, 0.0, 0.2]  # Initial condition
const Xf = X₀                    # final condition

ocp = @def begin
    t ∈ [0, tf], time
    X = (x, θ, v, ω) ∈ R⁴, state
    u ∈ R, control

    X(0) == X₀
    X(tf) == Xf

    Ẋ(t) == f(X(t), u(t))

    ∫(u(t)^2) → min
end

init = @init ocp begin
    X(t) := X₀
    u(t) := 0.0
end

# ERROR: CTBase.Exceptions.ExtensionError("missing dependencies to access Ipopt{CPU} options metadata", (:NLPModelsIpopt,), "Ipopt metadata", "Load NLPModelsIpopt extension first: using NLPModelsIpopt")
# Afficher "Load NLPModelsIpopt extension first: using NLPModelsIpopt" plutôt ?
sol = solve(ocp; display=true, grid_size=30, init=init, backend=:manual)

# 2. Explicitly define your strategies
disc = OptimalControl.Collocation(grid_size=100)

# 3. Specify ForwardDiffAD as the backend for the ADNLP modeler
# If the option `backend` is not natively recognized in the metadata as an ADNLPModels type, 
# you can use `bypass()` (or `force()`) to skip strict strategy validation and pass it directly.
mod = OptimalControl.ADNLP(backend=bypass(:default))

# 4. Define the solver
sol = OptimalControl.Ipopt(max_iter=1000, print_level=0)

# 5. Solve the problem
result = solve(ocp; discretizer=disc, modeler=mod, solver=sol)


# ================================
# 6. Extraction & Plotting
# ================================
tsol = time_grid(sol)
Xsol = state(sol).(tsol)
usol = control(sol).(tsol)

X_mat = reduce(hcat, Xsol)
q_sol = X_mat[1:2, :]'
dq_sol = X_mat[3:4, :]'

p1 = plot(tsol, q_sol, label=["x" "θ"], title="Configuration")
p2 = plot(tsol, dq_sol, label=["v" "ω"], title="Velocities")
p3 = plot(tsol, usol, label="u", title="Control", linetype=:steppost)

plot(p1, p2, p3, layout=(3, 1), size=(800, 700))