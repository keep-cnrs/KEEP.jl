α = (t, sol) -> sol(t, idxs=1)
τ = (t, sol) -> sol(t, idxs=2)
a_vec = (t, sol) -> begin
    α_ = α(t, sol)
    a = sol.prob.p.l
    a * SA[cos(α_), sin(α_), 0]
end
θ = (t, sol) -> begin
    τ_ = τ(t, sol)
    (; θ0, Δθ) = sol.prob.p
    θ = θ0 + Δθ * sin(2τ_)
end
φ = (t, sol) -> begin
    τ_ = τ(t, sol)
    (; φ0, Δφ) = sol.prob.p
    φ_ = φ0 + Δφ * sin(τ_)
end
r_hat = (t, sol) -> begin
    θ_ = θ(t, sol)
    φ_ = φ(t, sol)
    SA[sin(θ_)*cos(φ_), sin(θ_)*sin(φ_), cos(θ_)]
end
r = (t, sol) -> begin
    r_hat_ = r_hat(t, sol)
    a_hat_ = a_vec(t, sol)
    dotprod = r_hat_' * a_hat_
    l = sol.prob.p.r
    a = sol.prob.p.l
    dotprod + sqrt(dotprod^2 + l^2 - a^2)
end
r_vec = (t, sol) -> begin
    r_hat_ = r_hat(t, sol)
    r_ = r(t, sol)
    r_hat_ * r_
end
dα = (t, sol) -> sol(t, idxs=3)
dτ = (t, sol) -> sol(t, idxs=4)
ddα = (t, sol) -> sol(t, Val{1}, idxs=3)
T = (t, sol) -> begin
    dα_ = dα(t, sol)
    TF.torque_function(dα_, sol.prob.p) * dα_
end
F_l_vec = (t, sol) -> begin
    r_hat_ = r_hat(t, sol)
    a_vec_ = a_vec(t, sol)
    l_tension = r_hat_ ⋅ (SA[0, 0, 1] × a_vec_)
    A_tension = -p.I_eq / l_tension * r_hat_ * SA[1, 0]'
    b_tension = -T(t, sol) / l_tension * r_hat_
    F_tension = A_tension[:, 1] * ddα(t, sol) + b_tension
end
F_l = (t, sol) -> norm(F_l_vec(t, sol))
P = (t, sol) -> sol(t, Val{1}, idxs=5)
