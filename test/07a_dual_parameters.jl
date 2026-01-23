using Test

using ForwardDiff
using ComponentArrays

using KEEP.PointMassPara

## Test that we can take the derivative of an objective function (wich will later result from the integration of the dynamics) with respect to the parameters

function obj(tf, vbp::VBPara)
    (; r, I_eq, v_ref) = vbp
    return tf * I_eq * r^2 * v_ref^2
end

function obj(optim_vars::AbstractArray, vbp::VBPara)
    DType = eltype(optim_vars)
    vbp = DType.(vbp)
    (; tf, optim_params) = optim_vars
    vbp[keys(optim_params)] .= optim_params
    return obj(tf, vbp)
end

optim_vars = ComponentArray(tf=2, optim_params=(r=π, I_eq=exp(1)))
vbp = build_vbpara()

# Tests des fonctions (valeur)
tf = 5
@test obj(tf, vbp) == tf * vbp.I_eq * vbp.r^2 * vbp.v_ref^2
@test obj(optim_vars, vbp) == optim_vars.tf * optim_vars.optim_params.I_eq * optim_vars.optim_params.r^2 * vbp.v_ref^2

# Test de la dérivée
∇obj = ForwardDiff.gradient(x -> obj(x, vbp), optim_vars)
@test ∇obj.tf == optim_vars.optim_params.I_eq * optim_vars.optim_params.r^2 * vbp.v_ref^2
@test ∇obj.optim_params.r == 2 * optim_vars.tf * optim_vars.optim_params.I_eq * optim_vars.optim_params.r * vbp.v_ref^2
@test ∇obj.optim_params.I_eq == optim_vars.tf * optim_vars.optim_params.r^2 * vbp.v_ref^2
