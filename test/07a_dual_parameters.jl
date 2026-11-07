using Test

using ForwardDiff
using ComponentArrays

using KEEP.PointMassPara

## Test that we can take the derivative of an objective function (wich will later result from the integration of the dynamics) with respect to the parameters

function obj(tf, vbp::VBPara)
    (; l, I, w0) = vbp
    return tf * I * l^2 * w0^2
end

function obj(optim_vars::AbstractArray, vbp::VBPara)
    DType = eltype(optim_vars)
    vbp = DType.(vbp)
    (; tf, optim_params) = optim_vars
    vbp[keys(optim_params)] .= optim_params
    return obj(tf, vbp)
end

optim_vars = ComponentArray(tf=2, optim_params=(l=π, I=exp(1)))
vbp = build_vbpara()

# Tests des fonctions (valeur)
tf = 5
@test obj(tf, vbp) == tf * vbp.I * vbp.l^2 * vbp.w0^2
@test obj(optim_vars, vbp) == optim_vars.tf * optim_vars.optim_params.I * optim_vars.optim_params.l^2 * vbp.w0^2

# Test de la dérivée
∇obj = ForwardDiff.gradient(x -> obj(x, vbp), optim_vars)
@test ∇obj.tf == optim_vars.optim_params.I * optim_vars.optim_params.l^2 * vbp.w0^2
@test ∇obj.optim_params.l == 2 * optim_vars.tf * optim_vars.optim_params.I * optim_vars.optim_params.l * vbp.w0^2
@test ∇obj.optim_params.I == optim_vars.tf * optim_vars.optim_params.l^2 * vbp.w0^2
