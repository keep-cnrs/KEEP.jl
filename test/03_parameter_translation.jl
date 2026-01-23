using Test
using ComponentArrays
using Logging

using KEEP.PointMassPara

CA = ComponentArray

# Test creation
p = build_para(I_eq = π)
@test p.I_eq ≈ π


# Build a random Para
para = build_para()
para[keys(para)] .= rand(1)
@test :S ∈ keys(para)

# Convert it to VBPara (minimal number of parameters, not normalized)
para_to_vbpara = build_vbpara(para)
@test :S ∉ keys(para_to_vbpara)

# Convert it back to Para, specifying the S, ρ_air, d_l (which is lost information in para_to_vbpara)
para_to_vbpara_to_para = build_para(para_to_vbpara; S=para.S, ρ_air=para.ρ_air, d_l=para.d_l)
@test :S ∈ keys(para_to_vbpara_to_para)

# Check equality between para and para_to_vbpara_to_para
# All numerical fields should equal 0
@test all(para_to_vbpara_to_para .≈ para)

# test lmt
L, M, T = rand(3)
p_lmt = build_para(l=L, m=M, v_ref=L/T)
@test all(lmt(p_lmt) .≈ (L, M, T))
@test all(lmt(build_vbpara(p_lmt)) .≈ (L, M, T))


# Use the `normalize` function on a `VBPara`
normalized_vbpara = normalize_vbpara(para_to_vbpara)
@test all(lmt(normalized_vbpara) .== 1)

# Create a new `VBPara` with a modified value
vbpara_bis = copy(para_to_vbpara)
vbpara_bis.r = 55
@test vbpara_bis.r === 55.

# Use the `vbpara_set_params` function
p = build_para()
vbp = build_vbpara(p)
r, I_eq = rand(2)

p_updated = build_para(r=r, I_eq=I_eq)
vbp_updated_nt = vbpara_set_params(vbp, (; r=r, I_eq=I_eq))
vbp_updated_ca = vbpara_set_params(vbp, CA(r=r, I_eq=I_eq))

@test all(build_para(vbp_updated_nt) .≈ p_updated)
@test all(build_para(vbp_updated_ca) .≈ p_updated)

# `vbpara_set_params` with a specified surface S
S = rand()
p_kw = build_para(S=S)
vbp_kw = build_vbpara(p_kw)

p_kw_updated = build_para(S=S, r=r, I_eq=I_eq)
vbp_updated_nt_kw = vbpara_set_params(vbp, (; r=r, I_eq=I_eq); S=S)
vbp_updated_ca_kw = vbpara_set_params(vbp, CA(r=r, I_eq=I_eq); S=S)

# C_L and C_D
@test all(build_para(vbp_updated_nt_kw, S=S) .≈ p_kw_updated)
@test all(build_para(vbp_updated_ca_kw, S=S) .≈ p_kw_updated)
