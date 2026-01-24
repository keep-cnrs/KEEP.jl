# Type-stable: Para{T<: Real} l::T...

module PointMassPara

using ComponentArrays

export build_para, build_vbpara, lmt, normalize_vbpara, vbpara_set_params, get_char_time
export VBPara, Para
export NB_LINES

# Type aliases for clarity (still ComponentArrays)
const VBPara = ComponentArray
const Para = ComponentArray

# number of lines
const NB_LINES = 4

"""
Dimensionful parameters.
 - 10d: used for simulation
 - other: used for translation only.
"""
DEFAULT_PARA = ComponentArray(;
    l=2.0, # m
    m=6.0, # kg
    v_ref=9.0, # m s^-1
    g=9.81, # m s^-2
    h_ref=10.0, # m
    n_wind=7, # 1
    ρ_air=1.225, # kg m^-3
    r=31.0, # m
    S=10.0, # m^2
    C_L=0.6, # 1, Should be between 0.7 and 1 "Aerodynamic characterization of a soft kite by in situ flow measurement"
    C_D=0.6 / 7, # 1, Should be 0.2 instead of ~0.085 (L/D between 4 and 6)
    d_l=3e-3, # m
    C_D_l=0.0, # 1, already taken into account in C_D
    ρ_l=707.0, # kg m^-3
    I_eq=1300.0, # kg m^2
    Ωmax=deg2rad(360), # rad s^-1
    Ωmin=deg2rad(30), # rad s^-1
    Ωlim=deg2rad(720), # rad s^-1
    Cmax=15_000.0, # kg m^2 s^-2
    torque_slope=4_000.0, # kg m^2 s^-1 rad^-1
    θ0=deg2rad(65), # rad, 0 for horizontal and 90 at zenith
    φ0=deg2rad(0), # rad
    Δθ=deg2rad(6), # rad
    Δφ=deg2rad(50), # rad
)

function build_para(; kwargs...)::Para
    p = copy(DEFAULT_PARA)
    for (k, v) in kwargs
        p[k] = v
    end
    return p
end

function build_vbpara(p::Para=DEFAULT_PARA)::VBPara
    L, M, T = lmt(p)
    return ComponentArray(;
        l=p.l,
        m=p.m,
        v_ref=p.v_ref,
        g=p.g * T^2 / L,
        h_ref=(p.h_ref / L),
        r=(p.r / L),
        c_L=p.S * p.ρ_air * p.C_L / 2 * L / M,
        f=(p.C_L / p.C_D),
        c_D_l=p.d_l * NB_LINES * p.ρ_air * p.C_D_l / 6 * L^2 / M,
        m_l=π * p.d_l^2 * NB_LINES * p.ρ_l * p.g / 8 * T^2 / M,
        I_eq=(p.I_eq / (L^2 * M)),
        Cmax=p.Cmax * T^2 / (M * L^2),
        Ωmin=(p.Ωmin * T),
        Ωmax=(p.Ωmax * T),
        Ωlim=(p.Ωlim * T),
        torque_slope=p.torque_slope * T / (M * L^2),
        n_wind=p.n_wind,
        θ0=p.θ0,
        φ0=p.φ0,
        Δθ=p.Δθ,
        Δφ=p.Δφ,
    )
end

function build_para(
    p::VBPara; S=DEFAULT_PARA.S, ρ_air=DEFAULT_PARA.ρ_air, d_l=DEFAULT_PARA.d_l
)::Para
    L, M, T = lmt(p)

    S_normalized = S / L^2
    ρ_air_normalized = ρ_air * L^3 / M
    d_l_normalized = d_l / L

    C_L = 2 * p.c_L / (S_normalized * ρ_air_normalized)
    C_D = C_L / p.f
    C_D_l = 6 * p.c_D_l / (d_l_normalized * NB_LINES * ρ_air_normalized)
    ρ_l = 8 * p.m_l / (π * d_l_normalized^2 * NB_LINES * p.g) * M / L^3
    return build_para(;
        l=p.l,
        m=p.m,
        v_ref=p.v_ref,
        g=p.g * L / T^2,
        h_ref=(p.h_ref * L),
        n_wind=p.n_wind,
        ρ_air=ρ_air_normalized * M / L^3,
        r=(p.r * L),
        S=S_normalized * L^2,
        C_L=C_L,
        C_D=C_D,
        d_l=d_l_normalized * L,
        C_D_l=C_D_l,
        ρ_l=ρ_l,
        I_eq=(p.I_eq * M * L^2),
        Cmax=p.Cmax * M * L^2 / T^2,
        Ωmin=(p.Ωmin / T),
        Ωmax=(p.Ωmax / T),
        Ωlim=(p.Ωlim / T),
        torque_slope=p.torque_slope * M * L^2 / T,
        θ0=p.θ0,
        φ0=p.φ0,
        Δθ=p.Δθ,
        Δφ=p.Δφ,
    )
end

"""
Compute the characteristic length, mass and time for the set of parameters `p`.
"""
function lmt(p)
    return (p.l, p.m, p.l / p.v_ref)
end

"""
Initial guess of the period.

Obtained by a gross approximation of the cycle period. During one cycle, the kite goes :
 - once left to right (Δφ) →
 - once right to left (Δφ) ←
 - twice down to up (2Δθ) ↑↑
 - twice up to down (2Δθ) ↓↓

 ┌──▶──╥──◀──┐
 ↑     ⤋     ↑
 └──◀──╨──▶──┘

The length of this approximated path is r * (2Δφ + 4Δθ), and we assume its speed is constant equal to the speed of the wind at the reference height, v_ref.
"""
function get_char_time(p)
    return p.r / p.v_ref * (2p.Δφ + 4p.Δθ)
end

@doc raw"""
Normalize a vbpara (NOT a para), i.e. `l`, `m`, and `v_ref` are set to one.

To normalize a para, use:
```julia
p = build_para()
p_normalized = build_para(normalize_vbpara(build_vbpara(p)))
```
"""
function normalize_vbpara(p::VBPara)::VBPara
    p = copy(p)
    p[[:l, :m, :v_ref]] .= one(eltype(p))
    return p
end

"""
vbp : VBPara -> Para -> set new values -> VBPara

kwargs: kwargs to build_para (parameters that are lost when building VBPara, eg. S, ρ_air, d_l)"""
function vbpara_set_params(
    vbp::VBPara, new_params::ComponentArray=ComponentArray(); kwargs...
)::VBPara
    # p = build_para(vbp; kwargs...)
    # p[keys(new_params)] = new_params

    # It is unsatisfying to not put `kwargs` in the call to build_para. But it fixes the kwargs test.
    p = build_para(vbp)
    p[keys(new_params)] = new_params
    p[keys(kwargs)] = collect(values(kwargs))
    vbp = build_vbpara(p)
    return vbp
end

function vbpara_set_params(vbp, new_params::NamedTuple; kwargs...)
    vbpara_set_params(vbp, ComponentArray(new_params); kwargs...)
end

end  # module
