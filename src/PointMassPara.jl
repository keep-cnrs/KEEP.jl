# Type-stable: Para{T<: Real} a::T...

module PointMassPara

using ComponentArrays

export build_para, build_vbpara, lmt, normalize_vbpara, vbpara_set_params, get_char_time
export optimal_torque_slope
export VBPara, Para
export NB_LINES

# Type aliases for clarity (still ComponentArrays)
const VBPara = Any
const Para = Any

# number of lines
const NB_LINES = 4

"""
Dimensionful parameters.
 - 10d: used for simulation
 - other: used for translation only.
"""
DEFAULT_PARA = ComponentArray(;
    a=2.0, # m
    m=6.0, # kg
    w0=9.0, # m s^-1
    g=9.81, # m s^-2
    z0=10.0, # m
    n_wind=7, # 1
    ρ_a=1.225, # kg m^-3
    l=31.0, # m
    S=10.0, # m^2
    C_L=0.6, # 1, Should be between 0.7 and 1 "Aerodynamic characterization of a soft kite by in situ flow measurement"
    C_D=0.6 / 7, # 1, Should be 0.2 instead of ~0.085 (L/D between 4 and 6)
    d_l=3e-3, # m
    C_D_l=0.0, # 1, already taken into account in C_D
    ρ_l=707.0, # kg m^-3
    I=1300.0, # kg m^2
    Ωmax=deg2rad(360), # rad s^-1
    Ωmin=deg2rad(30), # rad s^-1
    Ωlim=deg2rad(720), # rad s^-1
    Cmax=15_000.0, # kg m^2 s^-2
    b=4_000.0, # kg m^2 s^-1 rad^-1
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

function awebox_para()::Para
    return p = build_para(;
        a=2.0,
        l=50.0,
        m=36.8,
        S=3.0,
        I=2000,
        b=1500,
        C_D_l=1.2,
        ρ_l=1464,
        d_l=2e-3,
        z0=10.0,
        w0=9.0,
        n=1 / 0.15,
    )
end

function build_vbpara(p::Para=DEFAULT_PARA)::VBPara
    L, M, T = lmt(p)
    return ComponentArray(;
        a=p.a,
        m=p.m,
        w0=p.w0,
        g=p.g * T^2 / L,
        z0=p.z0 / L,
        l=p.l / L,
        c_L=p.S * p.ρ_a * p.C_L / 2 * L / M,
        f=p.C_L / p.C_D,
        c_D_l=p.d_l * NB_LINES * p.ρ_a * p.C_D_l / 6 * L^2 / M,
        m_l=π * p.d_l^2 * NB_LINES * p.ρ_l * p.g / 8 * T^2 / M,
        I=p.I / (L^2 * M),
        Cmax=p.Cmax * T^2 / (M * L^2),
        Ωmin=p.Ωmin * T,
        Ωmax=p.Ωmax * T,
        Ωlim=p.Ωlim * T,
        b=p.b * T / (M * L^2),
        n_wind=p.n_wind,
        θ0=p.θ0,
        φ0=p.φ0,
        Δθ=p.Δθ,
        Δφ=p.Δφ,
    )
end

function build_para(
    p::VBPara; S=DEFAULT_PARA.S, ρ_a=DEFAULT_PARA.ρ_a, d_l=DEFAULT_PARA.d_l
)::Para
    L, M, T = lmt(p)

    S_normalized = S / L^2
    ρ_a_normalized = ρ_a * L^3 / M
    d_l_normalized = d_l / L

    C_L = 2 * p.c_L / (S_normalized * ρ_a_normalized)
    C_D = C_L / p.f
    C_D_l = 6 * p.c_D_l / (d_l_normalized * NB_LINES * ρ_a_normalized)
    ρ_l = 8 * p.m_l / (π * d_l_normalized^2 * NB_LINES * p.g) * M / L^3
    return build_para(;
        a=p.a,
        m=p.m,
        w0=p.w0,
        g=p.g * L / T^2,
        z0=p.z0 * L,
        n_wind=p.n_wind,
        ρ_a=ρ_a_normalized * M / L^3,
        l=p.l * L,
        S=S_normalized * L^2,
        C_L=C_L,
        C_D=C_D,
        d_l=d_l_normalized * L,
        C_D_l=C_D_l,
        ρ_l=ρ_l,
        I=p.I * M * L^2,
        Cmax=p.Cmax * M * L^2 / T^2,
        Ωmin=p.Ωmin / T,
        Ωmax=p.Ωmax / T,
        Ωlim=p.Ωlim / T,
        b=p.b * M * L^2 / T,
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
    return (p.a, p.m, p.a / p.w0)
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

The length of this approximated path is l * (2Δφ + 4Δθ), and we assume its speed is constant equal to the speed of the wind at the reference height, w0.
"""
function get_char_time(p)
    return p.l / p.w0 * (2p.Δφ + 4p.Δθ)
end

"""
    optimal_torque_slope(p = build_para(); sense = +, tol = 1e-3, maxiter = 5, rtol = 1e-4)

Torque slope `b` that satisfies the kinetic criterion `b = I ω`, with `ω` the
angular frequency of the limit cycle (`ω = 2π/period`, period in seconds from a
simulation). The period depends on `b`, so the value is a fixed point: solve the
cycle, set `b = I 2π/period`, solve again until `b` stops moving. Returns `b`
in the units of `p.b` and prints nothing.

```julia
p0 = build_para()
p0.b = optimal_torque_slope(p0)
vbp = build_vbpara(p0)
```

The limit cycle solver lives in `KEEP.LimitCycle`, which is included after this
module, hence the runtime lookup.
"""
function optimal_torque_slope(p=build_para(); sense=+, tol=1e-3, maxiter=5, rtol=1e-4)
    vbp = build_vbpara(p)
    L, M, T = lmt(vbp)
    b = p.b
    for _ in 1:maxiter
        period = parentmodule(@__MODULE__).LimitCycle.compute_limit_cycle(vbp; sense, tol).t[end]
        b = p.I * 2π / period
        abs(vbp.b * M * L^2 / T - b) <= rtol * b && break
        vbp.b = b * T / (M * L^2)  # physical b -> normalized
    end
    return b
end

@doc raw"""
Normalize a vbpara (NOT a para), i.e. `a`, `m`, and `w0` are set to one.

To normalize a para, use:
```julia
p = build_para()
p_normalized = build_para(normalize_vbpara(build_vbpara(p)))
```
"""
function normalize_vbpara(p::VBPara)::VBPara
    p = copy(p)
    p[[:a, :m, :w0]] .= one(eltype(p))
    return p
end

"""
vbp : VBPara -> Para -> set new values -> VBPara

kwargs: kwargs to build_para (parameters that are lost when building VBPara, eg. S, ρ_a, d_l)"""
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
    return vbpara_set_params(vbp, ComponentArray(new_params); kwargs...)
end

end  # module
