# L'initialisation de l'état dans l'ancien code (2*10d) supposait que α et dα étaient nuls.
# Pour tester la dynamique pour α et dα non-nuls, on intègre la dynamique puis on fait les tests à un temps quelconque.

using StaticArrays
using BenchmarkTools: @benchmark, @btime

import KEEP.PointMass10 as PM10
import KEEP.PointMass4 as PM4
using KEEP.PointMassPara

function main()
    p = build_para()
    vbp = build_vbpara(p)

    α = 0.
    τ = 1.
    dα = 0.
    dτ = 3.

    u0_4d = SA[α, τ, dα, dτ, 0.]
    u0_10d = PM10.init_u(τ, dτ, p)

    @code_warntype PM10.dynamics!(similar(u0_10d), u0_10d, p, 0.)
    @code_warntype PM4.dynamics(u0_4d, vbp, 0.)

    tf = 10

    PM10.integrate(u0_10d, tf, p)
    PM4.integrate(u0_4d, tf, vbp)

    @btime PM10.integrate($u0_10d, $tf, $p)
    @btime PM4.integrate($u0_4d, $tf, $vbp)

    return nothing
end

main()