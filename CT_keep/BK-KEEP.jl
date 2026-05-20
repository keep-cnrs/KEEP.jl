# cf. https://github.com/bifurcationkit/BifurcationKit.jl/issues/271#issuecomment-4478305759
# https://notes.inria.fr/XnYzM30iQFqRQB3xFyCpEw
#=
 - retrouver la courbe du vent non-optimisée de ECC
 - large intervalle pour I_eq vs ℓ pour voir si un cosinus entre les deux apparait (I_eq ∝ ℓ^2)
=#

# f(t, p, t) 

using KEEP.PointMassPara: build_vbpara
using KEEP.PointMass4: dynamics
using KEEP.LimitCycle: compute_limit_cycle

using BifurcationKit

vbp = build_vbpara()

f