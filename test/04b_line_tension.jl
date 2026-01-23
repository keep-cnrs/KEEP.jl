using Plots
using StaticArrays
using LinearAlgebra
using StrFormat

using KEEP: TAU0
import KEEP.PointMass4 as PM4
using KEEP.PointMassPara
using KEEP.LimitCycle
using KEEP.TorqueFunction
using KEEP.Visualisation: plot_trajectory_4D

# v_ref = 11.6896, que se passe-t-il ?
# le bras va plus loin
# la génératrice ne résiste pas assez ?

# v_ref = 11.6897: le bras a une très grande amplitude
# v_ref = 11.6896: le bras a une petite amplitude

v_ref = 11.6897
p = build_vbpara(build_para(v_ref=v_ref))

lc = compute_limit_cycle(SA[0, TAU0, 0, 0, 0], p)

t = range(extrema(lc.t)..., length=1001)
tension = PM4.compute_line_tension.(lc.(t), Ref(p))
display(plot(t, tension, label="", xlabel="t (s)", ylabel="Line tension (N)", title="Line tension as a function of time"))

