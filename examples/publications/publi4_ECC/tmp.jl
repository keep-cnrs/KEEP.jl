using StaticArrays
using ComponentArrays: ComponentArray as CA
import KEEP.PointMassPara as PMP
import KEEP.PointMass4 as PM4

p0 = PMP.build_para()
vbp0 = PMP.build_vbpara(CA(p0; torque_slope=1000, r=20, I=1000))
PM4.integrate(SA[1, 2., 2, 2, 0], 10, vbp0)