module KEEP


# Value of τ at the start of the cycle, must not be equal to 0 mod π
const TAU0 = 1e-10  # Can be changed to 3π/4
const DEFAULT_TOLERANCE = 1e-9

include("PointMassPara.jl")
include("TorqueFunction.jl")
include("Integrate.jl")
include("PointMass10.jl")
include("PointMass4.jl")
include("LimitCycle.jl")
include("SteadyState.jl")
include("Visualisation.jl")
include("Optimization.jl")

#=
No exports from base module.
Use `using keep.PointMassPara` for exports of submodules
`import keep.PointMass10 as PM10` for aliasing submodule name
and `import keep.LimitCycle: TAU0` for importing a specific function/constant (exported or not by the submodule)
=#

end
