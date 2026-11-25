# Capture the numeric reference snapshot from the current cache.
#   julia --project=applications/jtcam applications/jtcam/test/capture_reference.jl
# Run once after a successful FORCE_RECOMPUTE=1 paper run.

using Serialization
# Deserialising ODE solutions needs every module that defined their types.
using KEEP: PointMass10, PointMass4, PointMassPara, LimitCycle, Optimization
using OrdinaryDiffEqTsit5, OrdinaryDiffEqVerner
using BifurcationKit

const JTCAM = normpath(joinpath(@__DIR__, ".."))
include(joinpath(JTCAM, "paper_core.jl"))
include(joinpath(JTCAM, "paper_cache.jl"))

function numeric_snapshot(x, prefix="", out=Dict{String,Any}())
    if x isa NamedTuple
        for (k, v) in pairs(x)
            key = isempty(prefix) ? String(k) : "$prefix.$k"
            # Wall-clock timings and derived speedups are not deterministic.
            occursin("btime", key) && continue
            occursin("t_wall", key) && continue
            occursin("speedup", key) && continue
            occursin("ratio_cb", key) && continue
            numeric_snapshot(v, key, out)
        end
    elseif x isa Number
        out[prefix] = Float64(x)
    elseif x isa AbstractArray && eltype(x) <: Number
        out[prefix] = collect(Float64, x)
    end
    return out
end

path = paper_cache_path()
isfile(path) || error("no cache at $path: run paper.jl with FORCE_RECOMPUTE=1 first")
results = open(io -> deserialize(io), path, "r")
snapshot = numeric_snapshot(results)
out = joinpath(@__DIR__, "reference_results.jls")
open(io -> serialize(io, snapshot), out, "w")
println("wrote $out with $(length(snapshot)) numeric entries")
