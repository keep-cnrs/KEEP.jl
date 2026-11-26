# paper_cache.jl — hash-based results cache for the JTCAM paper.
#
# FORCE_RECOMPUTE=1 : always run every numerical study.
# FORCE_RECOMPUTE=0 : reuse a cached result bundle iff the key still matches.
#
# The key is deliberately ONE hash of concatenated hashes: the four compute
# files, the shared core, the whole KEEP src/ tree, the jtcam Project/Manifest
# and VERSION. We do not care *which* input changed, only that something did.
# Plotting and narrative code (paper.jl) is intentionally excluded, so a pure
# figure-style edit never re-runs an ODE or an optimization.

using SHA
using Serialization

const PAPER_COMPUTE_FILES = (
    "paper_core.jl",
    "Section_10D.jl",
    "Section_Modelling.jl",
    "Section_Numerical_analysis.jl",
    "Section_Optimization.jl",
)

# Hash one file's raw bytes into a 32-byte digest.
function _digest_file(path::AbstractString)
    open(path, "r") do io
        return sha256(io)
    end
end

# Hash every Julia source file under KEEP/src, sorted by path. Paths are part
# of the hash so a rename is detected as a change.
function _digest_keep_src()
    root = normpath(joinpath(@__DIR__, "..", "..", "src"))
    files = sort(
        String[
            joinpath(dp, f) for (dp, _, fs) in walkdir(root) for
            f in fs if endswith(f, ".jl")
        ],
    )
    buf = IOBuffer()
    for f in files
        write(buf, relpath(f, root))
        write(buf, _digest_file(f))
    end
    return sha256(take!(buf))
end

function paper_key()
    buf = IOBuffer()
    for f in PAPER_COMPUTE_FILES
        write(buf, _digest_file(joinpath(@__DIR__, f)))
    end
    write(buf, _digest_keep_src())
    write(buf, _digest_file(joinpath(@__DIR__, "Manifest.toml")))
    write(buf, _digest_file(joinpath(@__DIR__, "Project.toml")))
    write(buf, string(VERSION))
    return bytes2hex(sha256(take!(buf)))
end

function paper_cache_path(key=paper_key())
    return joinpath(@__DIR__, "tmp", "paper_cache", "results-$key.jls")
end

# `compute_all()` is defined in paper.jl once all Section files are included.
function load_or_compute()
    key = paper_key()
    path = paper_cache_path(key)
    forced = get(ENV, "FORCE_RECOMPUTE", "0") == "1"
    if !forced && isfile(path)
        open(path, "r") do io
            return deserialize(io)
        end
    end
    results = compute_all()
    mkpath(dirname(path))
    open(path, "w") do io
        return serialize(io, results)
    end
    return results
end
nothing
