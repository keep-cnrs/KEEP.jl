# test/runtests.jl

using Test
using ProgressMeter

const TEST_ROOT = @__DIR__
const TEST_DIRS = ["unit", "examples"]

const pattern = r"^\d{2}(\.?\d?)*?_(?<testname>[A-Z][^_]*?)\.jl"

# user-provided controls
const exclude = String[]          # e.g. ["01_BrokenTest.jl"]
const only = String[]          # e.g. ["unit/03_Core.jl"]

"""
    collect_tests()

Return a vector of absolute file paths matching `pattern`,
respecting `exclude` and `only`.
"""
function collect_tests()
    if !isempty(only)
        return [abspath(TEST_ROOT, f) for f in only]
    end

    files = String[]
    for d in TEST_DIRS
        dir = joinpath(TEST_ROOT, d)
        isdir(dir) || continue

        for f in readdir(dir)
            occursin(pattern, f) || continue
            f in exclude && continue
            push!(files, joinpath(dir, f))
        end
    end
    sort(files)
end

"""
    run_tests(files)
"""
function run_tests(files)
    p = Progress(length(files); desc="Running tests")
    for f in files
        @info "Including $(basename(f))"
        include(f)
        next!(p)
    end
end

# ---- entry point ----
@testset "Package tests" begin
    files = collect_tests()
    run_tests(files)
end
