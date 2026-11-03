using Pkg
Pkg.activate(@__DIR__)

using SafeTestsets

regex = r"(\d+)_(\w+)"
# This runner itself, and the sibling SMAIMODE paper script (run it directly).
const SKIP = ("0_ECC2026.jl", "1_SMAIMODE2026.jl")

for file in readdir(@__DIR__)
    file in SKIP && continue
    isnothing(match(regex, file)) && continue
    full_path = joinpath(".", file)
    eval(quote
        @safetestset $full_path begin
            # We use full_path to ensure include works from anywhere
            include($full_path)
        end
    end)
end
