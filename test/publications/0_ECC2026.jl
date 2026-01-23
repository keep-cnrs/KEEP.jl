using SafeTestsets

regex = r"(\d+)_(\w+)"
dir = "ECC2026"

for file in readdir(joinpath(@__DIR__, dir))
    isnothing(match(regex, file)) && continue
    full_path = joinpath(dir, file)
    eval(quote
        @safetestset $full_path begin
            # We use full_path to ensure include works from anywhere
            include($full_path)
        end
    end)
end
