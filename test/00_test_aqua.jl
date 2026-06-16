using Aqua
using Test
using KEEP

function test_aqua()
    @testset "Aqua.jl" begin
        Aqua.test_all(
            KEEP;
            ambiguities=false,
            #stale_deps=(ignore=[:SomePackage],),
            deps_compat=(ignore=[:LinearAlgebra, :Logging],),
            piracies=true,
        )
        # do not warn about ambiguities in dependencies
        Aqua.test_ambiguities(KEEP)
    end
end

# Uncomment to get an error
# test_aqua()