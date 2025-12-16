function test_aqua()
    @testset "Aqua.jl" begin
        Aqua.test_all(
            KEEP;
            ambiguities=false,
            #stale_deps=(ignore=[:SomePackage],),
            deps_compat=(ignore=[:LinearAlgebra, :Unicode, :Test, :BenchmarkTools],),
            piracies=true,
        )
        # do not warn about ambiguities in dependencies
        Aqua.test_ambiguities(KEEP)
    end
end
