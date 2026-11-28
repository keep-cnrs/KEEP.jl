using Test

using KEEP.PointMassPara
using KEEP.SteadyState

## The symmetry σ: (α, τ) ↦ (-α, τ - π) adds the images of the equilibria, it does not
## replace them: a set given with one member of a σ-pair comes back with both.
# (kept out of 05_steady_states.jl, which runtests.jl excludes)
@testset "process_steady_states adds the σ-images" begin
    q = [0.3, 1.0]
    qσ = [-q[1], q[2] - π]
    same(x, y) = isapprox(cos.(x), cos.(y)) && isapprox(sin.(x), sin.(y))
    pair = SteadyState.process_steady_states([q])
    @test length(pair) == 2
    @test any(x -> same(x, q), pair)
    @test any(x -> same(x, qσ), pair)
end

@testset "eight equilibria, closed under σ" begin
    ss = all_steady_states(build_vbpara())
    @test length(ss) == 8
    @test all(q -> any(x -> isapprox(cos.(x), cos.([-q[1], q[2] - π]); atol=1e-6) &&
                            isapprox(sin.(x), sin.([-q[1], q[2] - π]); atol=1e-6), ss), ss)
end
