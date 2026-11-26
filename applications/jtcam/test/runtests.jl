# JTCAM paper tests. Run with:
#   julia --project=applications/jtcam applications/jtcam/test/runtests.jl
#
# The golden PNGs in reference_figs/ were rasterised from the legacy PDFs
# (gs, png16m, 150 dpi) before the compute/plot split. Raw PDF bytes are not
# deterministic across processes, but their rasterisation is, so parity is
# checked on PNG bytes.

using Test
using SHA
using Serialization
# Deserialising ODE solutions needs every module that defined their types.
using KEEP: PointMass10, PointMass4, PointMassPara, LimitCycle, Optimization
using OrdinaryDiffEqTsit5, OrdinaryDiffEqVerner
using BifurcationKit

const JTCAM = normpath(joinpath(@__DIR__, ".."))
const FIGS = joinpath(JTCAM, "figs")
const REF = joinpath(@__DIR__, "reference_figs")
include(joinpath(JTCAM, "plot_helpers.jl"))
include(joinpath(JTCAM, "paper_core.jl"))
include(joinpath(JTCAM, "paper_cache.jl"))
include(joinpath(JTCAM, "paper_output.jl"))

function raster_pdf(pdf::AbstractString, out::AbstractString)
    run(
        Cmd([
            "gs",
            "-q",
            "-dNOPAUSE",
            "-dBATCH",
            "-sDEVICE=png16m",
            "-r150",
            "-sOutputFile=$out",
            pdf,
        ]),
    )
    return out
end

function raster_svg(svg::AbstractString, out::AbstractString)
    run(Cmd(["rsvg-convert", "--dpi-x=150", "--dpi-y=150", svg, "-o", out]))
    return out
end

sha256sum(path::AbstractString) = open(io -> bytes2hex(sha256(io)), path, "r")

function numeric_snapshot(x, prefix="", out=Dict{String,Any}())
    if x isa NamedTuple
        for (k, v) in pairs(x)
            key = isempty(prefix) ? String(k) : "$prefix.$k"
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

# A sweep can legitimately contain NaN (no cycle found); treat matching NaNs
# as equal rather than letting isapprox fail on them.
function approx_equal(a, b)
    isapprox(a, b; rtol=1e-6, atol=1e-9) && return true
    if a isa AbstractArray && b isa AbstractArray && size(a) == size(b)
        return all(
            i -> (isnan(a[i]) && isnan(b[i])) || isapprox(a[i], b[i]; rtol=1e-6, atol=1e-9),
            eachindex(a),
        )
    end
    return isnan(a) && isnan(b)
end

function golden_names()
    return [
        splitext(split(line)[2])[1] for
        line in readlines(joinpath(REF, "SHA256SUMS")) if !isempty(strip(line))
    ]
end

@testset "JTCAM paper" begin
    names = golden_names()

    @testset "figure parity ($(length(names)) figures)" begin
        mktempdir() do tmp
            for name in names
                pdf = joinpath(FIGS, "$name.pdf")
                @test isfile(pdf) || continue
                got = raster_pdf(pdf, joinpath(tmp, "$name.png"))
                @test sha256sum(got) == sha256sum(joinpath(REF, "$name.png"))
            end
        end
    end

    @testset "figure completeness" begin
        for name in names
            @test isfile(joinpath(FIGS, "$name.pdf"))
            @test isfile(joinpath(FIGS, "$name.svg"))
        end
    end

    @testset "SVG viewport scaling" begin
        raw = "<svg width=\"600\" height=\"400\" viewBox=\"0 0 2400 1600\">"
        scaled = scale_svg(raw; factor=2)
        @test occursin("width=\"1200\" height=\"800\"", scaled)
        @test occursin("viewBox=\"0 0 2400 1600\"", scaled)
        for name in names
            svg = read(joinpath(FIGS, "$name.svg"), String)
            @test occursin("<svg", svg)
            @test occursin("viewBox=", svg)
        end
    end

    @testset "numeric regression" begin
        refpath = joinpath(@__DIR__, "reference_results.jls")
        if isfile(refpath)
            ref = open(io -> deserialize(io), refpath, "r")
            cache = paper_cache_path()
            @test isfile(cache)
            results = open(io -> deserialize(io), cache, "r")
            got = numeric_snapshot(results)
            for (k, v) in ref
                @test haskey(got, k)
                haskey(got, k) && @test approx_equal(got[k], v)
            end
        else
            @warn "no reference_results.jls yet; run test/capture_reference.jl"
        end
    end
end
