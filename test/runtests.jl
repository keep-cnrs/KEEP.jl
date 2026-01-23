using Test
using SafeTestsets
using Underscores
using Logging
import Plots
import Makie
using ProgressMeter

struct PlotSilencer <: AbstractDisplay end

Base.display(d::PlotSilencer, p::Plots.Plot) = nothing

# Still displays Makie figures :(
Base.display(d::PlotSilencer, f::Makie.Figure) = nothing
Base.display(d::PlotSilencer, fap::Makie.FigureAxisPlot) = nothing
Base.display(d::PlotSilencer, s::Makie.Scene) = nothing

function with_no_plots(f)
    silencer = PlotSilencer()
    pushdisplay(silencer)
    try
        return f()
    finally
        popdisplay(silencer)
    end
end



const exclude = [
    "06c_phase_space_plot.jl",  # Uses Makie and is slow
    "05_steady_states.jl",  # SteadyState.all_steady_states throws a mysterious exception, cf. file
    # All slow and old
    "07b_slice_around_optimum.jl",
    "07c_no_constraint_optim.jl",
    "07d_sparse_optim.jl",
    "07e_sparse_optim_variants.jl",
]

const only = [
# "06c_phase_space_plot.jl",
]

function is_included(filename, pattern, exclude, only)
    # @show filename
    # @show only
    # @show isempty(only)
    # @show filename in only
    if !isempty(only)
        return (filename in only)  # if only is not empty, exclude if not in only
    end
    filename in exclude && return false  # else, is filename in excluded list ?
    return !isnothing(match(pattern, filename))  # Else exclude if no match
end

const pattern = r"^\d+[a-z]?_.+\.jl"

function main()
    all_files = @_ readdir(@__DIR__) |> sort
    included_files = @_ all_files |> filter(is_included(_, pattern, exclude, only), __)
    excluded_files = @_ all_files |> filter(!is_included(_, pattern, exclude, only), __)

    println("\tTested files:\n", @_ included_files |> join(__, "\n"), "\n\n")
    println("\tDiscarded files:\n", @_ excluded_files |> join(__, "\n"), "\n\n")

    silencer = PlotSilencer()
    pushdisplay(silencer)
    @time with_logger(SimpleLogger(Warn)) do
        @testset "KEEP" begin
            p = Progress(length(included_files); dt=-1, color=:green)
            next!(p)
            update!(p, 0, force=true)
            for filename in included_files
                p.desc = rpad(filename * " ...", 25)
                next!(p)
                eval(quote
                    @safetestset $filename begin
                        include($filename)
                    end
                end)
            end
        end
    end
    popdisplay(silencer)
    nothing
end

with_no_plots() do
    main()
end

