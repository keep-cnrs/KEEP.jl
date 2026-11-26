# Regenerates the figures (figs/*.pdf and *.png) and tables (tables/*.tex) of the paper.
#
#   julia --project=applications/jtcam-tex applications/jtcam-tex/main.jl [sections...] [--fast]
#
# sections: modelling, 10d, analysis, optimization (default: all of them).
# --fast: reduced iteration limits and grids, to test the workflow quickly
#         (the figures and tables are then NOT the ones of the paper).
ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")  # headless GR

const SECTIONS = ["modelling", "10d", "analysis", "optimization"]

function main(args=ARGS)
    fast = "--fast" in args
    wanted = filter(a -> !startswith(a, "--"), args)
    isempty(wanted) && (wanted = SECTIONS)
    unknown = setdiff(wanted, SECTIONS)
    isempty(unknown) || error("unknown section(s): $(join(unknown, ", ")); expected $(join(SECTIONS, ", "))")
    # section name => (source file, entry function)
    entry = Dict("modelling" => ("modelling.jl", :run_modelling),
        "10d" => ("numerics_10d.jl", :run_10d),
        "analysis" => ("analysis.jl", :run_analysis),
        "optimization" => ("optimization.jl", :run_optimization))
    for name in SECTIONS
        name in wanted || continue
        @info "=== section: $name (fast=$fast)"
        file, fname = entry[name]
        t = @elapsed begin
            include(joinpath(@__DIR__, "src", file))
            run = Base.invokelatest(getglobal, Main, fname)  # defined after main() started
            Base.invokelatest(run; fast)
        end
        @info "=== done: $name in $(round(t; digits=1)) s"
    end
    return nothing
end

main()
