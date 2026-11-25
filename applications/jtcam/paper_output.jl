# paper_output.jl — figure serialisation for the JTCAM paper.
#
# Every figure is written as vector SVG and PDF. The SVG viewport is doubled
# (viewBox untouched) so line widths, fonts and markers all scale together —
# no `{width=...}` attributes are needed in the Literate markdown.

# Scale the SVG viewport by an integer factor while leaving its coordinate
# system alone: strokes and text live in viewBox units and therefore scale
# proportionally, never becoming relatively thinner.
function scale_svg(svg::AbstractString; factor::Integer=2)
    m = match(r"width=\"(\d+)\" height=\"(\d+)\"", svg)
    m === nothing && return svg
    w, h = parse(Int, m.captures[1]), parse(Int, m.captures[2])
    replacement = "width=\"$(factor * w)\" height=\"$(factor * h)\""
    return replace(String(svg), m.match => replacement; count=1)
end

# Save both vector formats for a figure. PDF is the canonical rasterisation
# source for the golden-image tests; SVG drives the markdown rendering.
function save_fig(fig, name::AbstractString; factor::Integer=2)
    savefig(fig, joinpath(OUTPATH, "$name.pdf"))
    svg = sprint(io -> show(io, MIME("image/svg+xml"), fig))
    write(joinpath(OUTPATH, "$name.svg"), scale_svg(svg; factor))
    return nothing
end
nothing
