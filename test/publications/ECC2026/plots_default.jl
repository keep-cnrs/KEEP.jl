import Plots

const WIDTH_PT = 245.71811
const PT_PER_INCH = 72.27  # DPI of Latex
const DPI = 72  # DPI of plotting librairie, cannot be changed in Plots.jl
const FONTSIZE = 8
const LABELFONTSIZE = FONTSIZE
const TICKFONTSIZE = 6
const LEGENDFONTSIZE = TICKFONTSIZE
const ANNOTATIONFONTSIZE = FONTSIZE

const WIDTH_IN = WIDTH_PT / PT_PER_INCH
const FONT = "Computer Modern"

function plot_size(aspect_ratio=4/3; dpi=DPI)
    width_px = round(Int, WIDTH_IN * dpi)
    height_px = round(Int, width_px / aspect_ratio)
    return (width_px, height_px)
end

function my_defaults()
    Plots.default()
    Plots.default(
        dpi = DPI,
        fontfamily = FONT,
        size = plot_size(),
        formatter = :plain,
        markerstrokewidth=0,
        # markersize=5,
        label= "",
        framestyle = :grid,  # semi, zerolines, none
        labelfontsize=LABELFONTSIZE,
        tickfontsize=TICKFONTSIZE,
        legendfontsize=LEGENDFONTSIZE,
        annotationfontsize=ANNOTATIONFONTSIZE, # only work if no position/color/... is specified
        margin=0Plots.mm,
    )
end

my_defaults()

const PALETTE = Plots.palette(:auto)