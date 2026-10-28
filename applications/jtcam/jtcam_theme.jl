# jtcam_theme.jl
import Plots
using Plots.PlotMeasures

Plots.gr()

# --- Font Setup ---
function setup_font(font_filename="LinLibertine_R.ttf")
    paths = [
        joinpath(@__DIR__, font_filename),
        joinpath(expanduser("~/Library/Fonts"), font_filename),
        joinpath("/usr/share/fonts/truetype/linux-libertine", font_filename)
    ]

    idx = findfirst(isfile, paths)
    if idx !== nothing
        font_id = Plots.GR.loadfont(abspath(paths[idx]))
        if font_id >= 0
            Plots.gr_font_family["linux libertine"] = font_id
            return "linux libertine"
        end
    end

    @warn "$font_filename not found. Drop it in ~/Library/Fonts/ or beside this script. Falling back to 'times'."
    return "times"
end

const FONT = setup_font("LinLibertine_R.ttf")

# --- Dimensions ---
# LaTeX DPI is 72.27 pt/in; JTCAM textwidth is 155 mm (~441 pt), half-width is 75 mm (~213 pt)

# W = \showthe\textwidth = \showthe\columnwidth = \showthe\hsize
# WIDTH_FULL_INCH = W
const PT_PER_INCH = 72.27  # Latex constant
const DPI = 72  # Plots.jl constant

const WIDTH_PT = 415.41025

# --- Typography ---
const FONTSIZE = 12
const LABELFONTSIZE = FONTSIZE
const TITLEFONTSIZE = FONTSIZE
const TICKFONTSIZE = 10
const LEGENDFONTSIZE = 10
const ANNOTATIONFONTSIZE = 10

# --- Sizing Functions & Global Objects ---
"""By default, width is 1 * WIDTH_PT and height is 3/8 * WIDTH_PT (for two 4/3 plots)"""
function plot_size(width=1; height=3/8, dpi=DPI)
    width_in = width * WIDTH_PT / PT_PER_INCH
    height_in = height * WIDTH_PT / PT_PER_INCH
    width_px = round(Int, width_in * dpi)
    height_px = round(Int, height_in * dpi)
    return (width_px, height_px)
end

const SIZE_FULL = plot_size()
const SIZE_HALF = plot_size(0.5)

# --- Palette ---
const PALETTE = Plots.palette(:auto)
# Standard 8-color Okabe-Ito (Wong) palette
# const PALETTE = okabe_ito = [
#     colorant"#E69F00",  # Orange
#     colorant"#56B4E9",  # Sky Blue
#     colorant"#009E73",  # Bluish Green
#     colorant"#F0E442",  # Yellow
#     colorant"#0072B2",  # Blue
#     colorant"#D55E00",  # Vermilion
#     colorant"#CC79A7",  # Reddish Purple
#     colorant"#000000",  # Black
# ]

# const PALETTE = Plots.palette(:Dark2_8)

# JTCAM primary dark blue RGB(10, 40, 80) with harmonious accents
# const PALETTE = [
#     # Plots.RGB(10/255, 40/255, 80/255),  # Main accent color is too dark
#     Plots.RGB(0/255, 140/255, 185/255),  # Opting for the brighter one
#     Plots.RGB(0.85, 0.35, 0.15),
#     Plots.RGB(0.15, 0.60, 0.45),
#     Plots.RGB(0.70, 0.20, 0.50),
#     Plots.RGB(0.40, 0.40, 0.40),
# ]

# const PALETTE = [
#     RGB(0/255, 140/255, 185/255),  # 1. Cerulean Blue (Your anchor)
#     RGB(235/255, 105/255, 65/255),  # 2. Coral Orange
#     RGB(45/255, 170/255, 110/255),  # 3. Emerald Green
#     RGB(165/255, 90/255, 195/255),  # 4. Amethyst Violet
#     RGB(235/255, 175/255, 30/255),  # 5. Marigold Amber
#     RGB(225/255, 70/255, 125/255),  # 6. Raspberry Rose
#     RGB(30/255, 180/255, 175/255),  # 7. Bright Teal
#     RGB(90/255, 110/255, 195/255),  # 8. Cornflower Indigo
# ]


# --- Publication Defaults ---
function my_defaults()
    Plots.default()
    Plots.default(
        thickness_scaling=0.5,
        dpi=DPI,
        fontfamily=FONT,
        size=SIZE_HALF,
        palette=PALETTE,
        formatter=:plain,
        label="",

        # Font sizes
        labelfontsize=LABELFONTSIZE,
        tickfontsize=TICKFONTSIZE,
        legendfontsize=LEGENDFONTSIZE,
        titlefontsize=TITLEFONTSIZE,
        annotationfontsize=ANNOTATIONFONTSIZE,

        # Line & Marker aesthetics
        linewidth=3,
        markersize=7,
        markerstrokewidth=0,

        # Axes & Grid layout
        framestyle=:box,
        grid=true,
        gridalpha=0.1,
        gridlinewidth=1,
        minorgrid=false,
        tick_direction=:none,

        # Margins & Legends
        margin=0Plots.mm, # 4Plots.mm,
        legend=:topright,
        legend_background_color=Plots.RGBA(1.0, 1.0, 1.0, 0.85)
    )
end

const OUTPATH = mkpath(joinpath(@__DIR__, "figs"))
my_defaults()