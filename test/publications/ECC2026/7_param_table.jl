using PyFormattedStrings

import KEEP.PointMassPara as PMP

# | Parameter | Description | Value | Unit

struct Parameter
    key::Symbol
    sym::String
    descr::String
    unit::String
    value_override::Union{String,Nothing}
end


const LATEX_TO_UNICODE = Dict(
    raw"\ell" => "l",
    raw"\rho" => "ρ",
    raw"\theta" => "θ",
    raw"\Delta" => "Δ",
    raw"\phi" => "φ"
)

function eor_transform(sym)
    s = sym
    for (k, v) in LATEX_TO_UNICODE
        s = replace(s, k => v)
    end
    s = replace(s, raw"_" => "", raw"{" => "", raw"}" => "", raw"\\" => "")
    return lowercase(s)
end

function main()
    p = PMP.build_para()

    # symbol, descr, unit, optional value
    params = [
        Parameter(:l, "a", "Arm length", raw"\meter", nothing),
        Parameter(:r, raw"\ell", "Line length", raw"\meter", nothing),
        Parameter(:m, "m", "Kite mass", raw"\kg", nothing),
        Parameter(:I_eq, "I", "Arm inertia", raw"\kg\square\meter", nothing),
        Parameter(:ρ_air, raw"\rho_a", "Air density", raw"\kg\per\cubic\meter", nothing),
        Parameter(:d_l, raw"d_\ell", "Line diameter", raw"\meter", nothing),
        Parameter(:S, "S", "Kite surface", raw"\square\meter", nothing),
        Parameter(:C_L, "C_L", "Kite lift coefficient", "", nothing),
        Parameter(:C_D, "C_D", "Kite drag coefficient", "", "C_L / 7"),
        Parameter(:C_D_l, raw"C_{D,\ell}", "Line drag coefficient", "", nothing),
        Parameter(:ρ_l, raw"\rho_\ell", "Line density", raw"\kg\per\cubic\meter", nothing),
        Parameter(:g, "g", "Gravity", raw"\meter\per\square\second", nothing),
        Parameter(:h_ref, "z_0", "Wind reference altitude", raw"\meter", nothing),
        Parameter(:v_ref, "w_0", "Wind reference speed", raw"\meter\per\second", nothing),
        Parameter(:torque_slope, "b", "Arm braking coefficient", raw"\newton\meter\second\per\radian", nothing),
        Parameter(:θ0, raw"\theta_0", "Lemniscate mean polar angle", raw"\radian", f"{p.θ0 * 180/π:.0f}" * raw"\pi / 180"),
        Parameter(:Δθ, raw"\Delta\theta", "Lemniscate polar half-height", raw"\radian", f"{p.Δθ * 180/π:.0f}" * raw"\pi / 180"),
        Parameter(:Δφ, raw"\Delta\phi", "Lemniscate azimuthal half-width", raw"\radian", f"{p.Δφ * 180/π:.0f}" * raw"\pi / 180"),
    ]


    sort!(params, by=x -> eor_transform(x.sym))
    table_rows = map(params) do p_item
        val = isnothing(p_item.value_override) ? f"{p[p_item.key]:.9g}" : p_item.value_override
        unit_col = isempty(p_item.unit) ? "" : raw"\si{" * p_item.unit * "}"
        f"    ${p_item.sym}$ & {p_item.descr} & ${val}$ & {unit_col} \\\\ "
    end


    table_str = raw"""
\begin{table}[thpb]
\centering
\caption{Parameter Description and Values}
\label{tab:params}
\begin{tabular}{l p{3.2cm} r l}
    \toprule
    {Symbol} & {Description} & {Value} & {Unit} \\
    \midrule
""" * join(table_rows, "\n") * raw""" 
    \bottomrule
\end{tabular}
\end{table}"""

    println("%% BEGIN TABLE 1\n\n", table_str, "\n\n%% END TABLE 1")
end

main()