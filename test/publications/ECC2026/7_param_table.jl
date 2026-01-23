using ComponentArrays: ComponentArray as CA
using PyFormattedStrings

import KEEP.PointMassPara as PMP

# | Parameter | Value unit | Description |
# | --- | --- | --- |
# | α0 | 0 rad | Initial angle |

function lt(x, y)
    if (x[1] == '\\') == (y[1] == '\\')
        x < y
    else
        x[1] != '\\'
    end
end

function main()
    p = PMP.build_para()

    # symbol, descr, unit, optional value
    param_list = CA(
        l=("a", "Arm length", raw"\meter", nothing),
        r=("l", "Line length", raw"\meter", nothing),
        m=("m", "Kite mass", raw"\kg", nothing),
        I_eq=("I", "Arm inertia", raw"\kg\square\meter", nothing),
        ρ_air=(raw"\rho_a", "Air density", raw"\kg\per\cubic\meter", nothing),
        d_l=("d_l", "Line diameter", raw"\meter", nothing),
        S=("S", "Kite surface", raw"\square\meter", nothing),
        C_L=("C_L", "Kite lift coefficient", "1", nothing),
        C_D=("C_D", "Kite drag coefficient", "1", "C_D / 7"),
        C_D_l=("C_{D,l}", "Line drag coefficient", "1", nothing),
        ρ_l=(raw"\rho_l", "Line density", raw"\kg\per\cubic\meter", nothing),
        g=("g", "Gravity", raw"\meter\per\square\second", nothing),
        h_ref=("h_r", "Wind reference height", raw"\meter", nothing),
        v_ref=("w_r", "Wind reference speed", raw"\meter\per\second", nothing),
        torque_slope=("b", "Arm braking coefficient", raw"\newton\metre\second\per\radian", nothing),
        θ0=(raw"\theta_0", "Polar angle of the center of the lemniscate", raw"\radian", f"{p.θ0 * 180/π:.0f}" * raw"\pi / 180"),
        Δθ=(raw"\Delta\theta", "Polar angle half-amplitude of the lemniscate", raw"\radian", f"{p.Δθ * 180/π:.0f}" * raw"\pi / 180"),
        Δφ=(raw"\Delta\Phi", "Azimuthal angle half-amplitude of the lemniscate", raw"\radian", f"{p.Δφ * 180/π:.0f}" * raw"\pi / 180"),
    )

    table_str = raw"\begin{table}[thpb]
    \centering
    \caption{Parameter description and values}
    \begin{tabular}{l p{3.2cm} r l}
        \toprule
        {Symbol} & {Description} & {Value} & {Unit} \\
        \midrule
        "

    for key in sort(keys(param_list), by=key -> lowercase(param_list[key][1]), lt=lt)
        sym, descr, unit, optional_value = param_list[key]
        sym_str = raw"$" * sym * raw"$"
        val_str = raw"$" * ifelse(isnothing(optional_value), f"{p[key]:.9g}", optional_value) * raw"$"
        unit_str = raw"\unit{" * unit * "}"
        table_str *= sym_str * " & " * descr * " & " * val_str * " & " * unit_str * raw"\\
        "
    end

    table_str *= raw"\bottomrule
    \end{tabular}
    \end{table}
    "

    println("%% BEGIN TABLE\n\n", table_str, "\n\n%% END TABLE")
end

main()