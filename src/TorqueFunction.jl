module TorqueFunction

using Logging

export torque_function, set_torque_function!, reset_torque_function!
export DEFAULT_TORQUE, RATIONAL_TORQUE, LINEAR_TORQUE, CINFTY_TORQUE, CONTINUOUS_TORQUE, DISCONTINUOUS_TORQUE

"""
eq. 68 of `Fechner_et_al_2014_cf_KiteModels`, "Dynamic Model of a Pumping Kite Power System"
Soem more explainations here https://en.wikipedia.org/wiki/Induction_motor#Slip and #Torque

Torque of the generator in response of a rotation speed `dα`.
 - Has maximum of `Cmax` at `Ωmax`
 - is anti-symetric
 - decays to zero as x goes to infinity
 - is Cinfty

But :
 - Decay is in 1/x, while 1/x^2 or more should make more sense (such that power = dα * torque(dα) as dα -> infty = 0)
 - Doesn't have a horizontal tangent at x=0
 - Doesn't use Ωmin — or another coefficinet — to make the function steeper
"""
function torque_function_rational(dα, p)
    a = 2p.Cmax / p.Ωmax
    b = 1 / p.Ωmax^2
    return a * dα / (1 + b * dα^2)
end

"""piece-wise with discontinuous decay"""
function torque_function_discontinuous(dα, p)
    @warn "Using legacy non-smooth torque function `torque_function_legacy` instead of `torque_function`" maxlog = 1
    abs_out = begin
        if p.Ωmin < abs(dα) < p.Ωmax
            p.Cmax * (abs(dα) - p.Ωmin) / (p.Ωmax - p.Ωmin)
        else
            zero(dα)
        end
    end
    return sign(dα) * abs_out
end

"""piece-wise with continuous decay"""
function torque_function_continuous(dα, p)
    @warn "Using legacy non-smooth torque function `torque_function_continuous` instead of `torque_function`" maxlog = 1
    (; Ωmin, Ωmax, Ωlim, Cmax) = p
    abs_out = begin
        if Ωmin < abs(dα) < Ωmax
            (Cmax / (Ωmax - Ωmin)) * (abs(dα) - Ωmin)
        elseif Ωmax < abs(dα) < Ωlim
            Cmax * ((abs(dα) - Ωlim) / (Ωmax - Ωlim))^10
        else
            zero(dα)
        end
    end
    return sign(dα) * abs_out
end

"""hyperbolic sin/cos based version"""
function torque_function_Cinfty(dα, p)
    k = 1  # steepness parameter
    x = dα / p.Ωmax
    return p.Cmax * (k * sinh(k) * x) / (cosh(k * x) - (cosh(k) - k * sinh(k)))
end


"""
A simplified version of `torque_function_rational` where `dα` stays in the linear zone.

Can be used to determine `Ωmax` in post-processing.
"""
function torque_function_linear(dα, p)
    return p.torque_slope * dα
end

@enum TorqueFunctionChoice begin
    RATIONAL_TORQUE
    LINEAR_TORQUE
    CINFTY_TORQUE
    CONTINUOUS_TORQUE
    DISCONTINUOUS_TORQUE
end

const DEFAULT_TORQUE = LINEAR_TORQUE

const CURR_TORQUE_REF = Ref(DEFAULT_TORQUE)

"""
function set_torque_function(choice::TorqueFunctionChoice)

    Change the torque function used in the simulation and return the new function.
"""
function set_torque_function!(choice::TorqueFunctionChoice)
    CURR_TORQUE_REF[] = choice
    @info "Torque function changed to $choice"
end

function get_torque_function()
    return CURR_TORQUE_REF[]
end

function reset_torque_function!()
    set_torque_function!(DEFAULT_TORQUE)
end

function torque_function(dα, p, choice::TorqueFunctionChoice=CURR_TORQUE_REF[])
    if choice == RATIONAL_TORQUE
        return torque_function_rational(dα, p)
    elseif choice == DISCONTINUOUS_TORQUE
        return torque_function_discontinuous(dα, p)
    elseif choice == CONTINUOUS_TORQUE
        return torque_function_continuous(dα, p)
    elseif choice == CINFTY_TORQUE
        return torque_function_Cinfty(dα, p)
    elseif choice == LINEAR_TORQUE
        return torque_function_linear(dα, p)
    end
    error("Unsupported TorqueFunctionChoice enum value: $choice")
end

end  # module