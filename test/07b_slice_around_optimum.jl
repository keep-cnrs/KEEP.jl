using Test
using StrFormat
using Plots
using ForwardDiff: derivative
using Setfield
using Printf
using Logging

import KEEP.PointMass4 as PM4
using KEEP.PointMassPara
using KEEP.LimitCycle
using KEEP.TorqueFunction
using KEEP.Optimization

default(lw=3, formatter=:plain, label="")

tol = 1e-5
set_torque_function!(LINEAR_TORQUE)
p0 = build_para()
syms = [:r, :I_eq, :torque_slope]

lower, upper = make_bounds(p0, syms)
stats, model = optimize(p0, syms, lower, upper, tol=tol);

#=
Optim r, Ieq, torque_slope
check boites + α/dα + tension (positive et < 100 kg/m^2) + sensibilité autour du max
(+ nombre d'équilibres/orbites périodiques + coefficients de floquet, si jamais ça change!)

Comparaison avec les les paramètres calculés par Alain

paramétrer par gradient de vent -> chercher des lignes plus grandes ?, par vitesse du vent -> énergie en v^3 ?
Borne inf/sup normalisée à 0, 1 et courbe des variables par couleur + au dessus production d'énergie
=#

L, M, T = lmt(p0)
para_dims, state_dims, iterate_dims, power_dim = compute_dims(p0, syms)
(; solution, objective) = stats

solution = stats.solution

# Vérifier que les contraintes sont respectées
constraints_neg_log_error = -log10(sum(abs, (model.c!(similar(model.meta.ucon), solution))))
@test constraints_neg_log_error > -log10(tol)

# Vérifier que l'on tourne dans le bon sens : dτ0 > 0
@test solution[3] > 0

# Objectif : avant/après
power_dim = M * L^2 * T^-3
initial_guess = model.meta.x0
initial_obj = model.f(initial_guess) * power_dim
final_obj = model.f(solution) * power_dim

# Où se trouvent la solution dans la boite des contraintes ?
solution_dict = Dict(
    :Symbol => syms,
    :InitialValue => initial_guess[end-length(syms)+1:end] .* para_dims,
    :Lower => lower,
    :Solution => solution[end-length(syms)+1:end] .* para_dims,
    :Upper => upper
)
df = solution_dict # Placeholder name to keep using df in the formatted string if needed, or update the formatted string to just print the dict nicely.
# Actually, let's just create a nice string representation

function print_table(syms, vals...)
    # simple table printing
    header = ["Symbol", "Initial", "Lower", "Solution", "Upper"]
    data = hcat(string.(syms), [map(x -> @sprintf("%.4g", x), v) for v in vals]...)
    # ... simple logic to format this
    return data
end
# Simplify: Just print the arrays
solution_str = "
Symbol      : $(syms)
InitialValue: $(initial_guess[end-length(syms)+1:end] .* para_dims)
Lower       : $(lower)
Solution    : $(solution[end-length(syms)+1:end] .* para_dims)
Upper       : $(upper)"

@info f"
[Solution validity]
Status: \%s(stats.status)
-log(norm(c)) > -log(tol) : \%.2f(constraints_neg_log_error) > \%.2f(-log10(tol)) : \%s(constraints_neg_log_error > -log10(tol))
dτ0 < 0 : \%.2f(solution[3]) < 0 : \%s(solution[3] < 0)

[Objective improvement]
Initial objective : \%.2f(initial_obj) W
Final objective   : \%.2f(final_obj) W \%+.2f(100 * (final_obj - initial_obj) / initial_obj)%
Initial objective : %.2f(initial_obj) W
Final objective   : %.2f(final_obj) W %+.2f(100 * (final_obj - initial_obj) / initial_obj)%

[Solution]
%s(df)"

#Contraintes d'état -π/2 < α(t) < π/2, 0 < T(t) < 100S

# Variation du gradient de vent
# include("utils.jl");

function variation(var_sym, var_vals, p0, syms, lower, upper; kwargs...)
    objs = fill(NaN, length(var_vals))
    params = fill(NaN, length(var_vals), length(syms))
    prev_optimal_params = p0[syms]

    @time with_logger(NullLogger()) do
        for (i, val) in enumerate(var_vals)
            p0_loc = @set p0[var_sym] = val
            # p0_loc[syms] = prev_optimal_params
            para_dims, _, _, power_dim = compute_dims(p0_loc, syms)
            stats, model = optimize(p0_loc, syms, lower, upper, tol=tol)
            solution = stats.solution * ifelse(stats.status == :first_order, 1, NaN)
            objs[i] = model.f(solution) * power_dim
            prev_optimal_params .= solution[end-length(syms)+1:end] .* para_dims
            params[i, :] = prev_optimal_params
        end
    end

    plot(ylabel="Objective (W)")
    plot!(var_vals, objs, c=:red)
    vline!([p0[var_sym]], c=:black, lw=1)
    # hline!([0], alpha=0)
    fig_obj = plot!()
    normalized_params = (params .- lower') ./ (upper' - lower')
    plot(xlabel=string(var_sym), ylabel="Normalized parameters", yticks=([0, 1], ["LB", "UB"]), ylims=(0, 1))
    plot!(var_vals, normalized_params, label=permutedims(string.(syms)))
    fig_params = vline!([p0[var_sym]], c=:black, lw=1)

    figs = reshape([fig_obj, fig_params], 2, 1)
    return objs, params, plot(figs..., layout=size(figs))
end

# 10 seconds
_, _, fig = variation(:n_wind, 3:0.2:10, p0, syms, lower, upper, tol=1e-3)
display(fig)

# 10 seconds
_, _, fig = variation(:r, 40:2:80, p0, [:I_eq], [1000], [10000], tol=1e-3)
display(fig)

# Replace @test in optimize by if ... @warning ... end
# 10 seconds
_, _, fig = variation(:I_eq, 1000:1000:10000, p0, [:r], [10], [100], tol=1e-5)
display(fig)

function average_power(p; tol)
    lc = compute_limit_cycle(build_vbpara(p); sense=-, save_everystep=true, tol=tol)
    return lc.u[end][end] / lc.t[end]
end

function landscape_slices(var_syms, N_vals, p0, syms, lower, upper; kwargs...)
    # take one symbol
    # for each value in lb ub, find the  limit cycle
    # store the energy in a n_symb * N_vals matrix
    objs = fill(NaN, N_vals, length(var_syms))
    @time with_logger(NullLogger()) do
        for (i, var_sym) in enumerate(var_syms)
            for (j, val) in enumerate(range(lower[i], upper[i], length=N_vals))
                p0_loc = @set p0[var_sym] = val
                L, M, T = lmt(p0_loc)
                power_dim = M * L^2 * T^-3
                objs[j, i] = average_power(p0_loc; tol=tol) * power_dim
            end
        end
    end
    return objs
end

objs = landscape_slices(syms, 100, p0, syms, lower, upper, tol=tol)
plot(xlabel="Normalized parameters", ylabel="Objective (W)")
plot!(objs', label=permutedims(string.(syms)), c=palette(:tab10))
vline!([(lower[i] + upper[i]) / 2 for i in 1:length(syms)], c=:black, lw=1)
display(plot!())
