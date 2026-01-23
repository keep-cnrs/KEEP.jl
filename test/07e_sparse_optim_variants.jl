#=
X :
    - 1:5N+5 : u
    - tf
    - p

c = [
    α0 - αf
    τ0 - TAU0
    dα0 - dαf
    dτ0 - dτf
    W0 - 0
    u(i+1) - u(i) - tf/N * 1/2 * (f(u(i)) + f(u(i+1)))
    τf - τ0 - 2π
]
=#

using SplitApplyCombine
using Setfield
using Test
using ADNLPModels, NLPModelsIpopt
using NonlinearSolve
using StrFormat

using KEEP: TAU0
import KEEP.PointMass4 as PM4
using KEEP.PointMassPara
using KEEP.LimitCycle

#=
input : all normalized
output : all normalized
=#
function optimize(vbp, N, syms, lb_vbp, ub_vbp)
    global X0, lb, ub, U_end, sol1
    U_end = 5N + 5
    get_u(X, i) = @view X[5i-4:5i]
    get_tf(X) = X[U_end+1]
    get_p(X) = @view X[U_end+2:end]

    function objective(X, syms)
        N = get_N(X, syms)
        return get_u(X, N + 1)[end] / get_tf(X, N)
    end

    """cons is of size U_end + 1 = 5N + 6"""
    function cons!(cons, X, syms)
        U = [get_u(X, i) for i in 1:N+1]
        C = [get_u(cons, i) for i in 1:N+1]
        Δt = get_tf(X) / N

        # Precompute all dynamics in cache C
        # then overwrite with trapezoid residuals
        vbp_loc = @set vbp[syms] = get_p(X)
        for i in 1:N+1
            C[i] .= PM4.dynamics(U[i], vbp_loc, 0)
        end
        for i in N+1:-1:2
            @. C[i] = U[i] - U[i-1] - Δt / 2 * (C[i] + C[i-1])
        end
        α0, τ0, dα0, dτ0, W0 = U[1]
        αf, τf, dαf, dτf, _ = U[end]
        C[1][1] = α0 - αf
        C[1][2] = τ0 - TAU0
        C[1][3] = dα0 - dαf
        C[1][4] = dτ0 - dτf
        C[1][5] = W0 - 0
        cons[U_end+1] = τf - τ0 + 2π
        return cons
    end

    @info "Finding limit cycle for initial guess"
    lc = first([lc for lc in all_limit_cycles(vbp; save_everystep=true) if first(lc.u)[4] < 0])

    # [u(t) for t in range(0, T, length=100); T; params]
    X0 = [
        flatten([lc(t) for t in range(0, lc.t[end], length=N + 1)])...,
        lc.t[end],
        vbp[syms]...
    ]

    @test get_u(X0, 1) == lc(0)
    @test get_u(X0, N + 1) == lc(lc.t[end])
    @test get_tf(X0) == lc.t[end]
    @test get_p(X0) == collect(vbp[syms])

    # function make_init_nl_func(U_end, X0)
    #     """want to find its root for 0 residual initial guess"""
    #     function init_nl_func(res, X_partial, p)
    #         X = @set X0[1:U_end+1] = X_partial
    #         cons!(res, X, syms)
    #         return res
    #     end
    # end

    # @info "Altering the initial guess such that it verifies the constraints"
    # X_partial = X0[1:U_end+1]
    # nlp_init = NonlinearProblem(make_init_nl_func(U_end, X0), X_partial, ())
    # Δt1 = @elapsed sol1 = solve(nlp_init, NewtonRaphson())
    # X1 = @set X0[1:U_end+1] = sol1.u
    # @show sum(abs, sol1.prob.u0 - sol1.u)
    # @info f"Initial guess error
    # Before adjustment: \%.2e(sum(abs, cons!(ones(5N+6), X0, syms)))
    # After adjustment:  \%.2e(sum(abs, cons!(ones(5N+6), X1, syms)))
    # Difference: \%.2e(sum(abs, X0 - X1))
    # \%.0f(sol1.stats.nsteps) steps in \%.2f(Δt1)
    # "


    αmin, αmax = -π, π
    τmin, τmax = TAU0 - 2π, TAU0
    dαmax = 10maximum([get_u(X0, i)[3] for i in 1:N+1])
    dαmin = -dαmax
    dτmin, dτmax = minmax(0, 10argmax(abs, [get_u(X0, i)[4] for i in 1:N+1]))
    Wmin, Wmax = 0, 10get_u(X0, N + 1)[5]

    ε = 1e-1
    lb_state = [αmin, τmin, dαmin, dτmin, Wmin] .- ε
    lb = [repeat(lb_state, outer=N + 1); 0; lb_vbp[syms]]

    ub_state = [αmax, τmax, dαmax, dτmax, Wmax] .+ ε
    ub = [repeat(ub_state, outer=N + 1); 10lc.t[end]; ub_vbp[syms]]
    @test all(lb .<= X0 .<= ub)
    @test all(lb_state .<= ub_state)


    # sum(abs, cons!(ones(5N+6), X, syms)[5:5:N+5]) / N  # Doit être en N^2, ordre des trapèzes

    model = ADNLPModel!(X -> X[5N+6] / X[5N+7],
        X0, lb, ub,
        (c, X) -> cons!(c, X, syms), zeros(5N + 6), zeros(5N + 6); minimize=false, backend=:generic)

    @info "Starting solve"
    stats = ipopt(model; tol=1e-3, max_wall_time=60.)
    return stats
end

# N   10   20   30  40
# t   1.8  6.4  16  55
# ram 0.5  1.5  3.5 12
N = 30
p = build_para()
vbp = build_vbpara(p)
syms = [:r, :I_eq]
mult = 3
lb_vbp = @set vbp[syms] = vbp[syms] / mult
ub_vbp = @set vbp[syms] = vbp[syms] * mult
stats = optimize(vbp, N, [:r, :I_eq], lb_vbp, ub_vbp)

build_para(@set vbp[syms] = stats.solution[end-length(syms)+1:end])[syms]

X_part = X0[1:U_end+1]
res_func = make_res(U_end, X0)
res_func(similar(X_part), X_part)

nlp_init = NonlinearProblem(res_func, X_part)
solve(nlp_init)

n = 10
make_func(n) = (A = rand(n, n); b = rand(n); (res, u, p) -> (res .= A * u + b))
nlp = NonlinearProblem(make_func(n), ones(n), ())
t = @elapsed sol = solve(nlp; store_trace=Val(true))
sum(abs, sol.resid)

sol