module Integrate

using OrdinaryDiffEqTsit5
import LinearAlgebra: norm

"""
Integrate the dynamics from u0 at time t=0 for a given time horizon tf.
p is either a para or a vbpara"""
function _integrate(dynamics, u0, tf, p, alg=Tsit5(); save_everystep, tol, kwargs...)
    ode_prob = ODEProblem(dynamics, u0, tf, p)
    return solve(ode_prob, alg; save_everystep=save_everystep, abstol=tol, reltol=tol, internalnorm=(u, t) -> norm(u), kwargs...)
end


end