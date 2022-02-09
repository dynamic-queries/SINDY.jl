include("Parsers.jl")
include("Utils.jl")
using DataInterpolations

abstract type AbstractDiff  end

struct AnalyticalDeriv <: AbstractDiff end
struct FiniteDiff <: AbstractDiff end
struct TotalVariationalDerivativative <: AbstractDiff end

function differentiate(sol,o::TotalVariationalDerivativative)
    # interp = LagrangeInterpolation(munge(sol.u),sol.t)
    interp = QuadraticInterpolation(munge(sol.u),sol.t)
    y_reconstruct = interp.(sol.t)
    v = [DataInterpolations.derivative(interp,k) for k in sol.t]
    w = munge(v)
    if length(sol.u[1]) == 3
        fig = plot(w[1,:],w[2,:],w[3,:])
    elseif length(sol.u[1]) == 2
        fig = plot(w[1,:],w[2,:])
    end
    return v
end

function differentiate(sol,o::FiniteDiff)
    q = sol.u[2:end] .- sol.u[1:end-1]
    q = q./(sol.t[2] - sol.t[1])
    w = munge(q)
    if length(sol.u[1]) == 3
        fig = plot(w[1,:],w[2,:],w[3,:])
    elseif length(sol.u[1]) == 2
        fig = plot(w[1,:],w[2,:])
    end
    display(fig)
    return w
end

function differentiate(sol::Matrix,t,o::TotalVariationalDerivativative)
    interp = QuadraticInterpolation(sol,t)
    y_reconstruct = interp.(t)
    v = [DataInterpolations.derivative(interp,k) for k in t]
    w = munge(v)
    if size(sol,1) == 3
        fig = plot(w[1,:],w[2,:],w[3,:])
    elseif size(sol,1) == 2
        fig = plot(w[1,:],w[2,:])
    end
    display(fig)
    return v
end

function differentiate(sol::Matrix,t,o::FiniteDiff)
    q = sol[:,2:end] .- sol[:,1:end-1]
    q = q./(t[2] - t[1])
    w = q
    if size(sol,1) == 3
        fig = plot(w[1,:],w[2,:],w[3,:])
    elseif size(sol,1) == 2
        fig = plot(w[1,:],w[2,:])
    end
    display(fig)
    return w
end


#--------------- The analytical derivative is computed for testing---------------------#

function differentiate(sol,ds::LorenzSystem,o::AnalyticalDeriv)
    function lorenz(u) # Modify when changing parameters
        u1 = 10*(u[2] - u[1])
        u2 = u[1]*(28 - u[3]) - u[2]
        u3 = u[1]*u[2] - 8/3 * u[3]
        return [u1 u2 u3]
    end
    deriv = []
    for i = 1:length(sol.u)
        append!(deriv,lorenz(sol.u[i]))
    end
    deriv = reshape(deriv,(3,length(sol.u)))
    fig = plot(deriv[1,:],deriv[2,:],deriv[3,:],title="Lorenz - Derivative")
    display(fig)
    return deriv

end

function differentiate(sol,ds::LotkaVolterra,o::AnalyticalDeriv)
    function lk(u) # Modify when changing parameters
        u1 = 0.7*u[1] - 0.3*u[1]*u[2]
        u2 = -0.3*u[2] + 0.4*u[1]*u[2]
        return [u1 u2]
    end
    deriv = []
    for i = 1:length(sol.u)
        append!(deriv,lk(sol.u[i]))
    end
    deriv = reshape(deriv,(2,length(sol.u)))
    fig = plot(deriv[1,:],deriv[2,:],title="Lotka Volterra - Derivative")
    display(fig)
    return deriv
end
