abstract type AbstractDiff  end

struct AnalyticalDeriv <: AbstractDiff end
struct FiniteDiff <: AbstractDiff end
struct TotalVariationalDerivativative <: AbstractDiff end

function differentiate(sol,o::TotalVariationalDerivativative)
    # interp = LagrangeInterpolation(munge(sol.u),sol.t)
    interp = QuadraticInterpolation(munge(sol.u),sol.t)
    y_reconstruct = interp.(sol.t)
    v = [DataInterpolations.derivative(interp,k) for k in sol.t]
    temp = munge(v)
    display(plot(temp[1,:],temp[2,:],temp[3,:]))
    return v
end


function differentiate(sol,o::FiniteDiff)
    q = sol.u[2:end] .- sol.u[1:end-1]
    q = q./(sol.t[2] - sol.t[1])
    w = munge(q)
    fig = plot(w[1,:],w[2,:],w[3,:])
    display(fig)
    return w
end

function differentiate(sol,ds::LorenzSystem,o::AnalyticalDeriv)
    function lorenz(u) # Modify when changing parameters
        u1 = 10*(u[2] - u[1])
        u2 = u[1]*(12 - u[3]) - u[2]
        u3 = u[1]*u[2] - 8/3 * u[3]
        return [u1 u2 u3]
    end
    deriv = []
    for i = 1:length(sol.u)
        append!(deriv,lorenz(sol.u[i]))
    end
    deriv = reshape(deriv,(3,length(sol.u)))
    fig = plot(deriv[1,:],deriv[2,:],deriv[3,:])
    display(fig)
    return deriv
end
