abstract type DynamicalSystem end

struct LorenzSystem <: DynamicalSystem end
struct LotkaVolterra <: DynamicalSystem end
struct NSECylinder <: DynamicalSystem end
struct HelicopterData <: DynamicalSystem end

function datagen(data::LorenzSystem)
    function noise(u)
        for i = 1:length(u)
            for j = 1:length(u[i])
                u[i][j] += 0.1*randn()
            end
        end
        return u
    end

    function lorenz!(du,u,p,t)
        du[1] = p[1]*(u[2]-u[1])
        du[2] = u[1]*(p[2]-u[3]) - u[2]
        du[3] = u[1]*u[2] - p[3]*u[3]
    end

    u₀ = [1,0,0]
    p = [10,28,8/3] # Change d_lorenz
    tspan = (0,100)
    t = 0.01:0.01:100
    prob = ODEProblem(lorenz!,u₀,tspan,p)
    sol = solve(prob,Tsit5(),saveat=t)
    # noisy = noise(sol.u)
    fig = plot(sol,vars=(1,2,3))
    display(fig)
    return sol
end
