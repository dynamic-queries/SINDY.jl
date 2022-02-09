using DifferentialEquations
using Plots
using DataFrames
using CSV

abstract type DynamicalSystem end

struct ODE1 <: DynamicalSystem end
struct ODE3 <: DynamicalSystem end
struct LorenzSystem <: DynamicalSystem end
struct LotkaVolterra <: DynamicalSystem end
struct NSECylinder <: DynamicalSystem end
struct HelicopterData <: DynamicalSystem
    filename::String
    HelicopterData(file) = new(file)
end

function noise(u::Vector{Vector{T}}) where T
    for i = 1:length(u)
        for j = 1:length(u[i])
            u[i][j] += 0.1*randn()
        end
    end
    return u
end

function datagen(ds::ODE1)
    function ode1!(du,u,p,t)
        du[1] = -0.1*u[1] + 2*u[2]
        du[2] = -2*u[2] - 0.1*u[2]
    end

    u₀ = [1,1]
    p = [0.7,0.3,0.3,0.4]
    tspan = (0,20)
    t = 0.01:0.01:20
    prob = ODEProblem(ode1!,u₀,tspan,p)
    sol = solve(prob,Tsit5(),saveat=t)
    # noisy = noise(sol.u)
    fig = plot(sol,vars=(1,2),title="ODE System Linear")
    display(fig)
    return sol
end

function datagen(ds::ODE3)
    function ode1!(du,u,p,t)
        du[1] = -0.1*u[1]^3 + 2*u[2]^3
        du[2] = -2*u[2]^3 - 0.1*u[2]^3
    end

    u₀ = [1,1]
    p = [0.7,0.3,0.3,0.4]
    tspan = (0,20)
    t = 0.01:0.01:20
    prob = ODEProblem(ode1!,u₀,tspan,p)
    sol = solve(prob,Tsit5(),saveat=t)
    # noisy = noise(sol.u)
    fig = plot(sol,vars=(1,2),title="ODE System Non Linear")
    display(fig)
    return sol
end

function datagen(data::LorenzSystem)
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
    fig = plot(sol,vars=(1,2,3),title="Lorenz Attractor")
    display(fig)
    return sol
end

function datagen(ds::LotkaVolterra)
    function lv!(du,u,p,t)
        α,β,γ,δ = p
        du[1] = α*u[1] - β*u[1]*u[2]
        du[2] = - γ*u[2] + δ*u[1]*u[2]
    end

    u₀ = [1,1]
    p = [0.7,0.3,0.3,0.4]
    tspan = (0,20)
    t = 0.01:0.01:20
    prob = ODEProblem(lv!,u₀,tspan,p)
    sol = solve(prob,Tsit5(),saveat=t)
    # noisy = noise(sol.u)
    fig = plot(sol,vars=(1,2),title="Lotka Volterra")
    display(fig)
    return sol
end

function datagen(ds::HelicopterData)
    df = CSV.read(ds.filename,DataFrame,delim=",")
    data = Matrix(df)
    data = [data[:,3] data[:,4]]
    my_cgrad = cgrad([:blue, :green, :red])
    x = range(0,stop=131,length=size(data,1))
    fig = plot(data[:,1],data[:,2],lc=my_cgrad,line_z=x,title="Helicopter data")
    display(fig)
    return data
end
