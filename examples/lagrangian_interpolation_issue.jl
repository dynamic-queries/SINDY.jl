using DifferentialEquations
using DataInterpolations

function munge(A::Vector{Vector{Float64}})
    l1 = length(A)
    l2 = length(A[1])
    M = zeros(l2,l1)
    for  i = 1:l1
        M[:,i] = A[i]
    end
    M
end

function lorenz!(du,u,p,t)
    du[1] = p[1]*(u[2]-u[1])
    du[2] = u[1]*(p[2]-u[3]) - u[2]
    du[3] = u[1]*u[2] - p[3]*u[3]
end

u₀ = [1,0,0]
p = [10,12,8/3] # Change d_lorenz
tspan = (0,10)
t = 0:0.01:10
prob = ODEProblem(lorenz!,u₀,tspan,p)
sol = solve(prob,Tsit5(),saveat=t)

U = munge(sol.u)
T = sol.t

interpq = QuadraticInterpolation(U,T)
DataInterpolations.derivative(interpq,1)




DataInterpolations.derivative(interpq,10)



interpl = LagrangeInterpolation(U,T)
DataInterpolations.derivative(interpl,1)




DataInterpolations.derivative(interpl,10)
