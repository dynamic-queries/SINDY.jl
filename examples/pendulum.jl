using DataDrivenDiffEq
using LinearAlgebra
using ModelingToolkit
using OrdinaryDiffEq
using Plots
using Random
using Symbolics: scalarize

Random.seed!(1111) # For the noise

# Create a nonlinear pendulum
function pendulum(u, p, t)
    x = u[2]
    y = -9.81sin(u[1]) - 0.3u[2]^3 -3.0*cos(u[1]) - 10.0*exp(-((t-5.0)/5.0)^2)
    return [x;y]
end

u0 = [0.99π; -1.0]
tspan = (0.0, 15.0)
prob = ODEProblem(pendulum, u0, tspan)
sol = solve(prob, Tsit5(), saveat = 0.01)

# Create the data with additional noise
X = sol[:,:] + 0.1 .* randn(size(sol))
DX = similar(sol[:,:])

for (i, xi) in enumerate(eachcol(sol[:,:]))
    DX[:,i] = pendulum(xi, [], sol.t[i])
end

ts = sol.t

prob = ContinuousDataDrivenProblem(X, ts, GaussianKernel() ,
    U = (u,p,t)->[exp(-((t-5.0)/5.0)^2)], p = ones(2))

@variables u[1:2] c[1:1]
@parameters w[1:2]
u = scalarize(u)
c = scalarize(c)
w = scalarize(w)

h = Num[sin.(w[1].*u[1]);cos.(w[2].*u[1]); polynomial_basis(u, 5); c]

basis = Basis(h, u, parameters = w, controls = c)

λs = exp10.(-10:0.1:-1)
opt = STLSQ(λs)
res = solve(prob, basis, opt, progress = false, denoise = false, normalize = false, maxiter = 5000)

parameters(res)
result(res)
