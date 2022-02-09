using Test
using LinearAlgebra
include("../src/NumDiff.jl")
include("../src/Parsers.jl")

# Lorenz system
data = datagen(LorenzSystem())
v_ana = differentiate(data,LorenzSystem(),AnalyticalDeriv())
v_diff = differentiate(data,FiniteDiff())
v_vtd = differentiate(data,TotalVariationalDerivativative())
error_diff = norm(v_ana[:,1:end-1] .- v_diff[:,:])
error_tvd = norm(v_ana[:,:] .- munge(v_vtd))


# Lotka VOlterra
data = datagen(LotkaVolterra())
v_ana = differentiate(data,LotkaVolterra(),AnalyticalDeriv())
v_diff = differentiate(data,FiniteDiff())
v_vtd = differentiate(data,TotalVariationalDerivativative())
error_diff = norm(v_ana[:,1:end-1] .- v_diff[:,:])
error_tvd = norm(v_ana[:,:] .- munge(v_vtd))

# Helicopter data
data = permutedims(datagen(HelicopterData("src/data/heli.csv")))
t = Vector(0.01:0.01:size(data,2)*0.01)
v_tvd = munge(differentiate(data,t,TotalVariationalDerivativative()))






v_diff = differentiate(data,t,FiniteDiff())
