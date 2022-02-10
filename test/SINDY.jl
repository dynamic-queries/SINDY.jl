include("../src/SINDY.jl")
using LinearAlgebra

# Sample and test(integration tests) some of the functionalities for Types
ds = LorenzSystem()
ds1 = ODE1()
ds2 = ODE3()
opt = STLSQ(0.01)
lib = TrigBasis()
lib1 = PolyTrigBasis()

# Integration tests - lorenz attractor

data = datagen(LorenzSystem())
savefig("../figures/Traj_Lorenz.svg")
v = differentiate(data,TotalVariationalDerivativative())
savefig("../figures/Vel_Lorenz.svg")
b_type = PolyTrigBasis()
θ = basis(munge(data.u),b_type)
LinearAlgebra.cond(θ)
ξ = _optimize(θ,v,STLSQ(0.1))
_remake(ξ,b_type)
savefig("../figures/Remade_Lorenz.svg")

#TODO: The basis function is the only manual part of the implementation. This has to change.
#TODO: Include the forcing function expressions in the defn of the Dynamical System.
