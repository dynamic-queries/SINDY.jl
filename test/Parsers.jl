using Test
using BenchmarkTools

include("../src/Parsers.jl")

# Generates data for the Lorenz system.
ds = LorenzSystem()
data = datagen(ds)

# Generates data for a lotka volterra model
ds = LotkaVolterra()
data = datagen(ds)

# Generates data for a first order linear ODE
ds = ODE1()
data = datagen(ds)

# Generates data for a first order non linear ODE
ds = ODE3()
data = datagen(ds)

# Reads the data for the Helicopter Model
ds = HelicopterData("src/data/heli.csv")
data = datagen(ds)
