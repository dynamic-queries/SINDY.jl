module SINDY
    greet() = print("Bear with us while we install SINDY.jl - a library for sparse identification of Dynamical Systems.")

    using DifferentialEquations
    using DataInterpolations
    using Plots
    using ForwardDiff
    include("Parsers.jl")
    include("NumDiff.jl")
    include("Library.jl")
    include("Optimizers.jl")


    # TODO create a SINDY class that takes in data, the optimizer type ...
    # This needs to be automated.
end
