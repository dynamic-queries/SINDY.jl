module SINDY
    greet() = print("Bear with us while we install SINDY.jl - a library for sparse identification of Dynamical Systems.")
    include("NumDiff.jl")
    include("Library.jl")
    include("Optimizers.jl")
    include("Parsers.jl")
end
