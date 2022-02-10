using Test
using BenchmarkTools
include("../src/Optimizers.jl")

@testset "Least squares - simple" begin
    A = 1.0 .+ rand(3,3)
    b = rand(3)
    _optimize(A,b,LSTSQ())
end

@testset "Least squares - simple HD" begin
    A = 1.0 .+ rand(3,3)
    b = rand(2,3)
    _optimize(A,b,LSTSQ())
end

@testset "Thresholded Least squares - matrix" begin
    A = 1.0 .+ rand(3,3)
    b = rand(2,3)
    _optimize(A,b,STLSQ(1.0))
end
