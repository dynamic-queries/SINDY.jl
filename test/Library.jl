using Test
using BenchmarkTools
include("../src/Library.jl")

@testset "linear" begin
    x = rand(3,3)
    @test x == linear(x)
end

@testset "quadratic" begin
    # Check for p=3
    x = [[2,1,3] [1,2,3] [1,1,1]]
    c1 = x[:,1] .* x[:,1]
    c2 = x[:,1] .* x[:,2]
    c3 = x[:,1] .* x[:,3]
    c4 = x[:,2] .* x[:,2]
    c5 = x[:,2] .* x[:,3]
    c6 = x[:,3] .* x[:,3]
    C = [c1 c2 c3 c4 c5 c6]
    C1 = quadratic(x)
    @test C1 == C

    # Check for p=2
    x = [[1,2] [2,1]]
    c1 = x[:,1] .* x[:,1]
    c2 = x[:,1] .* x[:,2]
    c3 = x[:,2] .* x[:,2]
    C = [c1 c2 c3]
    C2 = quadratic(x)
    @test C2 == C
end

@testset "cubic" begin
    
end

@testset "transcedental" begin
    x = [[0,pi] [2*pi,3*pi]]
    c1 = sin.(x[:,1])
    c2 = sin.(x[:,2])
    c3 = cos.(x[:,1])
    c4 = cos.(x[:,2])
    c5 = tan.(x[:,1])
    c6 = tan.(x[:,2])
    C = [c1 c3 c5 c2 c4 c6]
    C2 = transcedental(x)
    @test size(C2) == (2,6)
    @test size(C) == (2,6)
    @test C == C2
end
