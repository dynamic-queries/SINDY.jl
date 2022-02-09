using Test
using BenchmarkTools
using LinearAlgebra
include("../src/Library.jl")

# Quadratic basis
@testset "quadratic" begin
    A = [1 2 3;4 5 6;7 8 10]
    B = quadratic(A,3,3)
    @test cond(A) < 1e6
    @test cond(B) < 1e6
end

# Cubic basis
@testset "cubic" begin
    A = [1 2;3 4]
    B = cubic(A,2,2)
    @test cond(B) <1e6

    A = [1 2 3;4 5 6;7 8 10]
    B = cubic(A,3,3)
    @test cond(B) < 1e6
end

# Test for a small unit matrix
@testset "ones" begin
    A = ones(3,3)
    o = PolynomialBasis()
    B = basis(A,o)
    @test cond(B) == Inf
end

# Test for an identity matrix
@testset "identity" begin
    A = [1 0 0;0 1 0;0 0 1]
    o = PolynomialBasis()
    B = basis(A,o)
    @test cond(B) < 1e6
end

# Condition number of a non square matrix
@testset "non square matrix" begin
    A = [1 0 0 0;0 1 0 0;0 0 0 1]
    o = PolynomialBasis()
    B = basis(A,o)
    @test cond(B) < 1e6
end

# Condition number of random non square matrix
@testset "random non square matrix" begin
    A = 1.0 .+ randn(3,2)
    o = PolynomialBasis()
    B = basis(A,o)
    @test cond(B) < 1e6
end

@testset "basis 2 -  correctness" begin
    A = [1 2;3 4]
    o = PolynomialBasis()
    B = basis(A,o)
    @test B[1,:] == [1,1]
    @test B[2,:] == [1,2]
    @test B[3,:] == [3,4]
end

@testset "basis 3 - correctness" begin
    A = [1 2 3;4 5 6;7 8 10]
    o = PolynomialBasis()
    B = basis(A,o)
    @test size(B,1) == 1 + 3 + 3*2 + 3*3
    @test size(B,2) == 3
    @test B[1,:] == [1,1,1]
    @test B[2,:] == [1,2,3]
    @test B[3,:] == [4,5,6]
    @test B[end,:] == [343,512,1000]
end

@testset "trig" begin
    A = [1 2;3 4]
    B = trig(A,2,2)
    B[1,:] ≈ [sin(1),sin(2)]
    B[2,:] ≈ [sin(3),sin(4)]
    B[3,:] ≈ [cos(1),cos(2)]
    B[4,:] ≈ [cos(3),cos(4)]
end


@testset "trig" begin
    A = [1 2 3; 4 5 6; 7 8 9]
    o = TrigBasis()
    B = basis(A,o)
    @test size(B,1) == 9
    @test size(B,2) == 3
    @test cond(B) < 1e6
end


@testset "PolyTrig basis" begin
    A = [1 2 3; 4 5 6; 7 8 9]
    o = PolyTrigBasis()
    B = basis(A,o)
    @test cond(B) < 1e6
end
