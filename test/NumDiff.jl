using Test
using LinearAlgebra
include("../src/NumDiff.jl")

function sample(f,t1::Vector,t2::Vector)
    Y = t2 * ones(length(t2))'
    X = t1 * ones(length(t2))'
    S = f(X,Y)
    S
end

@testset "Euler" begin
    f(x) = sin.(x)
    df(x) = cos.(x)
    t = Vector(-pi:0.01:pi)
    x = f(t)
    step = 0.01*ones(length(x)-1)
    v = zeros(length(step))
    v = gradient!(x,v,step,euler)
    display(norm(v - df(t)[1:end-1]))
    @test norm(v - df(t)[1:end-1]) < 1e-1
    using Plots
    plot(t,df(t))
    scatter!(t[1:end-1],v)
end

@testset "Polynomial" begin
    # Generate data
    f(x,y) = sin.(x*y)
    t1 = -pi:0.01:pi
    t2 = -pi:0.01:pi
    S = sample(f,Vector(t1),Vector(t1))

    # Define an interpolant : Linear Interpolation is bad. One needs a better way to do this.
    it = LinearInterpolation((t1,t2),S)

    # Check the accuray of interpolation
    xtest = -pi:1:pi
    ytest = -pi:1:pi
    R = sample(f,Vector(xtest),Vector(ytest))
    @test norm(S - it(t1,t2)) < 1e-10
    @test norm(R - it(xtest,ytest)) < 1
end
