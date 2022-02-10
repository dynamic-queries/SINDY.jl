include("../src/Library.jl")
using DataFrames
using CSV

function denoise!(sol)

end

function munge(A::Vector{Vector{Float64}})
    l1 = length(A)
    l2 = length(A[1])
    M = zeros(l2,l1)
    for  i = 1:l1
        M[:,i] = A[i]
    end
    M
end

function pprint(x::Matrix,btype::AbstractBasis)
    if size(x,2) == 3 && typeof(btype) == PolyTrigBasis
        labels = ["1","x","y","z","x²","xy","xz","y²","yz","z²",
                 "x³","x²y","x²z","xy²","y³","y²z","xz²","yz²","z³","sin(x)","sin(y)","sin(z)"
                 ,"cos(x)","cos(y)","cos(z)","tan(x)","tan(y)","tan(z)"]
    elseif size(x,2) == 3 && typeof(btype) == PolynomialBasis
        labels = ["1","x","y","z","x²","xy","xz","y²","yz","z²",
                 "x³","x²y","x²z","xy²","y³","y²z","xz²","yz²","z³"]

    elseif size(x,2) == 2 && typeof(btype) == PolyTrigBasis
        labels = ["1","x","y","x²","xy","y²",
                 "x³","x²y","xy²","y³","sin(x)","sin(y)"
                 ,"cos(x)","cos(y)","tan(x)","tan(y)"]
    elseif size(x,2) == 2 && typeof(btype) == PolynomialBasis
        labels = ["1","x","y","x²","xy","y²",
                 "x³","x²y","xy²","y³"]
    end
    df = DataFrame(basis=labels)
    df2 = DataFrame(x,:auto)
    print(hcat(df,df2))
end

struct LocalBasis <: AbstractBasis
    polyorder::Int
    type::AbstractBasis

end

function _remake(ξ,type::AbstractBasis)
    function f(du,u,p,t)

    end
end
