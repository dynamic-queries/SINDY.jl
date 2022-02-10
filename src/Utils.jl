include("../src/Library.jl")

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

struct LocalBasis <: AbstractBasis
    polyorder::Int

end


function _remake(ξ,type::AbstractBasis)

end
