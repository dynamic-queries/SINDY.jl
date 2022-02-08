abstract type AbstractOptimizer end

struct LSTSQ <: AbstractOptimizer
    abstol::Float64
    LSTSQ() = new(1e-16)
end

struct STLSQ <: AbstractOptimizer
    abstol::Float64
    λ::Float64
    STLSQ(thres::Float64) = new(1e-16,thres)
end

##TODO : SR3

function _optimize(θ::Matrix{T},v::Vector{Vector{T}},opt::LSTSQ) where T
    θ = permutedims(θ,(2,1))
    v = permutedims(munge(v))
    return θ\v
end

function _optimize(θ::Matrix{T},v::Matrix{T},opt::LSTSQ) where T
    θ = permutedims(θ,(2,1))
    v = permutedims(v)
    return θ\v
end

function _optimize(θ::Matrix{T},v::Vector{Vector{T}},opt::STLSQ) where T
    iter = 0
    abstol = 1e-8
    θ = permutedims(θ,(2,1))
    v = permutedims(munge(v))
    maxiter = maximum(collect(size(θ)))
    convlimit = abstol # Is the limit on the difference betweem two iterations.
    conv = Inf
    λ = opt.λ
    ξ = θ\v
    smallnums = abs.(ξ) .< λ
    bignums = @. !smallnums[:,1]
    x = similar(ξ)

    for i = 1:maxiter
        y = similar(ξ) # Iteration level least square estimate.
        ξ[smallnums] .= 0
        # θ[:,bignums] * ξ[bignums,i] = v[:,i]
        for i = 1:size(v,2)
            bignums .= @. !smallnums[:,i]
            ξ[bignums,i] .= θ[:,bignums] \ v[:,i]
        end
        smallnums .= abs.(ξ) .< λ
        conv = LinearAlgebra.norm2(y - ξ)
        display(conv)
    end
    ξ[smallnums] .= 0
    return ξ
end

function _optimize(θ::Matrix{T},v::Matrix{T},opt::STLSQ) where T
    iter = 0
    abstol = 1e-8
    θ = permutedims(θ,(2,1))
    v = permutedims(v)
    maxiter = maximum(collect(size(θ)))
    convlimit = abstol # Is the limit on the difference betweem two iterations.
    conv = Inf

    λ = opt.λ
    ξ = θ\v
    smallnums = abs.(ξ) .< λ
    bignums = @. !smallnums[:,1]
    x = similar(ξ)

    for i = 1:maxiter
        smallnums .= abs.(ξ) .< λ
        bignums .= @. !smallnums[:,1]
        y = similar(ξ) # Iteration level least square estimate.
        ξ[smallnums] .= 0
        # θ[:,bignums] * ξ[bignums,i] = v[:,i]
        for i = 1:size(v,2)
            bignums .= @. !smallnums[:,i]
            ξ[bignums,i] .= θ[:,bignums] \ v[:,i]
        end
        conv = LinearAlgebra.norm2(y - ξ)
        display(conv)
    end
    ξ[smallnums] .= 0
    return ξ
end
