abstract type AbstractBasis end

struct PolynomialBasis <: AbstractBasis end
struct TrigBasis <: AbstractBasis end
struct PolyTrigBasis <: AbstractBasis end

function n2_terms(p::Int)
    return Int(p*(p+1)/2)
end

function n3_terms(p::Int)
    return p*p
end

function linear(X::Array,p::Int,n::Int)
    return X
end


function quadratic(X::Array,p::Int,n::Int)
    np = n2_terms(p) # 6
    C = zeros(np,n)  # 6,3
    k = 1
    for i=1:p
        for j=i:p
            C[k,:] = X[i,:] .* X[j,:]
            k += 1
        end
    end
    return C
end

function cubic(X::Array,p::Int,n::Int)
    k = 1
    np = n3_terms(p)
    C = zeros(np,n)
    for i = 1:p
        for j = 1:p
            C[k,:] = X[i,:].^2 .* X[j,:]
            k += 1
        end
    end
    return C
end

function trig(X::Array,p::Int,n::Int)
    C = zeros(3*p,n)
    k = 1
    for i = 1:p
        C[k,:] .= sin.(X[i,:])
        k += 1
    end
    for i = 1:p
        C[k,:] .= cos.(X[i,:])
        k += 1
    end
    for i = 1:p
        C[k,:] .= tan.(X[i,:])
        k += 1
    end
    return C
end

function basis(X::Array,o::PolynomialBasis)
    s = size(X)
    ntsteps = s[2] #3
    p = s[1] #2
    # Since we consider unit, linear and quadratic, cubic terms terms.
    n = 1 + p + n2_terms(p) + n3_terms(p)
    θ = zeros(n,ntsteps)
    θ[1,:] = ones(ntsteps)
    k = 1+p
    l = k+n2_terms(p)
    m = l+n3_terms(p)
    θ[2:k,:] = linear(X,p,ntsteps)
    θ[k+1:l,:] = quadratic(X,p,ntsteps)
    θ[l+1:m,:] = cubic(X,p,ntsteps)
    return θ
end

function basis(X::Array,o::TrigBasis)
    s = size(X)
    p = s[1]
    n = s[2]
    θ = zeros(3*p,n)
    θ .= trig(X,p,n)
    return θ
end

function basis(X::Array,o::PolyTrigBasis)
    o1 = PolynomialBasis()
    o2 = TrigBasis()
    θ1 = basis(X,o1)
    θ2 = basis(X,o2)
    return [θ1;θ2]
end
