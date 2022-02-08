function n2_terms(p::Int)
    return Int(p*(p+1)/2)
end

function quadratic(X::Array)
    """
        Assume data is passed in the form x,y,z
        Solve this issue for an array of arbitrary size.
    """
    s = size(X) # 3,1001
    p = s[1]
    s = s[2]
    n = n2_terms(p) # 6
    C = []
    for i=1:p
        for j=i:p
            append!(C,X[i,:].*X[j,:])
        end
    end
    reshape(C,(n,s))
end

function cubic(X::Array)
    """
        Assume data is passed in the form x,y,z
        Solve this issue for an array of arbitrary size.
    """
    s = size(X) # 3,1001
    p = s[1]
    s = s[2]
    n = n2_terms(p) # 6
    C = []
    for i=1:p
        for j=i:p
            for k = j:p
                append!(C,X[i,:].*X[j,:] .* X[k,:])
            end
        end
    end
    # return C
    reshape(C,(10,s))
end




function linear(X::Array)
    """
        Assuming that the data passed is linear.
    """
    return X
end

function basis(X)
    s = size(X)
    ntsteps = s[2]
    nparams = s[1]
    # Since we consider unit, linear and quadratic terms.
    n = 1 + 3 + n2_terms(nparams) + 10 # You are better than this!
    θ = zeros(n,ntsteps)
    θ[1,:] = ones(ntsteps)
    θ[2:4,:] = linear(X)
    θ[5:10,:] = quadratic(X)
    θ[11:end,:] = cubic(X)
    return θ
end
