function linear(X::Array)
    """
        Assuming that the data passed is linear.
    """
    return X
end

function n2_terms(p::Int)
    return p*(p+1)/2
end

function quadratic(X::Array)
    """
        Assume data is passed in the form x,y,z
        Solve this issue for an array of arbitrary size.
    """
    s = size(X)
    p = s[2]
    s = s[1]
    n = n_terms(p)
    C = []
    for i=1:p
        for j=i:p
            append!(C,X[:,i].*X[:,j])
        end
    end
    reshape(C,(s,n))
end

function n3_terms(p::Int)

end

function cubic(X::Array)

end

function transcedental(X::Array)
    s = size(X)
    p = s[2]
    s = s[1]
    C = []
    for i = 1:p
        append!(C,sin.(X[:,i]))
        append!(C,cos.(X[:,i]))
        append!(C,tan.(X[:,i]))
    end
    return reshape(C,(s,3*p))
    # We cannot include e^{ix} since this makes the equation inhomogenous.
end

##TODO: Implement later
function bessel(X::Array)

end

##TODO: Implement later
function airy(X::Array)

end

function library(X::Array,F::Array)
    """
    Data : Array X is of the following form.

        --------> p
        |
        |
        |
tsteps

        For ease of access we denote X[:,1] as x , X[:,2] as y ....

    Idea :
        - Construct a library with the kind of basis specified
        - We need linear
        - Quadratic in all variables
        - Cubic
        - Transcedental
        - Bessel
        - Airy
        - Combinations of specific ones
    """

end
