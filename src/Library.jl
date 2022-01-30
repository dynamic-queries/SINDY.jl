function linear(X::Array)
    """
        Assuming that the data passed is linear.
    """
    return X
end

function quadratic(X::Array)
    """
        Assume data is passed in the form x,y,z
        Solve this issue for an array of arbitrary size.
    """
    c1 = X[:,1].^2
    c2 = X[:,2].^2
    c3 = X[:,3].^2
    c4 = X[:,1].*X[:,2]
    c5 = X[:,1].*X[:,3]
    c6 = X[:,2].*X[:,3]
    return [c1 c2 c3 c4 c5 c6]
end

function cubic(X::Array)

end

function transcedental(X::Array)

end

function bessel(X::Array)

end

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
