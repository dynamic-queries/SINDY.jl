using ForwardDiff
using Interpolations

@enum DiffType begin
    euler = 1
    polynomial = 2
    TVR = 3 # Total variational regularization
end

# This is the main function that dispatches the other functions based on the required algorithm
#-------------------------------------------------------------------------------------------------
function gradient!(X::Array,V::Array,T::Array,kind::DiffType)
    if kind == euler
        euler!(X,V,T)
    elseif kind == polynomial
        polynomial!(X,V)
    elseif kind == TVR
        TVR!(X,V,T)
    end
end
#-------------------------------------------------------------------------------------------------

# Collection of all algorithms.
#-------------------------------------------------------------------------------------------------

function euler!(X::Array,V::Array,T::Array)
    """
        X : Array{n,n}   : Samples from the time series.
        V : Array{n-1,n} : Computed velocity from the time series.
        T : Array{n-1,n} : This denotes the time step between measurements.
    """
    V = (X[2:end] .- X[1:end-1]) ./ T
end

function polynomial!(X::Array,V::Array)

end

function TVR!(X::Array,V::Array)


end
#-------------------------------------------------------------------------------------------------
