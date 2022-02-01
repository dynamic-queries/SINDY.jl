using ForwardDiff
using DataInterpolations

@enum DiffType begin
    euler = 1
    polynomial = 2 # This is also called Principal Derivative Analysis. See https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2631937/
    TVR = 3 # Total variational regularization
    # Include collaction here as well.
end

# This is the main function that dispatches the other functions based on the required algorithm. This is not the Julianic way of doing things.
#TODO : Implement types for each of the algorithm in order to exploit the multiple dispatch functionality of Julia.
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
        T : Array{1,n} : This denotes the time step between measurements.
    """
    dT = T[2:end] /- T[1:end-1]
    V = (X[2:end] .- X[1:end-1]) ./ dT ###TODO : Modify the tests appropriately.
end


##TODO : Include tests for polynomial approximation.
function polynomial!(X::Array,V::Array,T::Array)
    """
        - We approximate each component of x using a Lagrangian Interpolator
        - We maintain a list of interpolation objects for each significant dimension of the data.
        - Lagrangian Interpolation maybe unstable for some data.
        - The grid is regularly sampled => Runge Effect
        - Alternative is to embedd the regular grid on a Chebyshev grid. ?? How would one even do that!
    """
    s = size(X)
    p = s[2]
    s = s[1]
    Interps = []
    for i = 1:p
        intrep = DataInterpolations.LagrangianInterpolation(X[:,i],T)
        push!(Interps,intrep)
    end

    X_denoised = similar(X)
    for i = 1:p
        X_denoised[:,i] = Intreps[i].(T)
        V[:,i] = Interps[i].derivative.(T)
    end
    return X_denoised,V
end

function TVR!(X::Array,V::Array)

end
#-------------------------------------------------------------------------------------------------
