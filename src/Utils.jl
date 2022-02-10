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
    elseif size(x,2) == 2 && typeof(btype) == TrigBasis
        labels = ["sin(x)","sin(y)"
        ,"cos(x)","cos(y)","tan(x)","tan(y)"]
    end
    df = DataFrame(basis=labels)
    df2 = DataFrame(x,:auto)
    display(hcat(df,df2))
end

function _remake(ξ::Matrix{T},type::AbstractBasis) where T
    if size(ξ,2) == 3
        if typeof(type)==PolyTrigBasis
            basis = [(x,y,z)->1,(x,y,z)->x,(x,y,z)->y,(x,y,z)->z,
                    (x,y,z)->x^2,(x,y,z)->x*y,(x,y,z)->x*z,(x,y,z)->y^2,
                    (x,y,z)->y*z,(x,y,z)->z^2,(x,y,z)->x^3,(x,y,z)->x^2*y,
                    (x,y,z)->x^2*z,(x,y,z)->x*y^2,(x,y,z)->y^3,(x,y,z)->y^2*z,
                    (x,y,z)->x*z^2,(x,y,z)->y*z^2,(x,y,z)->z^3,(x,y,z)->sin(x),(x,y,z)->sin(y),
                    (x,y,z)->sin(z),(x,y,z)->cos(x),(x,y,z)->cos(y),(x,y,z)->cos(z),(x,y,z)->tan(x),
                    (x,y,z)->tan(y),(x,y,z)->tan(z)]
        elseif typeof(type) == PolynomialBasis
            basis = [(x,y,z)->1,(x,y,z)->x,(x,y,z)->y,(x,y,z)->z,
                    (x,y,z)->x^2,(x,y,z)->x*y,(x,y,z)->x*z,(x,y,z)->y^2,
                    (x,y,z)->y*z,(x,y,z)->z^2,(x,y,z)->x^3,(x,y,z)->x^2*y,
                    (x,y,z)->x^2*z,(x,y,z)->x*y^2,(x,y,z)->y^3,(x,y,z)->y^2*z,
                    (x,y,z)->x*z^2,(x,y,z)->y*z^2,(x,y,z)->z^3]
        end
        @assert size(ξ,1) == length(basis) "Non Equal basis indices"
        f1 = ξ[1,1]*basis[1]
        f2 = ξ[1,2]*basis[1]
        f3 = ξ[i,3]*basis[1]
        for i = 2:length(basis)
            f1 += basis[i]*ξ[i,1]
            f2 += basis[i]*ξ[i,2]
            f3 += basis[i]*ξ[i,3]
        end
        function f(du,u,p,t)
            du[1] = f1(u[1],u[2],u[3])
            du[2] = f2(u[1],u[2],u[3])
            du[3] = f3(u[1],u[2],u[3])
        end
        return f
    elseif size(ξ,2) == 2
        if typeof(type) == PolyTrigBasis
            basis = [(x,y,z)->1,(x,y,z)->x,(x,y,z)->y,(x,y,z)->x^2,(x,y,z)->x*y,
                    (x,y,z)->y^2,(x,y,z)->x^3,(x,y,z)->x^2*y,(x,y,z)->x*y^2,(x,y,z)->y^3,(x,y,z)->sin(x),
                    (x,y,z)->sin(y),(x,y,z)->cos(x),(x,y,z)->cos(y),(x,y,z)->tan(x),(x,y,z)->tan(y)]
        elseif typeof(type) == PolynomialBasis
            basis = [(x,y,z)->1,(x,y,z)->x,(x,y,z)->y,(x,y,z)->x^2,(x,y,z)->x*y,
                    (x,y,z)->y^2,(x,y,z)->x^3,(x,y,z)->x^2*y,(x,y,z)->x*y^2,(x,y,z)->y^3]
        end
        x_basis(x,y,z) = basis(ξ[:,1],x,y,z)
        y_basis(x,y,z) = basis(ξ[:,2],x,y,z)
        function f(du,u,p,t)
            du[1] = f1(u[1],u[2],u[3])
            du[2] = f2(u[1],u[2],u[3])
        end
        return f
    end
end
