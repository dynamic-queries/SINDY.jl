include("../src/SINDY.jl")
using LinearAlgebra
using DataDrivenDiffEq

# Sample and test(integration tests) some of the functionalities for Types
ds = LorenzSystem()
ds1 = ODE1()
ds2 = ODE3()
opt = STLSQ(0.01)
lib = TrigBasis()
lib1 = PolyTrigBasis()

# Integration tests - lorenz attractor
begin
        data = datagen(LorenzSystem())
        savefig("./figures/Traj_Lorenz.svg")
        v = differentiate(data,TotalVariationalDerivativative())
        savefig("./figures/Vel_Lorenz.svg")
        b_type = PolyTrigBasis()
        θ = basis(munge(data.u),b_type)
        LinearAlgebra.cond(θ)
        ξ = _optimize(θ,v,STLSQ(0.1))
        pprint(ξ,b_type)

        # This function needs to be defined manually from the output of pretty print.
        function lorenz_remake!(du,u,p,t)
                du[1] = -9.9756*u[1] + 9.97594*u[2]
                du[2] = 27.7373*u[1] - 0.9388*u[2] - 0.99*u[1]*u[3]
                du[3] = -2.65769*u[3] +  0.99615*u[1]*u[2]
        end

        u₀ = [1,0,0]
        p = [10,28,8/3]
        tspan = (0,100)
        t = 0.01:0.01:100

        prob = ODEProblem(lorenz_remake!,u₀,tspan,p)
        sol = solve(prob,Tsit5(),saveat=t)
        fig = plot(sol,vars=(1,2,3),title="Lorenz Attractor remade")
        display(fig)
        savefig("./figures/Remade_Lorenz.svg")
end


#TODO: The basis function is the only manual part of the implementation. This has to change.
#TODO: Include the forcing function expressions in the defn of the Dynamical System.

# Integration tests - Lotka Volterra

begin
        data = datagen(LotkaVolterra())
        savefig("./figures/Traj_Lotka.svg")
        v = differentiate(data,TotalVariationalDerivativative())
        savefig("./figures/Vel_Lotka.svg")
        b_type = PolynomialBasis()
        θ = basis(munge(data.u),b_type)
        LinearAlgebra.cond(θ)
        ξ = _optimize(θ,v,STLSQ(0.01))
        pprint(ξ,b_type)


        function lk!(du,u,p,t)
                du[1] = 0.7003*u[1] - 0.30014*u[1]*u[2]
                du[2] = -0.299953*u[2] + 0.399979*u[1]*u[2]
        end

        u₀ = [1,1]
        p = []
        tspan = (0,20)
        t = 0.01:0.01:20
        prob = ODEProblem(lk!,u₀,tspan,p)
        sol = solve(prob,Tsit5(),saveat=t)
        fig = plot(sol,vars=(1,2),title="Lotka Volterra remade")
        display(fig)
        savefig("./figures/Remade_Lotka.svg")
end

# Integration - Helicopter data

begin
        data = permutedims(datagen(HelicopterData("./src/data/heli.csv")))
        t = Vector(0.01:0.01:131)
        savefig("./figures/traj_heli.svg")
        plot(t,data[1,:],xlabel="t",ylabel="θ",title="θ")
        savefig("./figures/heli_θ.png")
        plot(t,data[2,:],xlabel="t",ylabel="ψ",title="ψ")
        savefig("./figures/heli_ψ.png")

        # v = differentiate(data,t,TotalVariationalDerivativative())
        # savefig("./figures/velo_heli.svg")
        # v,s = collocate_data(data,t)
        v = zeros(size(data))
        for i = 1:size(data,1)
                int = CubicSpline(data[i,:],t)
                v[i,:] = DataInterpolations.derivative.((int,),t)
        end

        b_type = TrigBasis()
        θ = basis(data,b_type)
        LinearAlgebra.cond(θ)
        ξ = _optimize(θ,v,STLSQ(0.001))
        pprint(ξ,b_type)


        # function helicopter!(du,u,p,t)
        #         du[1] = 0.0546156*sin(u[2]) -0.0431238*cos(u[1]) + 0.022303*cos(u[2]) + 0.0176159*tan(u[1])
        #         du[2] = 0.819696*sin(u[1]) -0.041357*sin(u[2]) + 0.17279*cos(u[1]) +  0.0282611*cos(u[2]) -0.469918*tan(u[1])
        # end

        # function helicopter!(du,u,p,t)
        #         du[1] = 0.0544514*sin(u[2]) -0.0432586*cos(u[1]) + 0.022315*cos(u[2]) + 0.0170897*tan(u[1])
        #         du[2] = 0.78517*sin(u[1]) -0.0393751*sin(u[2]) + 0.170829*cos(u[1]) +  0.0281283*cos(u[2]) -0.439536*tan(u[1])
        # end

        function helicopter!(du,u,p,t)
                du[1] = 0.00168582*sin(u[1]) + 0.0546215*sin(u[2]) -0.0431765 *cos(u[1]) + 0.022313*cos(u[2]) + 0.0160404*tan(u[1])
                du[2] = 0.797168*sin(u[1]) -0.0402192*sin(u[2]) + 0.171233*cos(u[1]) +  0.0281401*cos(u[2]) -0.450491*tan(u[1])
        end

        u₀ = [-0.67,0.67]
        p = []
        tspan = (0,131)
        t = 0.01:0.01:131
        prob = ODEProblem(helicopter!,u₀,tspan,p)
        sol = solve(prob,Tsit5(),saveat=t)

        # Validate the measurement


        temp = munge(sol.u)
        plot(sol.t,temp[1,:],xlabel="t",ylabel="θ",title="θ")
        savefig("./figures/cθ_learned.png")
        plot(sol.t,temp[2,:],xlabel="t",ylabel="ψ",title="ψ")
        savefig("./figures/cψ_learned.png")

        fig = plot(sol,vars=(1,2),title="Helicopter model remade",xlabel="ψ",ylabel="θ")
        display(fig)
        savefig("./figures/cRemade_heli.svg")
end
