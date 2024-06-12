using DifferentialEquations
using Plots; gr()
using Optimization
using OptimizationNLopt
using PolyChaos

# Variables
tSpan = (0.0,600)

μ = 0.027
k1 = 2.12*(10^(-9))
k2 = 7.13*(10^(-3))
k3 = 168
k33 = 168
k4 = 0.035
D = 0.0396
# D = 0
f = 10^(-3)

T = 5*(10^(6))
Vs = 1.25*(10^(5))
Vd = 0
Id = 0
Is = 0
Ic = 0

w = 0.025*[1,0.4,0.4,0.4,1,1]


# Function
function func(du, u, p, t)
    du[1] = μ*u[1] - k1*(u[5]+u[6])*u[1] + D*(T-u[1])
    du[2] = k1*u[6]*u[1] - (k1*u[5]-μ)*u[2] - D*u[2]
    du[3] = k1*u[5]*u[1] - (k1*u[6]+k2)*u[3] - D*u[3]
    du[4] = k1*(u[5]*u[2]+u[6]*u[3]) - k2*u[4] - D*u[4]
    du[5] = k3*u[3] - (k1*(u[1]+u[2]+u[3]+u[4])+k4+D)*u[5] 
    du[6] = k33*u[4] + f*k3*u[3] - (k1*(u[1]+u[2]+u[3]+u[4])+k4+D)*u[6]
end
function noise(du, u, p, t)
    du[1] = w[1]*u[1]
    du[2] = w[2]*u[2]
    du[3] = w[3]*u[3]
    du[4] = w[4]*u[4]
    du[5] = w[5]*u[5]
    du[6] = w[6]*u[6]
end
problem = SDEProblem(func, noise, [T;Id;Is;Ic;Vs;Vd], tSpan)


# State Estimation
# Bounds
k1_B = [1*10^(-9), 1*10^(-8)]
k2_B = [1*10^(-3), 1*10^(-2)]
k3_B = [100, 200]
k33_B = [100, 200]
k4_B = [0,0.1]
function estimationFunc(du, u, p, t)
    du[1] = μ*u[1] - p[1]*(u[5]+u[6])*u[1] + D*(T-u[1])
    du[2] = p[1]*u[6]*u[1] - (p[1]*u[5]-μ)*u[2] - D*u[2]
    du[3] = p[1]*u[5]*u[1] - (p[1]*u[6]+p[2])*u[3] - D*u[3]
    du[4] = p[1]*(u[5]*u[2]+u[6]*u[3]) - p[2]*u[4] - D*u[4]
    du[5] = p[3]*u[3] - (p[1]*(u[1]+u[2]+u[3]+u[4])+p[5]+D)*u[5] 
    du[6] = p[4]*u[4] + f*p[3]*u[3] - (p[1]*(u[1]+u[2]+u[3]+u[4])+p[5]+D)*u[6]
end
estimation = ODEProblem(estimationFunc, [T;Id;Is;Ic;Vs;Vd], tSpan)
# Dynamics
cbs = []
tStops = []
const estParams = [0.0,0.0,0.0,0.0,0.0,0.0]
function OptimizeParams_NLOPT(interval; maxtime=30) 
    i = 1
    while i*interval <= tSpan[2]
        timeIndex = i*interval
        preTimeIndex = (i-1)*interval
        function c(u,t,intregrator) 
            if (t == timeIndex)
                println("Optimizing at $(timeIndex)...")
            end
            t == timeIndex
        end
        function a!(intregrator)
            ti = time()
            map = Dict()
            function opt(x, p) 
                eSol = solve(estimation, Tsit5(), p=[x[1],x[2],x[3],x[4],x[5]], saveat=preTimeIndex:0.1:timeIndex)
                ev = eSol(timeIndex)
                tv  = intregrator.u
                err = tv-ev
                sq = err[1]^2 + err[2]^2 + err[3]^2 + err[4]^2 + err[5]^2
                map[[x[1],x[2],x[3],x[4],x[5]]] = sq
                return sq
            end
            x0 = [k1_B[1],k2_B[1],k3_B[1],k33_B[1],k4_B[1]]
            p = [1.0]
            funca = OptimizationFunction(opt)
            probl = Optimization.OptimizationProblem(funca, x0, p, lb = [k1_B[1],k2_B[1],k3_B[1],k33_B[1],k4_B[1]], ub = [k1_B[2],k2_B[2],k3_B[2],k33_B[2],k4_B[2]])
            sol = solve(probl, NLopt.G_MLSL(), local_method=NLopt.LN_NELDERMEAD(), maxtime=maxtime)
            if (estParams == [0,0,0,0,0,0]) 
                v = [sol.u..., map[sol.u]]
                estParams[1] = v[1]
                estParams[2] = v[2]
                estParams[3] = v[3]
                estParams[4] = v[4]
                estParams[5] = v[5]
                estParams[6] = v[6]
            else 
                tw = estParams[6] + map[sol.u]
                itw = (tw/estParams[6]) + (tw/(map[sol.u]))
                v = ((tw/estParams[6])/itw)*estParams + ((tw/(map[sol.u]))/itw)*([sol.u...,map[sol.u]])
                estParams[1] = v[1]
                estParams[2] = v[2]
                estParams[3] = v[3]
                estParams[4] = v[4]
                estParams[5] = v[5]
                estParams[6] = v[6]
            end
            println("OptimizeParams_NLOPT: Discrete Event at $(timeIndex). Elapsed time = $(time() - ti).")
        end
        cb = DiscreteCallback(c,a!) 
        push!(cbs,cb)
        push!(tStops, timeIndex)
        i = i + 1
    end
end

function OptimizeParams_PolyChaos(interval; n=5,d=6)
    i = 1
    while i*interval <= tSpan[2]
        timeIndex = i*interval
        preTimeIndex = (i-1)*interval
        function c(u,t,intregrator) 
            if (t == timeIndex)
                println("Optimizing at $(timeIndex)...")
            end
            t == timeIndex
        end
        function a!(intregrator)
            ti = time()
            uniform01 = Uniform01OrthoPoly(d)
            pce_uniform = [convert2affinePCE(k1_B[1],k1_B[2],uniform01),convert2affinePCE(k2_B[1],k2_B[2],uniform01),convert2affinePCE(k3_B[1],k3_B[2],uniform01),convert2affinePCE(k33_B[1],k33_B[2],uniform01),convert2affinePCE(k4_B[1],k4_B[2],uniform01)]
            samples = [evaluatePCE(pce_uniform[1], sampleMeasure(n,uniform01),uniform01),evaluatePCE(pce_uniform[2], sampleMeasure(n,uniform01),uniform01),evaluatePCE(pce_uniform[3], sampleMeasure(n,uniform01),uniform01),evaluatePCE(pce_uniform[4], sampleMeasure(n,uniform01),uniform01),evaluatePCE(pce_uniform[5], sampleMeasure(n,uniform01),uniform01)]
            lErr = Inf
            lP = []
            for k_1 in samples[1]
                for k_2 in samples[2]
                    for k_3 in samples[3]
                        for k_33 in samples[4]
                            for k_4 in samples[5]
                                eSol = solve(estimation, Tsit5(), p=[k_1,k_2,k_3,k_33,k_4], saveat=preTimeIndex:0.1:timeIndex)
                                ev = eSol(timeIndex)
                                tv  = intregrator.u
                                err = (tv[1] - ev[1])^2 + (tv[2] - ev[2])^2 + (tv[3] - ev[3])^2 + (tv[4] - ev[4])^2 + (tv[5] - ev[5])^2
                                if (err < lErr)
                                    lErr = err
                                    lP = [k_1,k_2,k_3,k_33,k_4,err]
                                end
                            end
                        end
                    end
                end
            end
            if (estParams == [0,0,0,0,0,0]) 
                v = lP
                estParams[1] = v[1]
                estParams[2] = v[2]
                estParams[3] = v[3]
                estParams[4] = v[4]
                estParams[5] = v[5]
                estParams[6] = v[6]
            else 
                tw = estParams[6] + lP[6]
                itw = (tw/estParams[6]) + (tw/lP[6])
                v = ((tw/estParams[6])/itw)*estParams + ((tw/lP[6])/itw)*lP
                estParams[1] = v[1]
                estParams[2] = v[2]
                estParams[3] = v[3]
                estParams[4] = v[4]
                estParams[5] = v[5]
                estParams[6] = v[6]
            end
            println("OptimizeParams_PolyChaos: Discrete Event at $(timeIndex). Elapsed time = $(time() - ti).")
        end
        cb = DiscreteCallback(c,a!) 
        push!(cbs,cb)
        push!(tStops, timeIndex)
        i = i + 1
    end
end


OptimizeParams_PolyChaos(50)
sol = solve(problem, SRIW1(), dt=1//(2^(1)), callback=CallbackSet(cbs...), tstops=tStops)

println(estParams)

# Plot
optSol = solve(estimation, Tsit5(), p=[estParams[1],estParams[2],estParams[3],estParams[4],estParams[5]])
graph = plot(dpi=600)
plot!(graph, sol, xaxis="Time (hr)", label=["T" "Id" "Is" "Ic" "Vs" "Vd"])
display(graph)
savefig(string(@__DIR__, "/figs/all_PolyChaos_estimation.png"))

graph0 = plot(dpi=600, xlabel="Time (hr)", ylabel="cells/mL")
plot!(graph0, sol.t, sol[1,:], label="T")
plot!(graph0, optSol.t, optSol[1,:], label="T*", linestyle=:dash)
display(graph0)
savefig(string(@__DIR__, "/figs/T_PolyChaos_estimation.png"))

graph1 = plot(dpi=600, xlabel="Time (hr)", ylabel="cells/mL")
plot!(graph1, sol.t, sol[2,:], label="Id")
plot!(graph1, sol.t, sol[3,:], label="Is")
plot!(graph1, sol.t, sol[4,:], label="Ic")
plot!(graph1, optSol.t, optSol[2,:], label="Id*", linestyle=:dash)
plot!(graph1, optSol.t, optSol[3,:], label="Is*", linestyle=:dash)
plot!(graph1, optSol.t, optSol[4,:], label="Ic*", linestyle=:dash)
display(graph1)
savefig(string(@__DIR__, "/figs/I_PolyChaos_estimation.png"))

graph2 = plot(dpi=600, xlabel="Time (hr)", ylabel="virions/mL")
plot!(graph2, sol.t, sol[5,:], label="Vs")
plot!(graph2, sol.t, sol[6,:], label="Vd")
plot!(graph2, optSol.t, optSol[5,:], label="Vs*", linestyle=:dash)
plot!(graph2, optSol.t, optSol[6,:], label="Vd*", linestyle=:dash)
display(graph2)
savefig(string(@__DIR__, "/figs/V_PolyChaos_estimation.png"))