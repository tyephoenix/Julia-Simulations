using DifferentialEquations
using Distributions
using StatsPlots; gr()
using Optimization
using OptimizationNLopt
using PolyChaos
using Random

Random.seed!(25)


# Vars
tSpan = (0.0,20.0)
trajectories = (20)^2
# Min; Max; Nom
k1 = (0.2789,0.8927,0.5858)
k2 = (0.1894,0.9331,0.56125) 
k3 = 3.264
k4 = 0.01591
V = 5
q1 = 2.75
q2 = 1
# Inlet (if app.); Steady state conc.; Initial conc.
Ca = (0.3,0.06044,0)
Cb = (7.35,1.72578,3.5)
Cc = (0.00957,0)
Cd = (0.07531,0.0025)


# ODE
function func(dx, x, θ, t) 
    dx[1] = (q1/V)*Ca[1] - x[1]*((q1+θ[3])/V) - x[1]*x[2]*(θ[1]+θ[2])
    dx[2] = (θ[3]/V)*Cb[1] - x[2]*((q1+θ[3])/V) - x[1]*x[2]*(θ[1]+θ[2]) - x[2]*x[3]*k3 - x[2]*x[4]*k4
    dx[3] = -x[3]*((q1+θ[3])/V) + x[1]*x[2]*θ[1] - x[2]*x[3]*k3
    dx[4] = -x[4]*((q1+θ[3])/V) + x[1]*x[2]*θ[2] - x[2]*x[4]*k4
end
problem = ODEProblem(func, [Ca[3]; Cb[3]; Cc[2]; Cd[2]], tSpan, [k1[3];k2[3];q2])


# Runtime
cbs = []
tStops = []
function FaultValve2(α)
    function a!(intregrator)
        intregrator.p[3] = intregrator.p[3] + α*((rand(Float64)*2)-1)
    end
    cb = PeriodicCallback(a!,1)
    push!(cbs, cb)
end
function BreakValve2(t) 
    function c(u,x,intregrator) 
        t == x
    end
    function a!(intregrator)
        intregrator.p[3] = 0
    end
    cb = DiscreteCallback(c,a!) 
    push!(cbs,cb)
    push!(tStops, t)
end
# FaultValve2(0.05)
# BreakValve2(12.5)

# Plot
function Graph() 
    graph = plot(title="Bilinear Reaction Network", xlabel="Time (s)", ylabel="Concentraion (L/mol)", dpi=600)
    uniform01 = Uniform01OrthoPoly(6)
    k1_pce_uniform = convert2affinePCE(k1[1],k1[2],uniform01)
    k1_samples = evaluatePCE(k1_pce_uniform, sampleMeasure(Int(floor(sqrt(trajectories))),uniform01),uniform01)
    k2_pce_uniform = convert2affinePCE(k2[1],k2[2],uniform01)
    k2_samples = evaluatePCE(k2_pce_uniform, sampleMeasure(Int(floor(sqrt(trajectories))),uniform01),uniform01)
    function prob(prob, i, repeat)
        c = floor(sqrt(trajectories))

        prob.p[1] = k1_samples[Int(floor(((i-1)/c))+1)]
        prob.p[2] = k2_samples[Int(mod(i-1, c)+1)]

        return prob
    end
    ensemble = EnsembleProblem(problem, prob_func = prob)
    sim = solve(ensemble, Tsit5(), EnsembleThreads(), trajectories=trajectories, callback=CallbackSet(cbs...), tstops=tStops)
    sum = EnsembleSummary(sim, tSpan[1]:0.1:tSpan[2])
    plot!(graph, sum, fillalpha=0.3, width=0, label=["" "" "" ""])

    function MaximizeOpt(j) 
        function optFunc(x, p) 
            sol = solve(problem, Tsit5(), p=[x[1],x[2],q2], saveat=tSpan[1]:0.1:tSpan[2], callback=CallbackSet(cbs...), tstops=tStops)
            v = sol.u[size(sol.t)[1]]
            return -v[j]
        end
        opt = OptimizationFunction(optFunc)
        optProb = Optimization.OptimizationProblem(opt, [k1[1],k2[2]], [1.0], lb = [k1[1],k2[1]], ub = [k1[2],k2[2]])
        sol = solve(optProb, NLopt.LN_NELDERMEAD())
        return sol.u
    end

    plot!(graph, solve(problem, Tsit5(), p=[MaximizeOpt(2)..., q2], callback=CallbackSet(cbs...), tstops=tStops), legend=:topright, label=["Ca" "Cb" "Cc" "Cd"])
    display(graph)
end
function Save(dir) 
    savefig(string(@__DIR__, "/figs/$dir.png"))
end

Graph()
# Save("stable-poly-chaos")