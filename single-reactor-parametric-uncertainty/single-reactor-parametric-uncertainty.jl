using DifferentialEquations
using Distributions
using StatsPlots; gr()
using Optimization
using OptimizationNLopt


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

k1_dist = Uniform(k1[1]*10, k1[2]*10)
k2_dist = Uniform(k2[1]*10, k2[2]*10)

# ODE
function func(dx, x, θ, t) 
    dx[1] = (q1/V)*Ca[1] - x[1]*((q1+θ[3])/V) - x[1]*x[2]*(θ[1]+θ[2])
    dx[2] = (θ[3]/V)*Cb[1] - x[2]*((q1+θ[3])/V) - x[1]*x[2]*(θ[1]+θ[2]) - x[2]*x[3]*k3 - x[2]*x[4]*k4
    dx[3] = -x[3]*((q1+θ[3])/V) + x[1]*x[2]*θ[1] - x[2]*x[3]*k3
    dx[4] = -x[4]*((q1+θ[3])/V) + x[1]*x[2]*θ[2] - x[2]*x[4]*k4
end
problem = ODEProblem(func,[Ca[3]; Cb[3]; Cc[2]; Cd[2]], tSpan, [k1[3];k2[3];q2])


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
    map = Dict()
    function prob(prob, i, repeat)
        c = floor(sqrt(trajectories))

        k1_x = floor(((i-1)/c))*((maximum(k1_dist) - minimum(k1_dist)) / c) + minimum(k1_dist)
        k1_p = pdf(k1_dist, k1_x)
        prob.p[1] = k1_x

        k2_x = (mod(i-1, c)+1)*((maximum(k2_dist) - minimum(k2_dist)) / c) + minimum(k2_dist)
        k2_p = pdf(k2_dist, k2_x)
        prob.p[2] = k2_x

        map[i] = (k1_x, k2_x)

        return prob
    end
    ensemble = EnsembleProblem(problem, prob_func = prob)
    sim = solve(ensemble, Tsit5(), EnsembleThreads(), trajectories=trajectories, callback=CallbackSet(cbs...), tstops=tStops)
    sum = EnsembleSummary(sim, tSpan[1]:0.1:tSpan[2])
    plot!(graph, sum, fillalpha=0.3, width=1, legend=:topright, label=["Ca" "Cb" "Cc" "Cd" "q2"])
    display(graph)

    function Maximize(j) 
        arr = sim
        l = size(arr)[2] - 2
        hV = 0
        tr = 0
        for i in 1:trajectories
            v = arr[j,l,i]
            if v > hV
                hV = v
                tr = i
            end
        end
        return map[tr]
    end
    println(Maximize(4))
end
function Save(dir) 
    savefig(string(@__DIR__, "/figs/$dir.png"))
end

Graph()
Save("broken-faulty")