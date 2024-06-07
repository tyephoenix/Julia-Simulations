using DifferentialEquations
using Plots; gr()
using Distributions
using Optimization
using OptimizationNLopt


# Vars
tSpan = (0.0,20.0)
trajectories = 100
# Pump 1 rate of flow of effluent : Tank 2 rate of flow of titrant
Pa(t) = 0
Pb(t,p) = 1.1*p

α = (0.5)
γ = (0,100)

Cmab = 1.345
Cimp = 3.5

# Distributions
impuritiesDist = JohnsonSU(50,5,3,1.5)
mabDist = JohnsonSU(65,3,3,1.5)
startingMix = Cmab / (Cmab + Cimp)
mix = MixtureModel([mabDist,impuritiesDist], [startingMix, 1-startingMix])
# erlang = Uniform(1, 100)

function loop(dx, x, θ, t)
    V = x[1] + x[2] + x[3]
    dx[1] = Pa(t)*Cmab
    dx[2] = Pa(t)*Cimp

    dx[3] = -θ[2]*Pb(t,0.2*x[3]) + 0.1*x[3]

    function Column(p,a)
        impProb = cdf(impuritiesDist,100)-cdf(impuritiesDist,p)
        imp = impProb * Cimp
        mabProb = cdf(mabDist,100)-cdf(mabDist,p)
        mab = mabProb * Cmab
        dx[1] = dx[1] - a*(x[1]-mab)
        dx[2] = dx[2] - a*(x[2]-imp)
    end 

    if (t >= 5)
        Column(θ[1], α[1])
    end
end
problem = ODEProblem(loop, [Cmab; Cimp;5], tSpan, [γ[1];0])


# Runtime
cbs = []
tStops = []
function RegulatePH()
    function c0(u,t,intregrator)
        u[3] - 7
    end
    function a0!(intregrator)
        intregrator.p[2] = 1
    end
    cb0 = ContinuousCallback(c0, a0!)
    push!(cbs,cb0)

    function c1(u,t,intregrator)
        u[3] - 5
    end
    function a1!(intregrator)
        intregrator.p[2] = 0
    end
    cb1 = ContinuousCallback(c1, a1!)
    push!(cbs,cb1)
end

RegulatePH()

# Plots
function GraphDist()
    graph = plot(xlims=(1,Inf), ylims=(0,20), title="Distribution of Particle Size (in)", xlabel="μm", ylabel="%")
    plot!(graph, 1:0.1:100, (t) -> 100*pdf(impuritiesDist,t), label="Impurities")
    plot!(graph, 1:0.1:100, (t) -> 100*pdf(mabDist,t), label="mAB")
    graph0 = plot(xlims=(1,Inf), ylims=(0,20), title="Distribution of Particle Size (in)", xlabel="μm", ylabel="%")
    plot!(graph0, 1:0.1:100, (t) -> 100*pdf(mix,t), label="Total")
    display(graph)
    display(graph0)
end

function GraphConcentration(γ; graph=plot(title="Concentration", label=["C(mAB)" "C(imp)" "pH"]))
    sol = solve(problem, Tsit5(), p=[γ,0], callback=CallbackSet(cbs...), tstops=0.0:0.1:25)
    plot!(graph, sol, label=["C(mAB)" "C(imp)" "pH"], xlabel="Time (t)", ylabel="mol")
    display(graph)
end

function GraphEnsemble()
    function prob(prob, i, repeat)
        l = (i-1)*((γ[2]-γ[1])/trajectories)+γ[1]
        prob.p[1] = l
        return prob
    end
    ensemble = EnsembleProblem(problem, prob_func = prob)
    sol = solve(ensemble, Tsit5(), callback=CallbackSet(cbs...), trajectories=trajectories,tstops=0.0:0.1:25)
    sum = EnsembleSummary(sol, tSpan[1]:0.1:tSpan[2])
    graph = plot(title="Optimal Concentration", label=["C(mAB)" "C(imp)" "pH"])
    plot!(graph, sum, width=0, legend=:topright,label=["" ""])
    display(graph)

    # Optimization
    function opt(x, p) 
        sol = solve(problem, Tsit5(), p=[x[1],0], saveat=0:0.1:20)
        v = sol.u[size(sol.t)[1]]
        return v[2] - v[1]
    end
    x0 = [γ[1]]
    funca = OptimizationFunction(opt)
    probl = Optimization.OptimizationProblem(funca, x0, [1.0], lb = γ[1], ub = γ[2])
    sol = solve(probl, NLopt.LN_NELDERMEAD())
    GraphConcentration(sol.u[1], graph=graph)
end

# ColumnEvent(10, 70)

GraphDist()
GraphEnsemble()
# GraphConcentration()
