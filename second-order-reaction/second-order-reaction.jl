using DifferentialEquations
using Plots; gr()


# Vars
tSpan = (0.0,10.0)
u0 = (1,1)
# Min,Max,True
k = [0.1,0.9,0.66]
Rab = 0.2
# Ensemble
trajectories = 1000

function func(du, u, p, t) 
    du[1] = -p[1]*u[1]*Rab
    du[2] = -(1-p[1])*u[2]*(1-Rab)
    du[3] = -(du[1] + du[2])
end
problem = ODEProblem(func, [u0[1]; u0[2]; 0], tSpan, [k[3]])


# Plot
p = plot(title="Second-Order Reaction", ylabel="Arbitrary Amount", dpi=600)
function Graph(x)
    sol = solve(problem, Tsit5(), p=[x])

    plot!(p, sol, legend=:topleft, xlabel="Time (s)", label=["A" "B" "Product"])
    display(p)
end
function Ensemble() 
    function prob(prob,i,repeat)
        j = (i-1)*((k[2]-k[1])/trajectories) + k[1]
        prob.p[1] = j
        return prob
    end
    ensemble = EnsembleProblem(problem, prob_func=prob)
    sim = solve(ensemble, Tsit5(), EnsembleThreads(), trajectories=trajectories)
    sum = EnsembleSummary(sim, tSpan[1]:0.1:tSpan[2])

    plot!(p, sum, fillalpha=0.3, xlabel="Time (s)", width=0, label=["" "" ""])
    display(p)
end
function Save(dir) 
    savefig(string(@__DIR__, "/figs/$dir.png"))
end


# Optimization
using Optimization
using OptimizationNLopt

function opt(x, p) 
    sol = solve(problem, Tsit5(), p=x[1], saveat=0:0.1:20)
    v = sol.u[size(sol.t)[1]]
    return v[1]+v[2]
end
x0 = [k[1]]
funca = OptimizationFunction(opt)
probl = Optimization.OptimizationProblem(funca, x0, [1.0], lb = k[1], ub = k[2])
sol = solve(probl, NLopt.LN_NELDERMEAD())

print(sol)
Ensemble()
Graph(sol.u[1])
Save("optimization")