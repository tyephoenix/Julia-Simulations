using DifferentialEquations
using QuadGK;
using Plots; gr()


# Vars
tSpan = (0.0,15.0)
α = (1,1)
β = 0
γ = (0.40,0.22)
v1_0 = 5
v2_0 = 10


# Functions
f1(t) = 2*t*t
f2(t) = 2

function func(du, u, p, t) 
    du[1] = f1(t) - α[1]*u[1] + β*α[2]*u[2] 
    du[2] = f2(t) + α[1]*u[1] - α[2]*u[2] - β*α[2]*u[2]
    du[3] = (1-β)*α[2]*u[2]
end

function noise(du, u, p, t) 
    du[1] = γ[1]*u[1]
    du[2] = γ[2]*u[2]
    du[3] = 0
end

dt = 1/(2^(4))
problem = SDEProblem(func, noise, [v1_0; v2_0; 0], tSpan)


# Plot
p = plot()
function PlotInputs() 
    plot!(p, range(0,15), f1, label="F1")
    plot!(p, range(0,15), f2, label="F2")
end

function PlotVolume()
    sol = solve(problem, EM(), dt=dt)

    plot!(p, sol, label=["Vol. Tank 1" "Vol. Tank 2" "Wash"])
end

# PlotInputs()
PlotVolume()