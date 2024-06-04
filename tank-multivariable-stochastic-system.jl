using DifferentialEquations
using Plots; gr()


# Vars
tSpan = (0.0,15.0)
α = (0.5,0.5)
β = 0
γ = (0.22,0.11)
u0 = (10,25)


# Functions
f1(t) = 0
f2(t) = t

function func(du, u, p, t) 
    du[1] = f1(t) - α[1]*u[1] + β*α[2]*u[2] 
    du[2] = f2(t) + α[1]*u[1] - α[2]*u[2]
    du[3] = (1-β)*α[2]*u[2]
end

function noise(du, u, p, t) 
    du[1] = γ[1]*u[1]
    du[2] = γ[2]*u[2]
    du[3] = 0
end

dt = 1/(2^(4))
problem = SDEProblem(func, noise, [u0[1]; u0[2]; 0], tSpan)


# Plot
p = plot()
function PlotInputs() 
    plot!(p, range(tSpan[1],tSpan[2]), f1, label="F1")
    plot!(p, range(tSpan[1],tSpan[2]), f2, label="F2")
end

function PlotVolume()
    sol = solve(problem, EM(), dt=dt)

    plot!(p, sol, label=["Vol. Tank 1" "Vol. Tank 2" "Wash"])
end

# PlotInputs()
PlotVolume()