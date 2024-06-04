using DifferentialEquations
using Plots; gr()


# Vars
tSpan = (0.0,35.0)
α = (0.5,0.5)
β = 0.5
u0 = (12.5,25)


# Functions
f1(t) = 0
f2(t) = 0

function func(du, u, p, t) 
    du[1] = f1(t) - α[1]*u[1] + β*α[2]*u[2] 
    du[2] = f2(t) + α[1]*u[1] - α[2]*u[2]
    du[3] = (1-β)*α[2]*u[2]
end

problem = ODEProblem(func, [u0[1]; u0[2]; 0; 0], tSpan)




# Plot
p = plot()
function PlotInputs() 
    plot!(p, range(tSpan[1],tSpan[2]), f1, label="F1")
    plot!(p, range(tSpan[1],tSpan[2]), f2, label="F2")
end

function PlotVolume()
    sol = solve(problem, Tsit5())

    plot!(p, sol, label=["Vol. Tank 1" "Vol. Tank 2" "Wash" "Total"])
end

# PlotInputs()
PlotVolume()