using DifferentialEquations
using QuadGK;
using Plots; gr()


# Vars
tSpan = (0.0,15.0)
a = (1,1)
b = 0
v1_0 = 5
v2_0 = 10


# Functions
f1(t) = 0
f2(t) = 0

function func(du, u, p, t) 
    du[1] = f1(t) - a[1]*u[1] + b*a[2]*u[2] 
    du[2] = f2(t) + a[1]*u[1] - a[2]*u[2] - b*a[2]*u[2]
    du[3] = (1-b)*a[2]*u[2]
end

function noise(du, u, p, t) 
    du[1] = 0.3
    du[2] = 0.3
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