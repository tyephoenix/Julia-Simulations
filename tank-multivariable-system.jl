using DifferentialEquations
using QuadGK;
using Plots; gr()


# Vars
tSpan = (0.0,15.0)
a = (0.5,0.3)
b = 1
v1_0 = 5
v2_0 = 10


# Functions
f1(t) = 2*t
f2(t) = 2

function func(du, u, p, t) 
    du[1] = f1(t) - a[1]*u[1] + b*a[2]*u[2] 
    du[2] = f2(t) - a[2]*u[2] - b*a[2]*u[2]
end

problem = ODEProblem(func, [v1_0; v2_0], tSpan)


# Plot
p = plot()
function PlotInputs() 
    plot!(p, range(0,15), f1, label="F1")
    plot!(p, range(0,15), f2, label="F2")
end


function PlotVolume()
    sol = solve(problem, Tsit5())

    plot!(p, sol)
end

PlotInputs()
PlotVolume()