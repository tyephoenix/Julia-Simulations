using DifferentialEquations
using Plots; gr()
# a = 1
# u0 = 5 
# f(u,p,t) = a*u
# tSpan = (0.0, 3.0)

# prob = ODEProblem(f, u0, tSpan, saveat=0.05)
# sol = solve(prob, Tsit5())

# Vars
tSpan = range(0.0,10.0)

a1 = 1
a2 = 1
b = 1

f1(t) = 0
f2(t) = 0

v1(t) = f1(t) + 5
p1(t) = a1*v1(t)

v2(t) = p1(t) + f2(t) + 10
p2(t) = a2*v2(t)

R(t) = b*p2(t)



p = plot()
function PlotInputs() 
    plot!(p, tSpan, f1, label="F1")
    plot!(p, tSpan, f2, label="F2")
end

function PlotPhi() 
    plot!(p, tSpan, p1, label="P1")
    plot!(p, tSpan, p2, label="P2")
end

# PlotInputs()
PlotPhi()

# function Tank2(in1, in2) 
#     # prob = ODEProblem(in1, h1, tSpan, saveat=0.05)
    
#     sol = solve(prob, Tsit5())
#     # plot(sol)
# end

# Tank2(f1, f2)