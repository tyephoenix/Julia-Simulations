using DifferentialEquations
using QuadGK;
using Plots; gr()

# Vars
tSpan = range(0.0,2.0)
a1 = 0.5
a2 = 0.01
b = 0

# Functions
f1(t) = 0
f2(t) = 0

v1_0 = 5
v1(t) = quadgk(f1, 0, t)[1] + v1_0
p1(t) = a1*v1(t)

v2_0 = 10
v2(t) = quadgk(p1, 0, t)[1] + quadgk(f2, 0, t)[1] + v2_0
p2(t) = a2*v2(t)

R(t) = b*p2(t)

# Tanks
t1(t) = v1(t) - quadgk(p1, 0, t)[1] + quadgk(R, 0, t)[1]
t2(t) = v2(t) - quadgk(p2, 0, t)[1]

# Plot
p = plot()
function PlotInputs() 
    plot!(p, tSpan, f1, label="F1")
    plot!(p, tSpan, f2, label="F2")
end

function PlotPhi() 
    plot!(p, tSpan, p1, label="P1")
    plot!(p, tSpan, p2, label="P2")
end

function PlotVolume()
    plot!(p, tSpan, t1, label="Vol. Tank 1")
    plot!(p, tSpan, t2, label="Vol. Tank 2")
end

# PlotInputs()
# PlotPhi()
PlotVolume()