using DifferentialEquations
using Plots; gr()


# Vars
tSpan = (0.0,25.0)
# Tank 1 rate of flow; Tank 2 rate of flow.
α = (0.9,0.5) 
# Amount returned to system.
β = 0.4
# Initial content in Tank 1; Tank 2; Feed Tank; pH in Tank 1.
u0 = (12.5,25,10,7)

# Functions
# Regulatory fill function 1 (pH)
f1(t) = 0.2t
# pH over time
fp(t) = 0.3

# Regulatory Fill function 2 (feed)
f2(t) = t*t
# Feed consumption function
fe(t) = t


# ODE
function func(du, u, p, t) 
    b = β * ((7-broadcast(abs, 7 - u[5]))/7)

    du[1] = p[3]*p[1]*f1(t) - α[1]*u[1] + b*α[2]*u[2] 
    du[2] = p[2]*p[1]*f2(t) + α[1]*u[1] - α[2]*u[2]
    du[3] = (1-b)*α[2]*u[2]
    du[4] = p[2]*p[1]*f2(t) - fe(t)
    du[5] = p[3]*p[1]*f1(t) - fp(t)
end

problem = ODEProblem(func, [u0[1]; u0[2]; 0; u0[3]; u0[4]], tSpan, [1;0;0])


# Dynamics
cbs = []
tStops = []
function blockFlow(x) 
    function c(u,t,intregrator) 
        t == x
    end
    function a!(intregrator)
        intregrator.p[1] = 0
        intregrator.p[2] = 0
    end
    cb = DiscreteCallback(c,a!) 
    push!(cbs,cb)
    push!(tStops, x)
end

function regulateFeed()
    function c0(u,t,intregrator)
        u[4] - (u0[3]/2)
    end
    function a0!(intregrator)
        intregrator.p[2] = 1
    end
    cb0 = ContinuousCallback(c0, a0!)
    push!(cbs,cb0)

    function c1(u,t,intregrator)
        u[4] - u0[3]
    end
    function a1!(intregrator)
        intregrator.p[2] = 0
    end
    cb1 = ContinuousCallback(c1, a1!)
    push!(cbs,cb1)
end

function regulatePH()
    function c0(u,t,intregrator)
        u[5] - 6
    end
    function a0!(intregrator)
        intregrator.p[3] = 1
    end
    cb0 = ContinuousCallback(c0, a0!)
    push!(cbs,cb0)

    function c1(u,t,intregrator)
        u[5] - 8
    end
    function a1!(intregrator)
        intregrator.p[3] = 0
    end
    cb1 = ContinuousCallback(c1, a1!)
    push!(cbs,cb1)
end

# Plot
graph = plot(xlims=tSpan, ylims=(0,Inf))
function PlotFeed()

end

function PlotInputs() 
    plot!(graph, range(tSpan[1],tSpan[2]), f1, label="F1")
    plot!(graph, range(tSpan[1],tSpan[2]), f2, label="F2")
end

function PlotVolume()
    sol = solve(problem, Tsit5(), callback=CallbackSet(cbs...), tstops=tStops)

    plot!(graph, sol.t, sol[1,:], label="Vol. Tank 1")
    plot!(graph, sol.t, sol[2,:], label="Vol. Tank 2")
    plot!(graph, sol.t, sol[3,:], label="Wash")
end

function PlotOutputs()
    sol = solve(problem, Tsit5(), callback=CallbackSet(cbs...), tstops=tStops)
    
    plot!(graph, sol.t, sol[4,:], label="Feed")
    plot!(graph, sol.t, sol[5,:], label="pH")
end


# Runtime
blockFlow(15)
regulateFeed()
regulatePH()

# PlotInputs()
PlotOutputs()
# PlotVolume()