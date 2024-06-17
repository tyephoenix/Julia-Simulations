using Plots; 
using DifferentialEquations
using Optimization
using OptimizationNLopt
using LinearAlgebra
using LaTeXStrings
specific_plot_kwargs = Dict{Symbol,Any}(
    # plotattr(:Plot)
    # plotattr(:Axis)
    :html_output_format => :png, # can try {:png,:svg} here
    :dpi => 300,                 # DPI=300 is minimum for good quality and print
    :grid => true,              # eliminate grid lines
    :gridstyle => :dash,         # if you actually want grid lines, choose from [:auto, :solid, :dash, :dot, :dashdot, :dashdotdot]
    :gridlinewidth => 0.5,       # set to 0 to make the grid and axes outlines disappear
    :framestyle => :box,         # add matlab-esque boundary box on plot
    :size => (450,450),          # typically (600,450) looks fine
    :topmargin => 0Plots.mm,
    :bottommargin => 0Plots.mm,
    :leftmargin => 0Plots.mm,
    :rightmargin => 5Plots.mm,
    :background_color_subplot => :white,
    :legend_background_color => nothing, # removes white background and adds translucence to legend, useful if you still want to see what's goin on behind a plot legend
    :legend_foreground_color => nothing, # removes outer black border of legend, just remove it alongside the white background of default legend to improve visibility
    :fontfamily => "computer modern", # DO NOT DEVIATE FROM THIS IF USING GR()
    :guidefontfamily => :match,       # match any fontfamily on title labels to axis labels
    :tickfontfamily => :match,       # match any fontfamily on title labels to ticks and tick labels
    :titlefontsize => 6, # for title fontsize
    :guidefontsize => 6, # for x and y labels
    :legendfontsize => 5, # for legend labels
    :tickfontsize => 6, # for axis ticks
    # :title => L"Sample text: $- \frac{{\hbar ^2}}{{2m}}\frac{{\partial ^2 \psi (x,t)}}{{\partial x^2 }} + U(x)\psi (x,t) = i\hbar \frac{{\partial \psi (x,t)}}{{\partial t}}$",
    # :xlabel => L"Sample text: $- \frac{{\hbar ^2}}{{2m}}\frac{{\partial ^2 \psi (x,t)}}{{\partial x^2 }} + U(x)\psi (x,t) = i\hbar \frac{{\partial \psi (x,t)}}{{\partial t}}$",
    # :ylabel => L"Sample text: $- \frac{{\hbar ^2}}{{2m}}\frac{{\partial ^2 \psi (x,t)}}{{\partial x^2 }} + U(x)\psi (x,t) = i\hbar \frac{{\partial \psi (x,t)}}{{\partial t}}$",

    # plotattr(:Series)
    # :label => L"Sample legend text: $a=\frac{{\hbar ^2}}{{2m}}$",
    :legend => :topleft,                # placement of legend, see `:(outer ?)(top/bottom ?)(right/left ?)`
    :linewidth => 1.0,
    :linestyle => :solid,               # see also [:auto, :solid, :dash, :dot, :dashdot, :dashdotdot]
    :linealpha => 1,
    # :linecolor => 1,
    # marker
    :markersize => 0.2,
    # :markershape => :star4,             # can choose from many, personal favorites = {:circle,:star4,:utriangle,:cross}
    :markeralpha => 1.0,
    :markercolor => nothing,
    # for that pesky marker outline on the markers
    :markerstrokewidth => 0.0,           # VERY IMPORTANT, make the marker outlines disappear -> 0
    :markerstrokestyle => :solid,        # the outline style on the marker
    :markerstrokecolor => nothing,       # can force this to be the same color as the marker, but still want outline
    :markerstrokealpha => 1,             # opacity of marker outline
)
gr(; specific_plot_kwargs...) # set gr backend






# Variables
tSpan = (0, 60)

Ca_0 = 20
Cb_0 = 30
Cc_0 = 0
k = 0.02


# True System
function reaction(du, u, p, t)
    du[1] = -k * u[1] * u[2]
    du[2] = -k * u[1] * u[2] + 0.5
    du[3] = k * u[1] * u[2]
end


# State Estimator
const L = [[0.023 0.0; 0.53 0.026; 0.44 0.50],0.0]
estX0 = [0,0,30]

cbs = []
tStops = []
function OptimizeObserver_NLOPT(interval; maxtime=30) 
    function observer(du, u, p, t)
        ym = [p[2]; p[3]]
        y_hat = [u[1] + u[2]; u[3]]
        mat = p[1] * (ym-y_hat)
    
        du[1] = -k * u[1] * u[2] + mat[1] 
        du[2] = -k * u[1] * u[2] + mat[2]
        du[3] = k * u[1] * u[2] + mat[3]
    end
    observer = ODEProblem(observer, estX0, tSpan)

    i = 1
    while i*interval <= tSpan[2]
        timeIndex = i*interval
        preTimeIndex = (i-1)*interval
        function c(u,t,intregrator) 
            if (t == timeIndex)
                println("Optimizing at $(timeIndex)...")
            end
            t == timeIndex
        end
        function a!(intregrator)
            ti = time()
            map = Dict()
            function opt(x, p) 
                ym = [intregrator.u[1] + intregrator.u[2], intregrator.u[3]]
                tL = [x[1] x[2]; x[3] x[4]; x[5] x[6]]
                tv  = intregrator.u
                eSol = solve(observer, Tsit5(), p=[tL, ym...], saveat=preTimeIndex:0.1:timeIndex)
                ev = eSol(timeIndex)
                err = tv-ev
                sq = err[1]^2 + err[2]^2 + err[3]^2
                if (L[2] != 0)
                    peSol = solve(observer, Tsit5(), p=[L[1], ym...], saveat=preTimeIndex:0.1:timeIndex)
                    pev = peSol(timeIndex)
                    peErr = tv-pev
                    peSq =  peErr[1]^2 + peErr[2]^2 + peErr[3]^2
                    map[tL] = [sq, peSq]
                else 
                    map[tL] = [sq, 0]
                end
                return sq
            end
            x0 = [0,0,0,0,0,0]
            p = [1.0]
            funca = OptimizationFunction(opt)
            probl = Optimization.OptimizationProblem(funca, x0, p, lb = [0,0,0,0,0,0], ub = [1,1,1,1,1,1])
            sol = solve(probl, NLopt.LN_NELDERMEAD(), maxtime=maxtime)
            nL = [sol.u[1] sol.u[2]; sol.u[3] sol.u[4]; sol.u[5] sol.u[6]]
            s = map[nL]
            if (L[2] == 0) 
                v = [sol.u..., s]
                L[1] = [v[1] v[2]; v[3] v[4]; v[5] v[6]]
                L[2] = v[7][1]
            else 
                tw = s[1] + s[2]
                itw = (tw/s[2]) + (tw/(s[1]))
                v = ((tw/s[2])/itw)*L[1] + ((tw/(s[1]))/itw)*(nL)
                L[1] = v
                acc = ((tw/s[2])/itw)*s[2] + ((tw/(s[1]))/itw)*(s[1])
                L[2] = acc
            end
            println("OptimizeObserver_NLOPT: Discrete Event at $(timeIndex). Elapsed time = $(time() - ti).")
        end
        cb = DiscreteCallback(c,a!) 
        push!(cbs,cb)
        push!(tStops, timeIndex)
        i = i + 1
    end
end


# Solving
# OptimizeObserver_NLOPT(2, maxtime=10)

problem = ODEProblem(reaction, [Ca_0, Cb_0, Cc_0], tSpan)
sol = solve(problem, Tsit5(), callback=CallbackSet(cbs...), tstops=tStops)

println(L[1])

function estimator(du, u, p, t)
    v = sol(t)
    ym = [v[1] + v[2]; v[3]]
    y_hat = [u[1] + u[2]; u[3]]
    mat = p[1] * (ym-y_hat)

    du[1] = -k * u[1] * u[2] + mat[1] 
    du[2] = -k * u[1] * u[2] + mat[2]
    du[3] = k * u[1] * u[2] + mat[3]
end
estProblem = ODEProblem(estimator, estX0, tSpan)
estSol = solve(estProblem, Tsit5(), p=[L[1]])

# Plot
function Save(dir) 
    savefig(string(@__DIR__, "/figs/$dir.png"))
end

graph = plot(sol, label=[L"$C_a$" L"$C_b$" L"$C_c$"], xlabel="Time (t)", ylabel="Concentration (mol)", title="Bilinear Reaction State Estimation")
plot!(graph, estSol, linestyle=:dash, label=[L"$C_a$ Estimation" L"$C_b$ Estimation" L"$C_c$ Estimation"], xlabel="Time (t)")
display(graph)

Save("state-estimation-100x")