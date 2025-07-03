#%%
include("Metrics.jl") # This file is a self-made one. Spacetime Metrics 
include("Accels.jl") # This file is a self-made one. Accelerations
include("RKS.jl") # This file is a self-made one. Runge-kutta methods
using Plots
theme(:juno)
using LaTeXStrings
#%%
function polar(x)
    r = x[2,:]; ph = x[3,:]
    X = @. r * cos(ph)
    Y = @. r * sin(ph)
    return X, Y
end
v0 = [1.0,-0.5,0.01]; x0 = [0.0,20.0,0.0]; 
N = Int(1e5); dt = 0.001;
xstart = -2x0[2]; _legend = :best; _left_margin = 5 * Plots.mm;
_xlabel = L"x \,\,[2M]"; _ylabel = L"y \,\,[2M]";
# _title = "Kerr spacetime (almost symmetric orbit, up spin)"; avalue = 0.5; Mvalue = 0.5;
# _title = "Kerr spacetime (almost symmetric orbit, no spin)"; avalue = 0; Mvalue = 0.5;
# _title = "Kerr spacetime (almost symmetric orbit, down spin)"; avalue = -0.5; Mvalue = 0.5;
################
xstart = -x0[2]; _legend = :bottomleft; _xlabel = L"x \,\,[M] \,(M = 1)"; _ylabel = L"y \,\,[M]";
_title = "Kerr spacetime (Larger mass, up spin)"; avalue = 0.5; Mvalue = 1; _left_margin = 0 * Plots.mm
# _title = "Kerr spacetime (Larger mass, down spin)"; avalue = -0.5; Mvalue = 1; _left_margin = 5 * Plots.mm
# _title = "Kerr spacetime (Larger mass, no spin)"; avalue = 0; Mvalue = 1; _left_margin = 5 * Plots.mm
################
# v0 = [1.0,-0.5,0.0001]; x0 = [0.0,200.0,0.0]; 
# xstart = -x0[2]; _legend = :bottomright;_xlabel = L"x \,\,[2M]"; _ylabel = L"y \,\,[2M]";# _left_margin = 5 * Plots.mm
# N = Int(1e5); dt = 0.01;
# _title = "Kerr spacetime (flyby, no spin)"; avalue = 0; Mvalue = 0.5;
# _title = "Kerr spacetime (flyby, up spin)"; avalue = 0.5; Mvalue = 0.5;
# _title = "Kerr spacetime (flyby, down spin)"; avalue = -0.5; Mvalue = 0.5;

# N = Int(1e5); dt = 0.001;
# N = Int(1e6); dt = 0.001;
################
# N = Int(5e6); dt = 0.3;
# N = Int(1e6); dt = 0.0001;
# v0 = [1.0,0,5.565782794195258e-7]; x0 = [0.0,2.3640692681965724e7,0.0]; 
println(_title)
T = N*dt
t = collect(0:dt:T)[1 : end-1]
#%%
args = (avalue,Mvalue)
@time v4, x4 = Rungekutta(N::Int64,dt::Float64,
    v0::Vector{Float64}, x0::Vector{Float64}, ; Accel=AccelKerr,Kerr = true,args=args)
ε = 0.1
@time v45, x45 = RKF45(N::Int64,dt::Float64,
    v0::Vector{Float64}, x0::Vector{Float64},
    ε::Float64; Accel=AccelKerr,Kerr = true,args=args)
# show(IOContext(stdout, :limit => true), "text/plain", x45)

@time v8, x8 = RK8(N::Int64,dt::Float64,
    v0::Vector{Float64}, x0::Vector{Float64}, ; Accel=AccelKerr,Kerr = true,args=args)
# %%
X4, Y4 = polar(x4)
X45, Y45 = polar(x45)
X8, Y8 = polar(x8) 
#%%
plot(X4, Y4, label = "RK4")
plot!(X45, Y45, label = "RKF45")
plot!(X8, Y8, label = "RK8")
plot!(aspect_ratio=:equal, 
    grid_ls=:dot,gridcolor=:white,gridwidth=1,gridalpha=0.3,
    dpi=400)
plot!(xlim=[xstart,x0[2]])

function myreplace(str::String)
    # str = "x$x"
    str = replace(str,"."=>"_")
    str = replace(str,"-"=>"m")
    return str
end
rkvalue = Mvalue + sqrt(Mvalue^2-avalue^2)
q = collect(0:0.01:2pi)
plot!(rkvalue * cos.(q),rkvalue * sin.(q), label = L"r = M+\sqrt{M^2-a^2}")

plot!(
    xlabel = _xlabel, ylabel = _ylabel,
    title = _title,
    legend=_legend,
    left_margin = _left_margin
    # annotation = (0, 0, L"f(x,a) = ax"),
)

file_name = (join(replace(split(_title, ""), " "=>"_", "(" =>"",  ")"=>"", ","=>""))
            * myreplace("a$(avalue)")
            * myreplace("M$(Mvalue)")
            * "N$N" * myreplace("dt$(dt)")
            * myreplace("v0X$(v0[1])Y$(v0[2])Z$(v0[3])")
            * myreplace("_x0X$(x0[1])Y$(x0[2])Z$(x0[3])"))
@time savefig(file_name)
println("\nDone")
