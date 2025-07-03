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
_title = "Schwarzschild spacetime (almost symmetric orbit)"
v0 = [1.0,-0.5,0.01]; x0 = [0.0,20.0,0.0]; 
N = Int(1e5); dt = 0.001;
_legend = :best
# _title = "Schwarzschild spacetime (flyby orbit)"
# v0 = [1.0,-0.5,0.0001]; x0 = [0.0,200.0,0.0]; 
# N = Int(1e6); dt = 0.001;
# _legend = :bottomright
################
# N = Int(5e6); dt = 0.3;
# N = Int(1e6); dt = 0.0001;
# v0 = [1.0,0,5.565782794195258e-7]; x0 = [0.0,2.3640692681965724e7,0.0]; 
println(_title)
T = N*dt
t = collect(0:dt:T)[1 : end-1]
#%%
rsvalue = 1
args = (rsvalue,)
# include("RK4.jl") # This file is a self-made one. 4th order Runge-kutta method.
@time v4, x4 = Rungekutta(N::Int64,dt::Float64,
    v0::Vector{Float64}, x0::Vector{Float64}, ; Accel=AccelSchwarzschild,Schwarzschild = true,args=args)
# include("RKF45.jl") 
ϵ = 0.1
@time v45, x45 = RKF45(N::Int64,dt::Float64,
    v0::Vector{Float64}, x0::Vector{Float64},
    ϵ::Float64; Accel=AccelSchwarzschild,Schwarzschild = true,args=args)
# show(IOContext(stdout, :limit => true), "text/plain", x)

@time v8, x8 = RK8(N::Int64,dt::Float64,
    v0::Vector{Float64}, x0::Vector{Float64}, ; Accel=AccelSchwarzschild,Schwarzschild = true,args=args)
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
plot!(xlim=[-2x0[2],x0[2]])

function myreplace(str::String)
    # str = "x$x"
    str = replace(str,"."=>"_")
    str = replace(str,"-"=>"m")
    return str
end
q = collect(0:0.01:2pi)
plot!(rsvalue * cos.(q),rsvalue * sin.(q), label = L"r = R_s")

plot!(
    xlabel =  L"x \,\,[R_s]", ylabel = L"y \,\,[R_s]",
    title = _title,
    legend=_legend,
    # annotation = (0, 0, L"f(x,a) = ax"),
)

file_name = ("Schwarzschild" * "N$N" * myreplace("dt$(dt)")
            * myreplace("v0X$(v0[1])Y$(v0[2])Z$(v0[3])")
            * myreplace("_x0X$(x0[1])Y$(x0[2])Z$(x0[3])"))
@time savefig(file_name)
println("\nDone")