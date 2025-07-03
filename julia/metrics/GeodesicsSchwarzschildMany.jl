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
_title = "Schwarzschild spacetime (many orbits)"
println(_title)

_legend = :outertopright
rsvalue = 1
args = (rsvalue,)
ϵ = 0.1
#%%

plot(
    aspect_ratio=:equal, 
    grid_ls=:dot,gridcolor=:white,gridwidth=1,gridalpha=0.3,
    dpi=400,
    )

num = 10
dphidot = 0.02 / num
N = Int(2e4); dt = 0.005;
T = N*dt
t = collect(0:dt:T)[1 : end-1]

for j in -num:num
    global v0 = [1.0,-0.5,dphidot*j]; global x0 = [0.0,20.0,0.0];
    
    # @time v4, x4 = Rungekutta(N::Int64,dt::Float64,
        # v0::Vector{Float64}, x0::Vector{Float64}, ; Accel=AccelSchwarzschild,Schwarzschild = true,args=args) 
    # @time v45, x45 = RKF45(N::Int64,dt::Float64,
        # v0::Vector{Float64}, x0::Vector{Float64},
        # ϵ::Float64; Accel=AccelSchwarzschild,Schwarzschild = true,args=args)
    @time v8, x8 = RK8(N::Int64,dt::Float64,
        v0::Vector{Float64}, x0::Vector{Float64}, ; Accel=AccelSchwarzschild,Schwarzschild = true,args=args)
    #%%
    # X4, Y4 = polar(x4)
    # X45, Y45 = polar(x45)
    X8, Y8 = polar(x8) 
    #%%
    # plot!(X4, Y4, label = "RK4")
    # plot!(X45, Y45, label = "RKF45")
    NUM = 2num + 1
    J = j + num + 1
    plot!(X8, Y8, label = "RK8 $j",color=RGB(1,J/NUM,0))#color="cyan"
    #%%

end
function myreplace(str::String)
    str = replace(str,"."=>"_")
    str = replace(str,"-"=>"m")
    return str
end
q = collect(0:0.01:2pi)
plot!(rsvalue * cos.(q),rsvalue * sin.(q), label = L"r = R_s",color="lightgreen")

plot!(
    xlabel =  L"x \,\,[R_s]", ylabel = L"y \,\,[R_s]",
    title = _title,
    legend=_legend,
    xlim=[-2x0[2],x0[2]],
    # annotation = (0, 0, L"f(x,a) = ax"),
)

file_name = ("SchwarzschildMany" * "N$N" * myreplace("dt$(dt)")
            * myreplace("v0X$(v0[1])Y$(v0[2])Z$(v0[3])")
            * myreplace("_x0X$(x0[1])Y$(x0[2])Z$(x0[3])"))
@time savefig(file_name)
println("\nDone")