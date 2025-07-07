#%%
include("RelativityTools.jl")
using .RelativityTools
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
function myreplace(str::String)
    str = replace(str,"."=>"_")
    str = replace(str,"-"=>"m")
    return str
end
function replace_velocity(v0, x0, j)
    g = SchwarzschildMetric(x0)
    psi = pi * j / num
    v0[2] = (-g[1]) * cos(psi)
    v0[3] = sqrt(-g[1]) * sin(psi) / x0[2]
    return v0
end
rsvalue = 1
args = (rsvalue,)
N = Int(2e4)
dt = 0.005
# vphi = 0.02; v0 = [1.0,-0.5,vphi]; x0 = [0.0,20.0,0.0]; default
v0 = [1.0, 0, 0]; x0 = [0.0, 1.5, 0.0];
function sub(v0, x0; num = 10, ϵ = 0.1)
    T = N*dt
    t = collect(0:dt:T)[1 : end-1]
    NUM = 2num + 1
    for j in -num:num
        J = j + num + 1
        # v0[3] = vphi * j / num
        v0 = replace_velocity(v0, x0, j)
        v4, x4 = Rungekutta(N::Int64, dt::Float64, v0::Vector{Float64}, x0::Vector{Float64} ; Accel=AccelSchwarzschild, Schwarzschild = true,args=args)
        X4, Y4 = polar(x4)
        plot!(X4, Y4, label = "RK4 $j", color=RGB(1,J/NUM,0)) #color="cyan"
        # @time v45, x45 = RKF45(N::Int64, dt::Float64,
            # v0::Vector{Float64}, x0::Vector{Float64},
            # ϵ::Float64; Accel=AccelSchwarzschild,Schwarzschild = true,args=args)
        # X45, Y45 = polar(x45)
        # plot!(X45, Y45, label = "RKF45")
        #%%
        # @time v8, x8 = RK8(N::Int64, dt::Float64,
        #     v0::Vector{Float64}, x0::Vector{Float64}, ; Accel=AccelSchwarzschild,Schwarzschild = true,args=args)
        # X8, Y8 = polar(x8) 
        # plot!(X8, Y8, label = "RK8 $j", color=RGB(1,J/NUM,0)) #color="cyan"
    end
    return 0
end

function main(v0, x0)
    _title = "Schwarzschild spacetime (many orbits)"
    println(_title)
    plot(
        title = _title,
        aspect_ratio=:equal, 
        grid_ls=:dot,gridcolor=:white,gridwidth=1,gridalpha=0.3,
        dpi=400,
        )
    q = collect(0:0.01:2pi)
    plot!(rsvalue * cos.(q), rsvalue * sin.(q), label = L"r = R_s", color="lightgreen")

    sub(v0, x0)

    _legend = :outertopright
    plot!(
        xlabel =  L"x \,\,[R_s]", ylabel = L"y \,\,[R_s]",
        legend=_legend,
        xlim=[-2x0[2],2x0[2]],
        ylim=[-2x0[2],2x0[2]],
        # annotation = (0, 0, L"f(x,a) = ax"),
    )

    # file_name = ("1SchwarzschildMany" * "N$N" * myreplace("dt$(dt)")
    #         * myreplace("v0X$(v0[1])Y$(v0[2])Z$(v0[3])")
    #         * myreplace("_x0X$(x0[1])Y$(x0[2])Z$(x0[3])"))
    file_name = "SchwarzschildMany-photon-surface"
    @time savefig(file_name)

    return 0
end

@time main(v0, x0)
