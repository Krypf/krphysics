#%%
include("RelativityTools.jl")
using .RelativityTools
using Plots
# theme(:juno)
using LaTeXStrings
import LinearAlgebra: norm

#%%

function replace_velocity(v0, x0, j, num)
    g = SchwarzschildMetric(x0)
    psi = pi * j / num
    v0[2] = (-g[1]) * cos(psi)
    v0[3] = sqrt(-g[1]) * sin(psi) / x0[2]
    return v0
end

function sub(N, dt, v0, x0; num = 10, ϵ = 0.1)
    NUM = 2num + 1
    for j in -num:num
        v0 = replace_velocity(v0, x0, j, num)
        v4, x4 = Rungekutta(N::Int64, dt::Float64, v0::Vector{Float64}, x0::Vector{Float64} ; Accel=AccelSchwarzschild, Schwarzschild = true,args=args)
        X4, Y4 = polar(x4)
        
        metric_array = metric_test(v4, x4; offset = true)
        println(norm(metric_array[1:end-10]))
        # plot curves
        J = j + num + 1
        # plot!(X4, Y4, label = "RK4 $j π / $num", color=RGB(1,J/NUM,0)) #color="cyan"
    end
    return nothing
end

function main(N, dt, v0, x0)
    _title = "Schwarzschild spacetime (orbits of photons)"
    println(_title)
    plot(
        title = _title,
        aspect_ratio=:equal, 
        grid_ls=:dot,gridcolor=:white,gridwidth=1,gridalpha=0.3,
        dpi=400,
        )
    q = collect(0:0.01:2pi)
    plot!(rsvalue * cos.(q), rsvalue * sin.(q), label = L"r = R_s", color="lightgreen")

    sub(N, dt, v0, x0)

    _legend = :outertopright
    plot!(
        xlabel =  L"x \,\,[R_s]", ylabel = L"y \,\,[R_s]",
        legend=_legend,
        xlim=[-2x0[2],2x0[2]],
        ylim=[-2x0[2],2x0[2]],
        # annotation = (0, 0, L"f(x,a) = ax"),
    )
    file_name = "SchwarzschildMany-photon-surface"
    @time savefig(file_name)

    return 0
end

rsvalue = 1
args = (rsvalue,)
N = Int(2e4)
dt = 0.005
T = N*dt
t = collect(0:dt:T)[1 : end-1]
v0 = [1.0, 0, 0]; x0 = [0.0, 1.5, 0.0];

@time main(N, dt, v0, x0)
