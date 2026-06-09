using OrdinaryDiffEq
using OrdinaryDiffEqLowOrderRK
using Printf

# y' = sin(x)。引数順は (u, p, t): u=y, p=パラメータ(なし), t=独立変数(=x)
f(u, p, t) = sin(t)

u0 = 0.0
xspan = (0.0, 2pi)
prob = ODEProblem(f, u0, xspan)

exact(x) = 1 - cos(x)

sol_rk4  = solve(prob, RK4();  dt = 0.01, adaptive = false)
sol_tsit = solve(prob, Tsit5(); reltol = 1e-8, abstol = 1e-8)

println("="^60)
@printf("%8s %12s %12s %12s %10s %10s\n",
        "x","exact","RK4","Tsit5","errRK4","errTsit5")
for x in range(0, 2pi; length = 9)
    ye = exact(x); y4 = sol_rk4(x); yt = sol_tsit(x)
    @printf("%8.4f %12.8f %12.8f %12.8f %10.2e %10.2e\n",
            x, ye, y4, yt, abs(y4-ye), abs(yt-ye))
end

xs = range(0, 2pi; length = 200)
@printf("max err RK4 (dt=0.01): %.3e\n", maximum(abs.(sol_rk4.(xs) .- exact.(xs))))
@printf("max err Tsit5 (1e-8) : %.3e\n", maximum(abs.(sol_tsit.(xs) .- exact.(xs))))
@printf("RK4 steps  : %d\n", length(sol_rk4.t))
@printf("Tsit5 steps: %d\n", length(sol_tsit.t))