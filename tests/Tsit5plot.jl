using OrdinaryDiffEq
using OrdinaryDiffEqLowOrderRK   # RK4 用（あなたの環境で必要だった）
using PyPlot
using Printf

# y' = sin(x), y(0)=0  → 解析解 y = 1 - cos(x)
f(u, p, t) = sin(t)
u0 = 0.0
xspan = (0.0, 2pi)
prob = ODEProblem(f, u0, xspan)
exact(x) = 1 - cos(x)

sol_rk4  = solve(prob, RK4();   dt = 0.01, adaptive = false)
sol_tsit = solve(prob, Tsit5(); reltol = 1e-8, abstol = 1e-8)

# 評価用の細かい x グリッド
xs = range(0, 2pi; length = 400)
ye = exact.(xs)
y4 = sol_rk4.(xs)
yt = sol_tsit.(xs)
err4 = abs.(y4 .- ye)
errt = abs.(yt .- ye)

# ---- 図1: 解の重ね描き ----
figure(figsize=(7,5))
plot(xs, ye, "k-",  lw=2,            label="exact  1 - cos(x)")
plot(xs, y4, "C0--", lw=1.5,         label="RK4 (dt=0.01)")
plot(xs, yt, "C3:",  lw=2,           label="Tsit5 (1e-8)")
xlabel("x"); ylabel("y(x)")
title("y' = sin(x): numerical vs exact")
legend(); grid(true, alpha=0.3)
tight_layout()
savefig("tests/sol_compare.png", dpi=150)

# ---- 図2: 数値誤差（片対数）----
figure(figsize=(7,5))
semilogy(xs, err4 .+ 1e-18, "C0-", lw=1.5, label="RK4 error")
semilogy(xs, errt .+ 1e-18, "C3-", lw=1.5, label="Tsit5 error")
xlabel("x"); ylabel("|numerical - exact|")
title("Absolute error vs x")
legend(); grid(true, which="both", alpha=0.3)
tight_layout()
savefig("tests/error_compare.png", dpi=150)

@printf("max err RK4   : %.3e  (steps %d)\n", maximum(err4), length(sol_rk4.t))
@printf("max err Tsit5 : %.3e  (steps %d)\n", maximum(errt), length(sol_tsit.t))
println("saved: tests/sol_compare.png, tests/error_compare.png")