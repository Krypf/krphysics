#%% Runge-kutta methods
function Rungekutta(N::Int64,dt::Float64,v0::Vector{Float64}, x0::Vector{Float64}; 
        Accel,Schwarzschild=false,Kerr=false,args::Tuple)# key word args 
    n = size(x0)[1]::Int64;
    v = zeros(n,N); v[:,1] = v0::Vector{Float64};# ω
    x = zeros(n,N); x[:,1] = x0::Vector{Float64};# θ

    if Schwarzschild    rs = args[1] end
    if Kerr             a, M = args[1], args[2] end
    for j in 1:N-1
        v0 = v[:,j]
        x0 = x[:,j]
        if (Schwarzschild == true && x0[2] <= rs)
            println("A particle is inside the Schwarzschild radius!")
            println("for loop broken")
            v = v[:,1:j-1]
            x = x[:,1:j-1]#[:,1:j-1]
            break
        end
        if (Kerr == true && ((x0[2] <= M+sqrt(M^2-a^2) + 0.001) || isnan(x0[1]))) 
            println("A particle is inside the event horizon in the Kerr space!")
            println("for loop broken")
            v = v[:,1:j-1]
            x = x[:,1:j-1]#[:,1:j-1]
            break
        end
        ################
        k1 = dt * v0;  # 傾きから変位を求める
        a1 = Accel(v0, x0, args)
        v1 = v0 + (0.5 * dt) * a1; # dt/2.0進んだところの傾き
        x1 = x0 + (0.5) * k1;
        ################
        k2 = dt * v1;
        a2 = Accel(v1, x1, args)#x0 + k1/2.0
        v2 = v0 + (0.5 * dt) * a2; # dt/2.0進んだ所でもう一つ傾きを求める
        x2 = x0 + (0.5) * k2;
        ################
        k3 = dt * v2;
        a3 = Accel(v2, x2, args)#x0 + k2/2.0
        v3 = v0 + dt * a3; # dt進んだところの傾きも求める
        x3 = x0 + k3;
        ################
        k4 = dt * v3;
        a4 = Accel(v3, x3, args)#x0 + k3
        ################
        x[:,j+1]::Vector{Float64} = x0 + (1.0/6.0) * (k1 + 2.0 *k2 + 2.0 *k3 + k4); # より精細にxを近似する
        ################
        v[:,j+1]::Vector{Float64} = v0 + (1.0/6.0 *dt) *(a1 + 2.0 *a2 + 2.0 *a3 + a4);
    end
    #println(v[end])
    return v, x
end
#%%
function newk(i,y0::Vector{Float64},A,k::Matrix{Float64},dt)
    ki = y0
    @simd for j in 1:(i-1)
        @inbounds @views(ki += dt * A[i,j] * (k[:,j]))# matrix 成分とベクトルの積和
    end
    # println(ki)
    return ki
end
function SCprodsum(S::Vector{Float64},M::Matrix{Float64},n::Int)
    _ans = zeros(n)
    Iterations = n #size(M,1)
    for i in 1:Iterations
        for j = eachindex(S)
            @inbounds (_ans[i] += S[j] * M[i,j])# scalar 係数と列ベクトルの和
        end
    end
    return _ans
end
#%
function RK4(N::Int64,dt::Float64,v0::Vector{Float64}, x0::Vector{Float64}; 
        Accel,Schwarzschild=false,Kerr=false,args::Tuple)# key word args 
    n = size(x0)[1]::Int64;# dimension
    v = zeros(n,N); v[:,1] = v0::Vector{Float64};# ω
    x = zeros(n,N); x[:,1] = x0::Vector{Float64};# θ
    # Butcher tableau
    A = ([
        0 0 0 0
        1/3 0 0 0
        -1/3 1 0 0
        1 -1 1 0
    ])#::Matrix{Float64}
    # b = Float64.([1/6, 1/3, 1/3, 1/6])
    b = Float64.([1/8, 3/8, 3/8, 1/8])
    # c = Float64.([0, 1/2, 1/2, 1])
    s = size(A,1) # the number of stages
    k = zeros(n,s)# [zeros(n) for _ in 1:s] #
    a = zeros(n,s)# [zeros(n) for _ in 1:s] #
    X = zeros(n,s)# [zeros(n) for _ in 1:s] #

    if Schwarzschild    rs = args[1] end
    if Kerr             avalue, M = args[1], args[2] end
    for j in 1:N-1
        v0 = v[:,j]
        x0 = x[:,j]
        if (Schwarzschild == true && x0[2] <= rs)
            println("A particle is inside the Schwarzschild radius!")
            println("for loop broken")
            v = v[:,1:j-1]
            x = x[:,1:j-1]#[:,1:j-1]
            break
        end
        if (Kerr == true && ((x0[2] <= M+sqrt(M^2-avalue^2) + 0.001) || isnan(x0[1]))) 
            println("A particle is inside the event horizon in the Kerr space!")
            println("for loop broken")
            v = v[:,1:j-1]
            x = x[:,1:j-1]#[:,1:j-1]
            break
        end
        ################
        @inbounds for i in 1:s
            k[:,i] = newk(i, v0, A ,a, dt)#v1 = v0 + dt * A[2,1] * a[1];
            X[:,i] = newk(i, x0, A ,k, dt)#x1 = x0 + dt * A[2,1] * k[1];
            a[:,i] = Accel(k[:,i], X[:,i], args)#
        end
        ################
        x[:,j+1]::Vector{Float64} = x0 + dt * SCprodsum(b, k, n); # より精細にxを近似する 桁落ち誤差あり
        v[:,j+1]::Vector{Float64} = v0 + dt * SCprodsum(b, a, n);
    end
    #println(v[end])
    return v, x
end

#%%
using LinearAlgebra
function RKF45(N::Int64,dt::Float64,v0::Vector{Float64}, x0::Vector{Float64},
    ϵ::Float64; Accel,Schwarzschild=false,Kerr=false,args::Tuple)# key word args 
    n = size(x0)[1]::Int64;# dimension
    v = zeros(n,N); v[:,1] = v0::Vector{Float64};# ω
    x = zeros(n,N); x[:,1] = x0::Vector{Float64};# θ
    # Butcher tableau
    A21 = 1/4;
    A31 = 3/32; A32 = 9/32;
    A41 = 1932/2197; A42 = -7200/2197; A43 = 7296/2197;
    A51 = 439/216; A52 = -8; A53 = 3680/513; A54 = -845/4104;
    A61 = -8/27; A62 = 2; A63 = -3544/2565; A64 = 1859/4104; A65 = -11/40;
    A = ([
        0.0 0.0 0.0 0.0 0.0 0.0
        A21 0.0 0.0 0.0 0.0 0.0 
        A31 A32 0.0 0.0 0.0 0.0
        A41 A42 A43 0.0 0.0 0.0
        A51 A52 A53 A54 0.0 0.0    
        A61 A62 A63 A64 A65 0.0
    ])#::Matrix{Float64}
    b1 = Float64.([16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55])
    b2 = Float64.([25/216, 0, 1408/2565, 2197/4104, -1/5, 0])
    s = size(A,1) # the number of stages
    k = zeros(n,s)# [zeros(n) for _ in 1:s] #
    a = zeros(n,s)# [zeros(n) for _ in 1:s] #
    X = zeros(n,s)# [zeros(n) for _ in 1:s] #

    if Schwarzschild    rs = args[1] end
    if Kerr             avalue, M = args[1], args[2] end
    for j in 1:N-1
        v0 = v[:,j]
        x0 = x[:,j]
        if (Schwarzschild == true && x0[2] <= rs)
            println("A particle is inside the Schwarzschild radius!")
            println("for loop broken")
            v = v[:,1:j-1]
            x = x[:,1:j-1]#[:,1:j-1]
            break
        end
        if (Kerr == true && ((x0[2] <= M+sqrt(M^2-avalue^2) + 0.001) || isnan(x0[1]))) 
            println("A particle is inside the event horizon in the Kerr space!")
            println("for loop broken")
            v = v[:,1:j-1]
            x = x[:,1:j-1]#[:,1:j-1]
            break
        end
        ################
        @inbounds for i in 1:s
            k[:,i] = newk(i, v0, A ,a, dt)#v1 = v0 + dt * A[2,1] * a[1];
            X[:,i] = newk(i, x0, A ,k, dt)#x1 = x0 + dt * A[2,1] * k[1];
            a[:,i] = Accel(k[:,i], X[:,i], args)#
        end
        ################
        x1 = SCprodsum(b1, k, n); # 5次精度
        x2 = SCprodsum(b2, k, n); # 4次精度
        D = dt * norm(x1 - x2)
        if D > ϵ * dt
            println("The step width is modified.")
            dtprime = dt * (ϵ * dt / D)^(1/4)
            @inbounds for i in 1:s
                k[:,i] = newk(i, v0, A ,a, dtprime)#v1 = v0 + dt * A[2,1] * a[1];
                X[:,i] = newk(i, x0, A ,k, dtprime)#x1 = x0 + dt * A[2,1] * k[1];
                a[:,i] = Accel(k[:,i], X[:,i], args)#
            end
            ################
            x1 = SCprodsum(b1, k, n); # 5次精度
            x2 = SCprodsum(b2, k, n); # 4次精度
        end
        if j%(div(N,10))==0
            println(D)
        end
        ################
        x[:,j+1]::Vector{Float64} = x0 + dt * x1; # より精細にxを近似する 桁落ち誤差あり
        v[:,j+1]::Vector{Float64} = v0 + dt * SCprodsum(b1, a, n);
    end
    #println(v[end])
    return v, x
end
#%%
function RK8(N::Int64,dt::Float64,
    v0::Vector{Float64}, x0::Vector{Float64}; Accel
    ,Schwarzschild=false,Kerr=false,args::Tuple)# key word args 
    n = size(x0)[1]::Int64;# dimension
    v = zeros(n,N); v[:,1] = v0::Vector{Float64};# ω
    x = zeros(n,N); x[:,1] = x0::Vector{Float64};# θ
    # Butcher tableau
    sqrt21 = sqrt(21)
    A21 = 1/2;
    A31 = 1/4; A32 = 1/4;
    A41 = 1/7; A42 = (-7-3*sqrt21)/98; A43 = (21+5*sqrt21)/49;
    ################
    A51 = (11*sqrt21)/84; A52 = 0; A53 = ( 18+4*sqrt21) /63; A54 = ( 21-sqrt21) /252;
    A61 = (5*sqrt21)/48; A62 = 0; A63 = ( 9+sqrt21) /36; A64 = ( -231+14*sqrt21) /360; A65 = ( 63-7*sqrt21) /80;
    A71 = (10*sqrt21)/42; A72 = 0; A73 = ( -432+92*sqrt21) /315; A74 = ( 633-145*sqrt21) /90; A75 = ( -504+115*sqrt21) /70; A76 = ( 63-13*sqrt21) /35;
    ################
    A81 = 1/14; A82 = 0; A83 = 0; A84 = 0; A85 = ( 14-3*sqrt21) /126; A86 = ( 13-3*sqrt21) /63; A87 = 1/9;
    A91 = 1/32; A92 = 0; A93 = 0; A94 = 0; A95 = ( 91-21*sqrt21) /576; A96 = 11/72; A97 = ( -385-75*sqrt21) /1152; A98 = ( 63+13*sqrt21) /128;
    AA1 = 1/14; AA2 = 0; AA3 = 0; AA4 = 0; AA5 = 1/9; AA6 = ( -733-147*sqrt21) /2205; AA7 =( 515+111*sqrt21) /504 ; AA8 = ( -51-11*sqrt21) /56; AA9 = ( 132+28*sqrt21) /245;
    AB1 = 0;    AB2 = 0; AB3 = 0; AB4 = 0; AB5 = ( -42+7*sqrt21) /18; AB6 = ( -18+28*sqrt21) /45; AB7 = ( -273-53*sqrt21) /72; AB8 = ( 301+53*sqrt21) /72; AB9 = ( 28-28*sqrt21) /45; ABA = ( 49-7*sqrt21) /18;
    A = ([
        0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 
        A21 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0  
        A31 A32 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 
        A41 A42 A43 0.0 0.0 0.0 0.0 0.0 0.0 0.0 
        A51 A52 A53 A54 0.0 0.0 0.0 0.0 0.0 0.0     
        A61 A62 A63 A64 A65 0.0 0.0 0.0 0.0 0.0 
        A71 A72 A73 A74 A75 A76 0.0 0.0 0.0 0.0
        A81 A82 A83 A84 A85 A86 A87 0.0 0.0 0.0
        A91 A92 A93 A94 A95 A96 A97 A98 0.0 0.0
        AA1 AA2 AA3 AA4 AA5 AA6 AA7 AA8 AA9 0.0
        AB1 AB2 AB3 AB4 AB5 AB6 AB7 AB8 AB9 ABA
    ])#::Matrix{Float64}
    b1 = Float64.([1/20, 0, 0, 0, 0, 0, 0, 49/180, 16/45, 49/180, 1/20])
    # b2 = Float64.([25/216, 0, 1408/2565, 2197/4104, -1/5, 0])
    s = size(A,1) # the number of stages
    k = zeros(n,s)# [zeros(n) for _ in 1:s] #
    a = zeros(n,s)# [zeros(n) for _ in 1:s] #
    X = zeros(n,s)# [zeros(n) for _ in 1:s] #

    if Schwarzschild    rs = args[1] end
    if Kerr             avalue, M = args[1], args[2] end
    for j in 1:N-1
        v0 = v[:,j]
        x0 = x[:,j]
        if (Schwarzschild == true && x0[2] <= rs)
            println("A particle is inside the Schwarzschild radius!")
            println("for loop broken")
            v = v[:,1:j-1]
            x = x[:,1:j-1]#[:,1:j-1]
            break
        end
        if (Kerr == true && ((x0[2] <= M+sqrt(M^2-avalue^2) + 0.001) || isnan(x0[1]))) 
            println("A particle is inside the event horizon in the Kerr space!")
            println("for loop broken")
            v = v[:,1:j-1]
            x = x[:,1:j-1]#[:,1:j-1]
            break
        end
        ################
        @inbounds for i in 1:s
            k[:,i] = newk(i, v0, A ,a, dt)#v1 = v0 + dt * A[2,1] * a[1];
            X[:,i] = newk(i, x0, A ,k, dt)#x1 = x0 + dt * A[2,1] * k[1];
            a[:,i] = Accel(k[:,i], X[:,i], args)#
        end
        ################
        x1 = SCprodsum(b1, k, n); #
        # x2 = SCprodsum(b2, k, n); #
        # D = dt * norm(x1 - x2)
        # if D > ϵ * dt
        #     println("The step width is modified.")
        #     dtprime = dt * (ϵ * dt / D)^(1/4)
        #     @inbounds for i in 1:s
        #         k[:,i] = newk(i, v0, A ,a, dtprime)#v1 = v0 + dt * A[2,1] * a[1];
        #         X[:,i] = newk(i, x0, A ,k, dtprime)#x1 = x0 + dt * A[2,1] * k[1];
        #         a[:,i] = Accel(k[:,i], X[:,i])#
        #     end
        #     ################
        #     x1 = SCprodsum(b1, k, n); # 5次精度
        #     x2 = SCprodsum(b2, k, n); # 4次精度
        # end
        # if j%(div(N,10))==0
        #     println(D)
        # end
        ################
        x[:,j+1]::Vector{Float64} = x0 + dt * x1; # より精細にxを近似する 桁落ち誤差あり
        v[:,j+1]::Vector{Float64} = v0 + dt * SCprodsum(b1, a, n);
    end
    #println(v[end])
    return v, x
end
if abspath(PROGRAM_FILE) == @__FILE__# Rungekutta4 and RKF45 test code
    v0 = [1.0,-0.5,0.01]; x0 = [0.0,20.0,0.0]; 
    N = Int(1e5); dt = 0.001; ϵ = 0.1;
    T = N*dt
    t = collect(0:dt:T)[1 : end-1]
    ################
    include("Accels.jl")
    ################
    println("N = ", N, "; dt = ", dt)
    println("Rungekutta4")
    @time v1, x1 = Rungekutta(N::Int64,dt::Float64,
        v0::Vector{Float64}, x0::Vector{Float64}, ; Accel=AccelSchwarzschild,Schwarzschild = true)
    println("RKF45")
    @time v2, x2 = RKF45(N::Int64,dt::Float64,
        v0::Vector{Float64}, x0::Vector{Float64}, 
        ϵ::Float64, ; Accel=AccelSchwarzschild,Schwarzschild = true)
        # display(v); println()
    println("RK8")
    @time v3, x3 = RK8(N::Int64,dt::Float64,
        v0::Vector{Float64}, x0::Vector{Float64}, ; Accel=AccelSchwarzschild,Schwarzschild = true)    
    # display(x)
    # println(maximum(v2 - v1))
    # println(maximum(x2 - x1))
    println(maximum(v3 - v1))
    println(maximum(x3 - x1))
end
