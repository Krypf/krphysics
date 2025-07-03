#%% modified for loop condition
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
#%%
j = 1
n = 3; N = 2;
x = zeros(n,N)
x[:,1] = [1,1.1,0]
x[:,2] = [1,1.1,3]


using LinearAlgebra
function RKF45plus(N::Int64,dt::Float64,v0::Vector{Float64}, x0::Vector{Float64},
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
        # Kerr, for loop modified!
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
