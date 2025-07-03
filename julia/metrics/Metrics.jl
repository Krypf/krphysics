using LinearAlgebra
function SchwarzschildMetric(x;eta00=-1,rs = 1,M=0.5)
    r = x[2]; θ = pi / 2;
    f = 1 - rs / r
    # g = Diagonal([-f, 1/f, r^2 ,(r * sin(θ))^2 ])
    g = Diagonal([-f, 1/f, (r * sin(θ))^2 ])
end
# SchwarzschildMetric.([x0,2x0])