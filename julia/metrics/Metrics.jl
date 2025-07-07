using LinearAlgebra
function SchwarzschildMetric(x; eta00=-1, rs = 1, M=0.5)
    r = x[2]; θ = pi / 2;
    f = 1 - rs / r
    # g = Diagonal([-f, 1/f, r^2 ,(r * sin(θ))^2 ])
    g = (- eta00) * Diagonal([-f, 1/f, (r * sin(θ))^2 ])
    return g
end

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

function metric_value(metric_matrix, vector1, vector2)
    n = length(vector1)
    A = SCprodsum(vector2, Matrix(metric_matrix), n)
    g12 = dot(vector1, A)
    return g12
end
function metric_test(vectors, positions, N)
    _ans = zeros(N)
    for j in 1:N
        g = SchwarzschildMetric(positions[:, j])
        v = vectors[:,j]
        a = metric_value(g, v, v)
        _ans[j] = a
    end
    return _ans
end