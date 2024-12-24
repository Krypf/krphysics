# 4-dimensional manifold:
M = Manifold(4, 'M')
U = M.open_subset('U')
chart.<λ, μ, ψ, χ> = U.chart() 
eU = chart.frame() 
var('e alpha')  # e と α を定数として扱う
def show_field_strength(A):
    print(A.display(eU)) # 微分形式 a を表示
    print(A.exterior_derivative().display(eU))
# Differential form of degree 1 on M
# Define the four-potential
A_typeA = M.diff_form(1, name='A_A')
A_typeB = M.diff_form(1, name='A_B')
A_typeC = M.diff_form(1, name='A_C')
A_typeD = M.diff_form(1, name='A_D')

# 微分形式の成分を設定
A_typeA[eU, :] = [0, 0,
    e * λ * μ * (μ * cos(alpha) + λ * sin(alpha)) / (λ^2 + μ^2),
    e * (λ * cos(alpha) - μ * sin(alpha)) / (λ^2 + μ^2)]
show_field_strength(A_typeA)

A_typeD[eU, :] = [0, 0, e * λ * cos(alpha), - e * μ * sin(alpha)]  # 成分を指定
# show_field_strength(A_typeD)