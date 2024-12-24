#%%
from sympy import symbols, cos, sin
from sympy.diffgeom import Manifold, Patch, CoordSystem, Differential

# 座標系の定義
M = Manifold('M', 4)  # 4次元空間
patch = Patch('P', M)
R4 = CoordSystem('R4', patch, ['λ', 'μ', 'ψ', 'χ'])

# 座標関数と微分形式を取得
λ, μ, ψ, χ = R4.coord_functions()
dλ, dμ, dψ, dχ = R4.base_oneforms()

# ポテンシャル A の定義
e = symbols('e')  # 定数 e
α = symbols('α')  # パラメータ α
f1 = (λ * μ * (μ  * cos(α) + λ * sin(α))) / (λ**2 + μ**2)
f2 = (λ * cos(α) - μ * sin(α)) / (λ**2 + μ**2)

# metric A
A = e * (f1 * dψ + f2 * dχ)

# 外微分 dA を計算
dA = Differential(A)
print(dA)

#%%

#%%
