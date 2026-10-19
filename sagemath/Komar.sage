# Set up the Kerr spacetime by Kerr.sage
load("~/krphysics/sagemath/Kerr.sage")

K = set_const()
g = K.metric_tensor()
g_inv = g.inverse()
r, theta, M, a = K.components()

# 1. 計量からLevi-Civita接続（Christoffel記号）を取得
# 1. get the Levi-Civita connection (Christoffel symbols) from the metric tensor
nabla = g.connection()
Gamma = nabla.coef()

# 2. Komar 積分を計算するため、∇^0 ξ^r_{(0)} の計算
# 2. compute ∇^0 ξ^r_{(0)} for the Komar integral
# ∇^0 ξ^r_{(0)} = g^{00} Γ^r_{00} + g^{0φ} Γ^r_{φ0}
# indices in SageMath: 0:t, 1:r, 2:theta, 3:phi
term1 = g_inv[0,0] * Gamma[1,0,0]
term2 = g_inv[0,3] * Gamma[1,3,0]

# sum
nabla_up0_xi_r = (term1 + term2).expr().simplify_rational()
print("∇^0 ξ^r_{(0)} =", nabla_up0_xi_r)

# 3. 無限遠 r -> ∞ での極限評価
# 3. evaluate the asymptotic behavior at infinity r -> ∞
# r = 1/u と置き換え、u -> 0 の周りでテイラー展開して主要項を取得
# substitute r = 1/u and expand around u = 0 to get the leading term
u = var('u')
asymptotic_term = nabla_up0_xi_r.subs(r == 1/u)
limit_term = asymptotic_term.series(u, 3).truncate()

print("r -> ∞ での漸近形:", limit_term.subs(u == 1/r))

# ------------

# phi component

# 1. compute ∇^0 ξ^r_{(φ)} for the Komar integral
# ∇^0 ξ^r_{(φ)} = g^{00} Γ^r_{0φ} + g^{0φ} Γ^r_{φφ}
# indices in SageMath: 0:t, 1:r, 2:theta, 3:phi
term1_phi = g_inv[0,0] * Gamma[1,0,3]
term2_phi = g_inv[0,3] * Gamma[1,3,3]

# simplify as a rational function
nabla_up0_xi_r_phi = (term1_phi + term2_phi).expr().simplify_rational()
print("∇^0 ξ^r_{(φ)} =", nabla_up0_xi_r_phi)

# 2. evaluate the asymptotic behavior at infinity r -> ∞
# u = var('u')
asymptotic_term_phi = nabla_up0_xi_r_phi.subs(r == 1/u)

# O(1/r^2) の主要項を取り出すため、u=0 の周りでテイラー展開
# Taylor expansion around u = 0 to get the leading term
limit_term_phi = asymptotic_term_phi.series(u, 3).truncate()

print("r -> ∞ での漸近形:", limit_term_phi.subs(u == 1/r))
