# %%
from sympy.core.symbol import Symbol
from sympy import symbols,sin,cos,latex,sqrt
from einsteinpy.symbolic import MetricTensor
from sympy.matrices.dense import diag
from sympy.tensor.array import Array

# %%
def SchwarzschildMetric(M=Symbol('M'),eta00=-1,rs = False):
    # c = G = 1
    x = symbols('t r theta phi')# x0 = t
    r = x[1]
    f = (1 - 2*M / r)# M is mass
    if rs is True:
        rs = Symbol('r_s')
        f = (1 - rs / r)# rs is the Schwarzschild radius
    metric = (-eta00) * diag([-f, 1/f, r**2, (r * sin(x[2]))**2],unpack=True)
    metric = Array(metric)
    m_obj = MetricTensor(metric, x)
    return m_obj
# %%
def KerrMetric(M=Symbol('M'),a=Symbol('a'),eta00=-1):
    x = symbols('t r theta phi')# x0 = t
    r = x[1]; θ = x[2];
    ρ = sqrt(r**2 + a**2 * cos(θ)**2) # power 注意
    Δ = r**2 - 2 * M * r + a**2 
    S = (r**2 + a**2 ) **2 - a**2 * Δ * sin(θ) **2 
    f0 = - (1 - 2 * M * r / ρ**2 )
    metric = (-eta00) * diag([f0, ρ**2 / Δ , ρ**2, S * sin(θ)**2 / ρ**2 ],unpack=True)
    metric[0,3] = eta00 * 4*M*a*r*sin(θ)**2 / ρ**2 / 2# the sign equals eta00. not minus one 
    metric[3,0] = metric[0,3]
    metric = Array(metric)
    m_obj = MetricTensor(metric, x)
    return m_obj
# %%
from sympy.abc import lamda
from IPython.core.display import display

def geodesics(μ, x, Γ, λ = lamda):
    M = [j for j in range(len(x))]
    term2 = 0
    for a in M:
        for b in M:
            term2 += Γ[μ][a][b] * x[a].diff(λ) * x[b].diff(λ)
    # eq = x[μ].diff(λ).diff(λ) + term2.simplify()
    LHS = x[μ].diff(λ).diff(λ)
    display(LHS)
    print(latex(LHS))
    print('=')
    return (-term2)

# %%
if __name__ == '__main__':
    g = SchwarzschildMetric(eta00=-1,rs=True)
    print("Schwarzschild metric"); display(g);
    g = KerrMetric()
    print("Kerr metric"); display(g);
# %%
