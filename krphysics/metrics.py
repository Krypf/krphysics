# %%
from sympy.core.symbol import Symbol
from sympy import symbols,sin,cos,latex,sqrt, Matrix
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
    # Sigma = Symbol('Sigma')
    Sigma = r**2 + a**2 * cos(θ)**2
    # Δ = Symbol('Delta')
    Δ = r**2 - 2 * M * r + a**2 
    # ρ = sqrt(Sigma) # power = half
    S = (r**2 + a**2 ) **2 - a**2 * Δ * sin(θ) **2 
    f0 = - (1 - 2 * M * r / Sigma )
    metric = (-eta00) * diag([f0, Sigma / Δ , Sigma, S * sin(θ)**2 / Sigma ],unpack=True)
    metric[0,3] = eta00 * 4*M*a*r*sin(θ)**2 / Sigma / 2# the sign equals eta00. not minus one 
    metric[3,0] = metric[0,3]
    metric = Array(metric)
    m_obj = MetricTensor(metric, x)
    return m_obj
# %%
from sympy.abc import lamda
from IPython.display import display

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
def block_determinant():
    g = KerrMetric()
    print("The KerrMetric is a stationary and axisymmetric solution to the Einstein field equations, describing a rotating black hole. The metric is given by:")
    display(Matrix(g.tensor()))
    print("The determinant of a block matrix:")
    display((g[0, 0] * g[3, 3] - g[0, 3]**2).simplify())

block_determinant()
# %%
# create_result() is a method of Result class, which is used to create a result file for the current script. It is used to store the results of the script in a structured way. The result file is created in the same directory as the script, and can also be used to document the code and explain the results.
# %%
# execute the script as a module
# python3 -m krphysics.metrics 
if __name__ == '__main__':
    from utils.result import Result # print_md, print_latex
    # Initialize it
    res = Result(script_name="krphysics/metrics.py")
    res.create_result()
    g = SchwarzschildMetric(eta00=-1,rs=True)
    res.print_md("Schwarzschild metric"); display(g);
    res.print_latex(g.tensor())

    g = KerrMetric()
    res.print_md("Kerr metric"); display(g);
    res.print_latex(g.tensor())
# %%
# %%

