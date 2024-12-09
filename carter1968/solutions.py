#%%
from sympy import symbols, Function, Matrix
from IPython.display import display
from einsteinpy.symbolic import MetricTensor, RicciTensor, EinsteinTensor
from sympy.tensor.array import Array
from sympy.matrices.dense import diag
#%%
# Defining the coordinates as symbols
coordinates = symbols('lambda mu chi psi')
λ, μ, chi, psi = coordinates
mu = μ
d = Function("d")
#%%
# Define constants and variables
from sympy.abc import e, m, n
from sympy.abc import q, p
Λ = symbols("Lambda")
# Define the function Delta_lambda
Delta_lambda = (Λ + e**2) * λ**2 - 2 * m * λ + n
# Define the function Delta_mu
Delta_mu = (Λ - e**2) * mu**2 + 2 * q * mu + p
#%%
metric_matrix = diag([
    1 / Delta_lambda, 1 / Delta_mu, Delta_mu, - Delta_lambda
    ],unpack=True)

m_obj = MetricTensor(Array(metric_matrix), coordinates)
display(m_obj)
#%%
# Ric = RicciTensor.from_metric(m_obj)
#%%
# G = EinsteinTensor.from_metric(m_obj)
# G = G.simplify()
# for j in range(4):
#     display(G[j])
#%%
import sys
import os
absolute_path = os.path.expanduser("~/krphysics")

# 自作モジュールがあるディレクトリを追加
sys.path.append(absolute_path)
from krphysics.ds_ads import *
# os.chdir("~/krphysics")
#%%
G = einstein_tensor(metric_matrix, coordinates, output=False)
#%%
from sympy.abc import alpha
F = zeros(4)
F[0, 3] = e * cos(alpha) 
F[3, 0] = - F[0, 3]
F[1, 2] = - e * sin(alpha) 
F[2, 1] = - F[1, 2]
#%%
g_inv = (metric_matrix).inv()
trace_term = (F.transpose() * g_inv * F * g_inv).trace() * metric_matrix / 4
four_pi_T = F * g_inv * F.transpose() - trace_term

#%%
LHS = simplify((Matrix(G) - (Λ * metric_matrix)))
LHS
#%%
RHS = simplify(2 * four_pi_T)
#%%
LHS - RHS