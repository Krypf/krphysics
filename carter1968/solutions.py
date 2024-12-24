#%%
from sympy import symbols, Function, Matrix
from IPython.display import display
from einsteinpy.symbolic import MetricTensor
# , RicciTensor, EinsteinTensor
from sympy.tensor.array import Array
from sympy.matrices.dense import diag
#%%
import sys
import os
absolute_path = os.path.expanduser("~/krphysics")
# 自作モジュールがあるディレクトリを追加
sys.path.append(absolute_path)
from krphysics.ds_ads import *
#%%
# Define constants
from sympy.abc import e, m, n
from sympy.abc import q, p
from sympy.abc import h
from sympy.abc import alpha
Λ = Symbol("Lambda")
# Defining the coordinates as symbols
coordinates = symbols('lambda mu chi psi')
λ, μ, chi, psi = coordinates
mu = μ
# d = Function("d")
#%%
class fields_typeA:
    def metric(matrix_type=True):
        lambda_ = λ
        # 計量の成分を計算する準備
        Z = lambda_**2 + mu**2
        # \Delta_lambda と \Delta_mu の定義
        Delta_lambda = Λ * lambda_**4 / 3 + h * lambda_**2 - 2 * m * lambda_ + p + e**2
        Delta_mu = Λ * mu**4 / 3 - h * mu**2 + 2 * q * mu + p

        g_lambda_lambda = Z / Delta_lambda
        g_mu_mu = Z / Delta_mu

        dchi_dpsi_term = sp.Matrix([
            [Delta_mu, 0],
            [0, -Delta_lambda]
        ])

        mixing = sp.Matrix([
            [-lambda_**2, 1],
            [mu**2, 1]
        ])
        mixed_terms = mixing.T @ dchi_dpsi_term @ mixing

        g_chi_chi = mixed_terms[0, 0] / Z
        g_psi_psi = mixed_terms[1, 1] / Z
        g_chi_psi = mixed_terms[0, 1] / Z

        g_matrix = sp.Matrix([
            [g_lambda_lambda, 0, 0, 0],
            [0, g_mu_mu, 0, 0],
            [0, 0, g_chi_chi, g_chi_psi],
            [0, 0, g_chi_psi, g_psi_psi]
        ])
        return g_matrix

    def field_strength():
        F = zeros(4)# electromagnetic tensor (field strength)
        
        F[0, 3] = -(e*λ**2*μ**2*cos(alpha) - e*μ**4*cos(alpha) - 2*e*λ*μ**3*sin(alpha))/(λ**4 + 2*λ**2*μ**2 + μ**4) 
        F[0, 2] = - (e*λ**2*cos(alpha) - e*μ**2*cos(alpha) - 2*e*λ*μ*sin(alpha))/(λ**4 + 2*λ**2*μ**2 + μ**4)
        F[1, 3] =  + (2*e*λ**3*μ*cos(alpha) + e*λ**4*sin(alpha) - e*λ**2*μ**2*sin(alpha))/(λ**4 + 2*λ**2*μ**2 + μ**4)
        F[1, 2] = - (2*e*λ*μ*cos(alpha) + e*λ**2*sin(alpha) - e*μ**2*sin(alpha))/(λ**4 + 2*λ**2*μ**2 + μ**4)
        for i in range(4):
            for j in range(4):
                if i > j:
                    F[i, j] = - F[j, i]
        return F


#%%

class fields_typeD:
    def metric(matrix_type=True):
        # Define the function Delta_lambda
        Delta_lambda = (Λ + e**2) * λ**2 - 2 * m * λ + n
        # Define the function Delta_mu
        Delta_mu = (Λ - e**2) * mu**2 + 2 * q * mu + p
        metric_matrix = diag([
            1 / Delta_lambda, 1 / Delta_mu, Delta_mu, - Delta_lambda
            ],unpack=True)
        if matrix_type:
            return metric_matrix
        else:
            m_obj = MetricTensor(Array(metric_matrix), coordinates)
            display(m_obj)
            return m_obj

    def field_strength():
        F = zeros(4)# electromagnetic tensor (field strength)
        F[0, 3] = e * cos(alpha) 
        F[3, 0] = - F[0, 3]
        F[1, 2] = - e * sin(alpha) 
        F[2, 1] = - F[1, 2]
        return F
    
#%%
def energy_momentum_tensor(fields):
    metric_matrix = fields.metric()
    F = fields.field_strength()
    g_inv = (metric_matrix).inv()
    trace_term = (F.transpose() * g_inv * F * g_inv).trace() * metric_matrix / 4
    four_pi_T = F * g_inv * F.transpose() - trace_term
    # multiplied by 4 \pi
    return four_pi_T

def main(fields, coordinates):
    metric_matrix = fields.metric()
    
    G = einstein_tensor(metric_matrix, coordinates, output=False)
    LHS = simplify((Matrix(G) - (Λ * metric_matrix)))
    display(LHS)
    four_pi_T = energy_momentum_tensor(fields)
    RHS = simplify(2 * four_pi_T)
    display(LHS - RHS)
    
    return 0
#%%
if __name__ == '__main__':
    # main(fields_typeD, coordinates)
    main(fields_typeA, coordinates)
#%%
