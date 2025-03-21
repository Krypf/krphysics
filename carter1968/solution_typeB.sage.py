

# This file was *autogenerated* from the file carter1968/solution_typeB.sage
from sage.all_cmdline import *   # import sage library

_sage_const_0 = Integer(0); _sage_const_1 = Integer(1); _sage_const_2 = Integer(2); _sage_const_3 = Integer(3); _sage_const_4 = Integer(4)# Function to define the metric tensor
def define_metric(M, Z, Delta_lamda, Delta_mu, l, mu):
    # Define the metric tensor
    g = M.metric()

    # Define the metric components
    g[_sage_const_0 , _sage_const_0 ] = Z / Delta_lamda                # g_lambda_lambda
    g[_sage_const_1 , _sage_const_1 ] = Z / Delta_mu                   # g_mu_mu
    g[_sage_const_2 , _sage_const_2 ] = -Delta_lamda / Z               # g_chi_chi
    g[_sage_const_3 , _sage_const_3 ] = (Z - _sage_const_4  * l**_sage_const_2  * mu**_sage_const_2  * Delta_lamda) * Delta_mu  # g_psi_psi
    g[_sage_const_2 , _sage_const_3 ] = g[_sage_const_3 , _sage_const_2 ] = -_sage_const_2  * l * mu * Delta_lamda / Z  # g_chi_psi = g_psi_chi

    return g

# Example: Call the function with appropriate arguments

M = Manifold(_sage_const_4 , 'M', structure='Lorentzian')
# Define coordinates
C = M.chart(names=('lamda', 'mu', 'chi', 'psi',)); (lamda, mu, chi, psi,) = C._first_ngens(4)

Z, Delta_lamda, Delta_mu, l, mu = var('Z Delta_lamda Delta_mu l mu')
g = define_metric(M, Z, Delta_lamda, Delta_mu, l, mu)
print(g.display())

__time__ = cputime(); __wall__ = walltime(); ginv = g.inverse(); print("Time: CPU {:.2f} s, Wall: {:.2f} s".format(cputime(__time__), walltime(__wall__)))
print(ginv.display())

print("========")
# 変数の定義
var('lamda h m e q p Lamda')

# 数式の定義
l = _sage_const_0 
Z = lamda**_sage_const_2  + l**_sage_const_2 
Lamda = _sage_const_0 
Delta_lamda = Lamda * (_sage_const_1 /_sage_const_3  * lamda**_sage_const_4  + _sage_const_2  * l**_sage_const_2  * lamda**_sage_const_2  - l**_sage_const_4 )               + h * (lamda**_sage_const_2  - l**_sage_const_2 ) - _sage_const_2  * m * lamda + e**_sage_const_2 
Delta_mu = -h * mu**_sage_const_2  + _sage_const_2  * q * mu + p

# g_at_Lamda_zero
g = define_metric(M, Z, Delta_lamda, Delta_mu, l, mu)
print(g.display())

__time__ = cputime(); __wall__ = walltime(); ric = g.ricci(); print("Time: CPU {:.2f} s, Wall: {:.2f} s".format(cputime(__time__), walltime(__wall__)))
print(ric.display())

