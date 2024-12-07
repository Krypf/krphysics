# %%
from sympy import *
from sympy.abc import a,b,c,u,r,t,H
v = Symbol('varphi') # v instead of phi

from IPython.core.display import display

from einsteinpy.symbolic import MetricTensor, ChristoffelSymbols, RiemannCurvatureTensor,RicciTensor, RicciScalar,EinsteinTensor

# %%
def one_sheet_hyperboloid(u,v):
    x = a * cosh(u) * cos(v)
    y = b * cosh(u) * sin(v)
    z = c * sinh(u)
    return Matrix([x,y,z]).subs(b,a).subs(c,a)
def two_sheet_hyperboloid(u,v):
    x = a * sinh(u) * cos(v)
    y = b * sinh(u) * sin(v)
    z = c * cosh(u)
    return Matrix([x,y,z]).subs(b,a).subs(c,a)
# %%
def J(f,args=(u,v)):# jacobian matrix
    x1, x2 = args[0], args[1]
    X = Matrix([x1,x2])
    return f(x1,x2).jacobian(X)

def jacobian_matrix(fc,x,args=None):# jacobian matrix
    if args is None:
        args = x
    return fc(x).jacobian(Matrix(args))

# surface definition end
# %%
def eta(n):# Minkowski metric
    x = eye(n)
    x[-1,-1] = -1
    return x

def flat_metric(args, time_zero = True, _sign = 1):
    r, s = args[0], args[1]
    # the number of plus one = r
    # the number of minus one = s
    n = r + s
    x = eye(n)
    
    for j in range(s):
        if time_zero:
            k = j
        else:
            k = -(j+1)
        
        x[k,k] = - 1
    return x * _sign # sign of the Clifford algebra
# flat_metric(1,2) = diag(-1,-1,1)
# %%

def g1(f):# metric
    j = J(f)
    return simplify((j.T @ j))

def metric_from_jacobian(fc,x,_sign,__time_zero = True, args_jac=None):
    jac = jacobian_matrix(fc,x,args_jac)
    # jac.T is the transpose of j
    return simplify((jac.T @ flat_metric(_sign,time_zero = __time_zero) @ jac))

_sign_of_ads2 = (1,2)
_sign_of_ds2 = (2,1)
# %%
def scalar(metric_matrix, x, output=True):
    metric_matrix = Array(metric_matrix)
    m_obj = MetricTensor(metric_matrix, x)
    # print('metric'); display(m_obj);
    Ric = RicciTensor.from_metric(m_obj)
    R = RicciScalar.from_riccitensor(Ric)
    
    if output == True:
        print("Ricci tensor")
        display(Ric.tensor().simplify())
        
        print("Ricci scalar")
        display(R.simplify())
    
    return Ric, R

def einstein_tensor(metric, x, output=True):
    ################
    Ric, R = (scalar(metric, x, output=True))
    Ric, R = Ric.simplify(), R.simplify()
    ################
    G = Ric - tensorproduct(R / 2,metric)
    G = G.simplify()
    if output == True:
        print("Einstein tensor")
        display(G)
    
    return (G)
# %%
import sympy as sp

def get_spherical_vars(n,angles_symbol='\\theta', subscript=0):
    # Define symbolic variables
    # angle
    if n == 1:
        return []
    if subscript == 0:
        spherical_vars = sp.symbols(angles_symbol + '_{0:%d}' % (n-1))
    elif subscript == 1:
        spherical_vars = sp.symbols(angles_symbol + '_{1:%d}' % (n))
    return spherical_vars

def spherical_to_cartesian(n,radius = sp.Symbol("r"), _angles_symbol = '\\theta', A = 1):
    # Define symbolic variables
    # angle
    if n == 1:
        return sp.Matrix([radius])
    spherical_vars = get_spherical_vars(n,angles_symbol=_angles_symbol)
        
    # Create the coordinate mappings
    cartesian_coords = []
    polar_angle = spherical_vars[0] / A # theta_0
    cartesian_coords.append(radius * sp.cos(polar_angle))
    product_sin = sp.sin(polar_angle)
    
    for j in range(1, n - 1):
        cartesian_coords.append(radius * sp.cos(spherical_vars[j] / A))
        cartesian_coords[j] *= product_sin
        product_sin *= sp.sin(spherical_vars[j] / A)
    
    cartesian_coords.append(radius * product_sin)

    print("Cartesian coordinates:")
    print(cartesian_coords)
    
    return sp.Matrix(cartesian_coords)

# %%
def spherical_to_de_Sitter(n,radius = sp.Symbol("a"), A = 1, __time_zero = True):
    # Define symbolic variables
    cartesian_coords_sphere = spherical_to_cartesian(n,radius=radius, A = A)
    u = sp.Symbol("u")
    u_overA = u / A
    cartesian_coords_sphere *= sp.cosh(u_overA)
    cartesian_coords_sphere = list(cartesian_coords_sphere)
    # add the last coordinate
    if __time_zero:
        cartesian_coords_sphere.insert(0, radius * sp.sinh(u_overA))
    else:
        cartesian_coords_sphere.append(radius * sp.sinh(u_overA))
    return sp.Matrix(cartesian_coords_sphere)

# %%
def get_cartesian_vars(n):
    return symbols('x1:'+str(n+1))
def cartesian_to_metric_sphere(n):
    x = get_cartesian_vars(n)
    # x sub n + 1
    x_np1_squared = sum([x[j]**2 for j in range(n)])
    return ([[KroneckerDelta(j,k) + x[j] * x[k] / (1 - x_np1_squared) for j in range(n)] for k in range(n)])
def cartesian_to_metric_hyperbola(n):
    x = get_cartesian_vars(n)
    # x sub n + 1
    x_np1_squared = sum([x[j]**2 for j in range(n)])
    return ([[KroneckerDelta(j,k) - x[j] * x[k] / (x_np1_squared - 1) for j in range(n)] for k in range(n)]) 
# %%
from itertools import product

def display_geodesic(matrix_metric, x):
    m_obj = MetricTensor(Array(matrix_metric), x)
    Gamma = ChristoffelSymbols.from_metric(m_obj).simplify()
    
    # variables
    # X = Function("X")(*x)
    s = Symbol("s")
    n = len(x)
    X = [Function("\gamma^"+str(j))(s) for j in range(n)]
    
    L = range(n)
    for mu in L:
        print("The geodesic equation with respect to", end=" ")
        display(x[mu])
        term1 = X[mu].diff(s).diff(s)
        for j,k in product(L,L):
            term1 += Gamma[mu,j,k] * X[j].diff(s) * X[k].diff(s)
        display(term1)
        print(latex(term1))

# %% coordinate systems
def ds_global(n, A = 1, __time_zero = True):
    return spherical_to_de_Sitter(n, A = A, __time_zero = __time_zero)

def ds_global_a(n, __time_zero = True):
    return spherical_to_de_Sitter(n, A = a, __time_zero = __time_zero)

def ds_static(n, A = a):
    # thetas = get_spherical_vars(n - 2)
    r = sp.Symbol("r")
    x_0 = sqrt(a**2 - r**2) * sinh(u / A)
    x_1 = sqrt(a**2 - r**2) * cosh(u / A)
    if n >= 3:
        x_2_n = spherical_to_cartesian((n - 2) + 1,radius = r)
        return Matrix([x_0, x_1] + list(x_2_n))
    else:
        return Matrix([x_0, x_1])
def ds_static_outside(n, A = a):
    # thetas = get_spherical_vars(n - 2)
    r = sp.Symbol("r")
    x_0 = sqrt(r**2 - a**2) * sinh(u / A)
    x_1 = sqrt(r**2 - a**2) * cosh(u / A)
    if n >= 3:
        x_2_n = spherical_to_cartesian((n - 2) + 1,radius = r)
        return Matrix([x_0, x_1] + list(x_2_n))
    else:
        return Matrix([x_0, x_1])

def ds_FLRW(n):
    # thetas = get_spherical_vars(n - 2)
    r = sp.Symbol("r")
    t = sp.Symbol("t")
    t /= a
    x_0 = a * (sp.sinh(t) + Rational(1 / 2) * r**2 * sp.exp(t))
    x_1 = a * (sp.cosh(t) - Rational(1 / 2) * r**2 * sp.exp(t))
    if n >= 3:
        x_2_n = spherical_to_cartesian((n - 2) + 1,radius = a * r * sp.exp(t))
        return Matrix([x_0, x_1] + list(x_2_n))
    else:
        return Matrix([x_0, x_1])

# from itertools import product

def laplacian_mfd(metric_matrix, x, f, D = None):
    _sum = 0
    # parameters
    n = len(x)
    indices = range(n)
    # metric
    g = metric_matrix
    # D = sqrt(abs(det(g)))

    for (j, k) in product(indices, indices):
        # print(j,k)
        factor1 = D * (g.inv())[j,k] * f.diff(x[k])
        _sum += factor1.diff(x[j]) / D
    return _sum
# %% poincare
def ds_metric_poincare(n):
    H = Symbol("H")
    args_flat = (1,n-1)
    
    return flat_metric(args_flat, time_zero=False) / (H * u)**2
def ds_metric_poincare_plus(n):
    # almost plus
    H = Symbol("H")
    args_flat = (n-1,1)
    
    return flat_metric(args_flat, time_zero=True) / (H * u)**2
def ds_tetrad_poincare(n):
    H = Symbol("H")
    
    return eye(n) / (H * u)
# %%
def get_christoffel(metric_matrix, x):
    m_obj = MetricTensor(Array(metric_matrix), x)
    Gamma = ChristoffelSymbols.from_metric(m_obj).simplify()
    return Gamma

################
# %% gamma matrices
from sympy.physics.quantum import TensorProduct

def pauli_matrices(mu: int):
    if   mu == 0:
        sigma = eye(2)# sigma0
    elif mu == 1:
        sigma = Matrix(([0,1],[1,0]))# sigma1
    elif mu == 2:
        sigma = Matrix(([0,-1j],[1j,0]))# sigma2
    elif mu == 3:
        sigma = Matrix(([1,0],[0,-1]))# sigma3
    
    return sigma

def gamma_matrices_dirac(mu):
    # Return gamma matrices in the Dirac representation
    n = 4
    if   mu == 0:
        gamma = TensorProduct(Matrix([[1,0],[0,-1]]), pauli_matrices(0))# gamma^0
    elif (mu == 1) or (mu == 2) or (mu == 3):
        gamma = TensorProduct(Matrix([[0,1],[-1,0]]) , pauli_matrices(mu))# gamma^mu
    elif mu == 5:
        gamma = 1j * prod([gamma_matrices_dirac(mu) for mu in range(n)])# gamma^5
    
    return gamma

def gamma_local_dirac(mu, tetrad):
    n = n_global
    tetrad_inv = tetrad**-1
    _ans = zeros(n,n)
    for a in range(n):
        _ans += gamma_matrices(a, rep="Dirac") * tetrad_inv[mu, a]
    return _ans

# %% spin structure in a curved spacetime

from typing import NamedTuple

class SpinStructure(NamedTuple):
    import sympy
    tetrad: sympy.matrices.dense.MutableDenseMatrix

    metric_matrix: sympy.matrices.dense.MutableDenseMatrix

    args: list# arguments x
    
    sign: tuple# the sign of the Clifford
    
    def dim(self) -> int:
        return len(self.args)
    # dim: int# dimension n -> auto
    
    # rep # gamma matrix -> function is not type

def gamma_local(_spin_structure, gamma_rep, mu):
    tetrad, n = _spin_structure.tetrad, _spin_structure.dim()
    
    tetrad_inv = tetrad**-1
    # display(tetrad_inv)
    _ans = zeros(n,n)
    for a in range(n):
        # pass
        _ans += gamma_rep(a) * tetrad_inv[mu, a]
    return _ans 

from clifford import anti_commutator

def confirm_gamma(_spin_structure, gamma_rep):
    tetrad, n = _spin_structure.tetrad, _spin_structure.dim()

    indices = range(n)
    g_inv = g**-1
    # local gamma commutation relation
    for mu, nu in product(indices, indices):
        if mu >= nu:
            print(mu, nu)
            x = anti_commutator(gamma_local(_spin_structure, gamma_rep, mu), gamma_local(_spin_structure, gamma_rep, nu))
            display(x)
            if mu == nu:
                display(x - 2 * g_inv[mu, nu] * eye(n))

# %%
def spin_connection_from_tetrad(_spin_structure, a, b, mu):
    tetrad, g, x = _spin_structure.tetrad, _spin_structure.metric_matrix, _spin_structure.args
    
    n = _spin_structure.dim()
    
    g_inv = g**-1 # auto-calculation when g is diagonal
    e_super = tetrad @ g_inv# e^{a}{}_{\mu} g^{\mu\nu}
    _sign = _spin_structure.sign
    e_sub = flat_metric(_sign, time_zero=False) @ tetrad
    _const = 1 / 2

    term1 = 0
    for nu in range(n):
        term1 += e_super[a, nu] * (tetrad[b, nu].diff(x[mu]) - tetrad[b, mu].diff(x[nu])) * _const
    term2 = 0
    for nu in range(n):
        term2 += e_super[b, nu] * (tetrad[a, nu].diff(x[mu]) - tetrad[a, mu].diff(x[nu])) * _const
    term3 = 0
    for rho in range(n):
        for sigma in range(n):
            for c in range(n):
                term3 += e_super[a, rho] * e_super[b, sigma] * (e_sub[c, sigma].diff(x[rho]) - e_sub[c, rho].diff(x[sigma])) * tetrad[c, mu] * _const
    
    return term1 - term2 - term3

# %%
def spin_connection_mu(_spin_structure, mu):
    n = _spin_structure.dim()
    _ans = zeros(n, n)
    indices = range(n)
    
    for a, b in product(indices, indices):
        _ans[a, b] = spin_connection_from_tetrad(_spin_structure, a, b, mu)
    return _ans

from clifford import commutator

def sigma_matrices_general(gamma_rep, mu, nu):
    _const = 1j / 4
    gamma_mu = gamma_rep(mu)
    gamma_nu = gamma_rep(nu)
    return _const * commutator(gamma_mu, gamma_nu)

def omega_mu(_spin_structure, gamma_rep, mu):
    n = _spin_structure.dim()
    _ans = zeros(n, n)
    indices = range(n)
    _const = 1j / 2
    # tetrad = _spin_structure.tetrad
    # metric_matrix = _spin_structure.metric_matrix
    # x = _spin_structure.args
    
    for a, b in product(indices, indices):
        # _ans =  spin_connection_components(metric_matrix, x, tetrad, a, b, mu) * sigma_matrices(a, b)
        _ans += spin_connection_from_tetrad(_spin_structure, a, b, mu) * sigma_matrices_general(gamma_rep, a, b)
        # ! not = but +=
    return _ans * _const
# %%

def diff_spin_connection(spinor, _spin_structure, gamma_rep, mu):
    x = _spin_structure.args
    return spinor.diff(x[mu]) + omega_mu(_spin_structure, gamma_rep, mu) * spinor

def dirac_operator(spinor, _spin_structure, gamma_rep, _const=1):
    # _const = 1 or 1j
    n = _spin_structure.dim()
    x = _spin_structure.args
    tetrad = _spin_structure.tetrad
    
    _ans = zeros(n,1)
    for mu in range(n):
        _ans += simplify(_const * gamma_local(_spin_structure, gamma_rep, mu) * diff_spin_connection(spinor, _spin_structure, gamma_rep, mu))
    return simplify(_ans)


if __name__ == '__main__':
    n = 4
    tetrad = e = ds_tetrad_poincare(n)
    metric_matrix = g = ds_metric_poincare(n)
    args = x = ([u] + list(get_cartesian_vars(n-1)))
    # rep = γ = gamma_matrices_dirac
    _sign = (1,3)
    _spin_structure = SpinStructure(e, g, x, _sign)
    gamma_rep = gamma_matrices_dirac
    # display(_spin_structure)
    # print(_spin_structure.dim())

    # display(gamma_local(_spin_structure, gamma_rep, 3))
    # confirm_gamma(_spin_structure, gamma_rep)
    
    # a, b, mu = 0, 1, 1; display(spin_connection_from_tetrad(_spin_structure, a, b, mu))
    # display(spin_connection_mu(_spin_structure, 1))
    # display(omega_mu(_spin_structure, gamma_rep, 1))
    
    spinor = Matrix([Function("\psi_"+str(j))(*x) for j in range(n)])
    # display(diff_spin_connection(spinor, _spin_structure, gamma_rep, mu))
    mu = 3
    
    spinor_diff = dirac_operator(spinor, _spin_structure, gamma_rep)
    spinor_diff2 = dirac_operator(spinor_diff, _spin_structure, gamma_rep)

    display(spinor)
    display(factor(spinor_diff))
    display(factor(spinor_diff2))
# %%
# n = 3; _sign = (n, 1);
# n = 4; _sign = (n, 1); # ds4
# n = 5; _sign = (n, 1);
if __name__ == '__main__':
    args_jac = [u,r] + list(get_spherical_vars(n - 2 + 1) )
    g = metric_from_jacobian(ds_static_outside, n,_sign, __time_zero = True, args_jac=args_jac);
    display(g)

# %%
