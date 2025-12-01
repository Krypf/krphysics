load("spacetime.sage")
from itertools import product

def eps_oneup(metric, g_inv):
    U = metric.domain()
    n = U.dimension()
    eps_above = U.tensor_field(1, 3, name='ε^i_j1j2j3')
    eps = metric.volume_form()
    for k in range(n):
        for i, j1, j2, j3 in product(range(n), repeat=4):
            # (k, j1, j2, j3) is a permutation
            s = 0
            for k in range(n):
                s += g_inv[i,k] * eps[k, j1, j2, j3]
            eps_above[i, j1, j2, j3] = simplify(s)
    return eps_above

def hodge_one_form(one_form, epsilon_tensor):
    z = epsilon_tensor.contract(0, one_form)
    return z

def eps_twoup(metric, g_inv):
    """
    Build ε^{i1 i2}{}_{j1 j2} = g^{i1 k1} g^{i2 k2} ε_{k1 k2 j1 j2}
    for a 4-dimensional metric (Lorentzian or Euclidean).

    INPUT:
        metric : SageManifolds metric tensor g_ij
        g_inv  : inverse metric tensor g^{ij}

    OUTPUT:
        eps_above : (2,2)-type tensor ε^{i1 i2}{}_{j1 j2}
    """
    U = metric.domain()
    n = U.dimension()
    eps_cov = metric.volume_form()   # ε_{i1 i2 i3 i4}
    eps_above = U.tensor_field(2, 2, name='ε^i1i2_j1j2')

    for i1, i2, j1, j2 in product(range(n), repeat=4):
        s = 0
        for k1 in range(n):
            for k2 in range(n):
                s += g_inv[i1, k1]*g_inv[i2, k2]*eps_cov[k1, k2, j1, j2]
        eps_above[i1, i2, j1, j2] = simplify(s)
    return eps_above

def hodge_two_form(two_form, epsilon_tensor):
    z = epsilon_tensor.contract(0, 1, two_form, 0, 1) / 2
    return z

# Minkowski space = Min
dd = Min.chart.coframe()
d0 = dd[0] # abbreviation
d1 = dd[1]
d2 = dd[2]
d3 = dd[3]
pp = Min.chart.frame()
p0 = pp[0]
metric, g_inv = eta, eta
y = eps_oneup(metric, g_inv)
result_one_form = {
    0: hodge_one_form(d0, y) == - d1.wedge(d2).wedge(d3),
    1: hodge_one_form(d1, y) == - d0.wedge(d2).wedge(d3),
    2: hodge_one_form(d2, y) == - d0.wedge(d3).wedge(d1),
    3: hodge_one_form(d3, y) == - d0.wedge(d1).wedge(d2)
    }
if all(result_one_form.values()):
    print("The results of the hodge star of one (1) forms are all true")

y2 = eps_twoup(metric, g_inv)
result_two_form = {
    (0, 1): hodge_two_form(dd[0].wedge(dd[1]), y2) == - (d2).wedge(d3),
    (0, 2): hodge_two_form(dd[0].wedge(dd[2]), y2) == - (d3).wedge(d1),
    (0, 3): hodge_two_form(dd[0].wedge(dd[3]), y2) == - (d1).wedge(d2),
    (1, 2): hodge_two_form(dd[1].wedge(dd[2]), y2) == (d0).wedge(d3),
    (2, 3): hodge_two_form(dd[2].wedge(dd[3]), y2) == (d0).wedge(d1),
    (3, 1): hodge_two_form(dd[3].wedge(dd[1]), y2) == (d0).wedge(d2),
}
if all(result_two_form.values()):
    print("The results of the hodge star of two (2) forms are all true")

# mistaken program
# a1 = d0.hodge_dual(eta, minus_eigenvalues_convention=True)
# a2 = d0.hodge_dual(eta)
# a3 = d1.hodge_dual(eta)
# print(a1.display())
# print(a2.display())
# print(a3.display())

# def hodge_dual_one_form(one_form, frame, metric):
#     eps = metric.volume_form()
#     p = one_form.degree()
#     n = metric.domain().dimension()
#     coeffs = one_form[:]
#     x = coeffs[0] * eps.contract(0, frame[0], 0)
#     for i in range(1, n):
#         x += coeffs[i] * eps.contract(0, frame[i], 0)
#     return x#.down(metric)

# _ans = hodge_dual_one_form(dd[0], pp, eta); print(_ans.display()) # false
# _ans = hodge_dual_one_form(dd[1], pp, eta); print(_ans.display())
# _ans = hodge_dual_one_form(dd[2], pp, eta); print(_ans.display())
# _ans = hodge_dual_one_form(dd[3], pp, eta); print(_ans.display())

