load("~/krphysics/sagemath/config.sage")
import time

def commutator_field(basis, i, j):
    x = basis[i].bracket(basis[j])
    return x.components()

def metric_value(i, j, metric = None, vectors = None):
    if metric is None:
        exit("Input the metric tensor.")
    if basis is None:
        exit("Input the basis vectors.")
    return metric(vectors[i], vectors[j])

def show_metric_values(metric, vectors):
    dim = metric.domain().dimension()
    n = len(vectors)
    if dim != n:
        exit(1)
    print("Show all metric values")
    print(metric.display())
    print(vectors.display())
    for i in range(n):
        for j in range(n):
            print((i, j), metric_value(i, j, metric = metric, vectors = vectors).display())

def first_christoffel_symbol(i, j, k, metric = None, basis = None):
    if metric is None:
        exit("Input the metric tensor.")
    if basis is None:
        exit("Input the basis vectors.")
    # Define the Christoffel symbols Γ_{ijk}
    term1 = (basis[i](metric(basis[j], basis[k])))  # e_i(g_{jk})
    term2 = (basis[j](metric(basis[k], basis[i])))  # e_j(g_{ki})
    term3 = (basis[k](metric(basis[i], basis[j])))  # e_k(g_{ij})
    return (term1 + term2 - term3) / 2

def structure_coefficient(mu, i, j, basis = None, forms = None):
    if basis is None:
        exit("Input the basis vectors.")
    if forms is None:
        exit("Input the forms (covariant vectors).")
    com = basis[i].bracket(basis[j])
    return forms[mu](com)

def lower_connection_term(i, L, j, metric, basis, forms, gauge = ""):
    dim = metric.domain().dimension()
    _ans = sum(
        metric_value(i, mu, metric = metric, vectors = basis) * structure_coefficient(mu, L, j, basis = basis, forms = forms)
        for mu in range(dim)
        )
    return _ans

def show_structure_coefficients(basis, forms):
    num_basis = len(basis)
    num_forms = len(forms)
    print("Display all first structure coefficients:")
    for i in range(num_basis):  # Loop over all i, j, mu
        for j in range(num_basis):
            for mu in range(num_forms):
                c_ijk = structure_coefficient(mu, i, j, forms=forms, basis=basis)
                print(f"c^{mu}_{i}{j} = {c_ijk.display()}")
    return 0

def show_first_christoffel_symbols(metric, basis):
    dim = metric.domain().dimension()
    n = len(basis)
    if dim != n:
        print("The length of the basis set is different from the dimension of a manifold.")
        exit(1)
    print("Display all first Christoffel symbols:")
    for i in range(n):  # Loop over all i, j, k
        for j in range(n):
            for k in range(n):
                Gamma_ijk = first_christoffel_symbol(i, j, k, metric=metric, basis=basis)
                print(f"Γ_{i}{j}{k} = {Gamma_ijk.display()}")
    return 0

# form solution.sage
# Levi-Civita connection
def lower_connection_coefficient(g, c, i, l, j, ortho_normal=True):
    """
    Computes the sum over mu of:
    gamma_(i, l, j) = 1/2 * (g_{i mu} * c^mu_{l j} - g_{j mu} * c^mu_{l i} + g_{l mu} * c^mu_{i j})

    Parameters:
    - i, l, j: lower indices (integers)
    - g: a 2D list or matrix representing the tensor g[i][mu]
    - c: a 3D list representing the tensor c[mu][a][b]

    Returns:
    - The symbolic or numeric result of summing gamma_(i, l, j) over mu
    """
    N = len(c)
    x = sum(
        (g[i][mu] * c[mu][l][j] - g[j][mu] * c[mu][l][i] + g[l][mu] * c[mu][i][j])
        for mu in range(N)
        )
    if ortho_normal:
        return (1/2) * x
    else:
        return 0
        # first_christoffel_symbol

def upper_connection_coefficient(g, c, i, l, j, g_inv = None):
    """
    gamma^{i}_(l, j) = sum of g^{i mu} gamma_(mu, l, j)
    """
    if g_inv == "eta":
        g_inv = g
    
    N = len(c)

    return sum(
        g_inv[i][mu] * lower_connection_coefficient(g, c, mu, l, j)
        for mu in range(N)
    )

def all_upper_coefficients(g, c, g_inv = None, domain = None):
    t0 = time.time()
    print("The all_upper_coefficients (Levi-Civita connection) has been calculated... ")

    # Initialize 3D list with zeros
    N = len(c)
    gamma = [[[0 for j in range(N)] for l in range(N)] for i in range(N)]

    for i in range(N):
        for l in range(N):
            for j in range(N):
                x = upper_connection_coefficient(g, c, i, l, j, g_inv = g_inv)
                gamma[i][l][j] = domain.scalar_field(x)
    t1 = time.time(); print(f"in {t1 - t0} seconds");

    return gamma


def Riemannian_curvature_component(basis, gamma, c, i, j, k, l):
    N = len(basis)
    # term1: e_k (gamma ^i_{lj}) + gamma ^μ_{lj} gamma ^i_{kμ}
    term1 = basis[k](gamma[i][l][j]) + sum(gamma[mu][l][j] * gamma[i][k][mu] for mu in range(N))
    
    # term2: - (e_l (gamma ^i_{kj}) + gamma ^μ_{kj} gamma ^i_{lμ})
    term2 = -(basis[l](gamma[i][k][j]) + sum(gamma[mu][k][j] * gamma[i][l][mu] for mu in range(N)))
    
    # term3: - c ^μ_{kl} gamma ^i_{μj}
    term3 = -sum(c[mu][k][l] * gamma[i][mu][j] for mu in range(N))
    
    return term1 + term2 + term3

def Riemannian_curvature(basis, gamma, c):
    t0 = time.time()
    print("The Riemannian_curvature has been calculated... ")
    # Create a 4D zero array
    N = len(basis)
    if N != len(c):
        exit(1)
    Riemann = [[[[0 for l in range(N)] for k in range(N)] for j in range(N)] for i in range(N)]
    for i in range(N):
        for j in range(N):
            for k in range(N):
                for l in range(k+1, N):  # k < l only, then fill antisymmetric part
                    x = Riemannian_curvature_component(basis, gamma, c, i, j, k, l)
                    Riemann[i][j][k][l] = x
                    Riemann[i][j][l][k] = -x
    t1 = time.time(); print(f"in {t1 - t0} seconds");
    return Riemann

def Ricci_tensor_component(Riem, i, j):
    Ric = 0
    N = len(Riem)
    for mu in range(N):
        Ric += Riem[mu][i][mu][j]
    return Ric

def Ricci_tensor(Riem):
    t0 = time.time()
    print("The Ricci_tensor has been calculated... ")
    N = len(Riem)
    Ric = [[0 for j in range(N)] for i in range(N)]
    for i in range(N):
        for j in range(i, N):  # i <= j only, then fill antisymmetric part
            Ric[i][j] = Ricci_tensor_component(Riem, i, j)
            Ric[j][i] = Ric[i][j]
    t1 = time.time(); print(f"in {t1 - t0} seconds");
    return Ric

def scalar_curvature(Ric, g_inv):
    t0 = time.time()
    print("The scalar_curvature has been calculated... ")

    R = 0
    N = len(Ric)
    if N != (g_inv).nrows() | N != (g_inv).nrows():
        exit("Size unmatched!")
    for i in range(N):
        for j in range(N):
            R += g_inv[i][j] * Ric[i][j]

    t1 = time.time(); print(f"in {t1 - t0} seconds");
    return R

def calculate_einstein_tensor(Ric, R, g, cosmological = True):
    # from solution.sage
    # curvatures given
    print("Returns the Einstein tensor...")
    Einstein_tensor = matrix(Ric) - 1 / 2 * R * g
    # show_two_tensor(Einstein_tensor, _dim, array=False)
    if cosmological:
        print("\tminus the cosmological term...")
        LHS = Einstein_tensor - (Lambda * g)
        return LHS
    return Einstein_tensor

# another form of functions for commutator algebra

def commutator_symbol(basis, i, j, f_sym):
    a = basis[i](basis[j](f_sym)) - basis[j](basis[i](f_sym))
    return a

def structure_constant_symbol(basis, i, j, chart, f_sym):
    # Compute the commutator vector field
    x = commutator_symbol(basis, i, j, f_sym)
    x_expr = x.expr()
    coeffs = []
    for k in range(dim):
        coeffs.append(x_expr.coefficient(diff(f_sym.expr(), chart[k])))
    return (coeffs)

def compute_structure_tensor(basis, chart, f_sym):
    """
    Computes the full structure constants tensor c[k][i][j],
    where c^k_{ij} satisfies c^k_{ji} = -c^k_{ij} (antisymmetry).

    Parameters:
    - basis: list of vector fields
    - chart: coordinate chart
    - f_sym: symbolic version of the frame

    Returns:
    - c: a 3D list [k][i][j] representing c^k_{ij}
    """
    dim = len(basis)
    # Initialize 3D list with zeros
    c = [[[0 for j in range(dim)] for i in range(dim)] for k in range(dim)]

    for i in range(dim):
        for j in range(i+1, dim):  # i < j only, then fill antisymmetric part
            coeffs = structure_constant_symbol(basis, i, j, chart, f_sym)
            for k in range(dim):
                c[k][i][j] = coeffs[k]
                c[k][j][i] = -coeffs[k]  # antisymmetry enforced here

    return c
# end commutator algebra

def second_fundamental_form_definition(gamma, normal_form, i, j):
    N = len(normal_form)
    if N != len(gamma):
        exit(1)
    return sum(gamma[alpha][i][j] * normal_form[alpha] for alpha in range(N))

def second_fundamental_form(gamma, normal_form, function = second_fundamental_form_definition):
    N = len(normal_form)
    if N != len(gamma):
        exit(1)
    K = [[0 for j in range(N)] for i in range(N)]
    for i in range(N):
        for j in range(N):
            K[i][j] = function(gamma, normal_form, i, j)
    return K

def trace_free_second(K, g, g_inv):
    N = len(K)
    trace_K = scalar_curvature(K, g_inv)
    return matrix(K) - (1 / N) * trace_K * g