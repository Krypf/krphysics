def show_vectors(vectors):
    print("Show vectors")
    for v in vectors:
        print(v.display())
def show_diffs(forms):
    print("Exterior derivatives of forms")
    for v in forms:
        print(diff(v).display())

def commutator_field(basis, i, j):
    x = basis[i].bracket(basis[j])
    return x.components()

def show_commutators(basis):
    n = len(basis)
    print(f"Commutators of a basis {basis}")
    for i in range(n):
        for j in range(n):
            if i < j:
                print((i, j), basis[i].bracket(basis[j]).display())

# another form of functions for commutator algebra

def commutator_symbol(basis, i, j, f_sym = f_sym):
    a = basis[i](basis[j](f_sym)) - basis[j](basis[i](f_sym))
    return a

def structure_constants(basis, i, j, chart=chart, f_sym=f_sym):
    # Compute the commutator vector field
    x = commutator_symbol(basis, i, j)
    x_expr = x.expr()
    coeffs = []
    for k in range(dim):
        coeffs.append(x_expr.coefficient(diff(f_sym.expr(), chart[k])))
    return (coeffs)

def compute_structure_tensor(basis, dim=dim, chart=chart, f_sym=f_sym):
    """
    Computes the full structure constants tensor c[k][i][j],
    where c^k_{ij} satisfies c^k_{ji} = -c^k_{ij} (antisymmetry).

    Parameters:
    - basis: list of vector fields
    - dim: dimension of the space
    - chart: coordinate chart
    - f_sym: symbolic version of the frame
    - basis: basis to be used for commutators

    Returns:
    - c: a 3D list [k][i][j] representing c^k_{ij}
    """
    # Initialize 3D list with zeros
    c = [[[0 for j in range(dim)] for i in range(dim)] for k in range(dim)]

    for i in range(dim):
        for j in range(i+1, dim):  # i < j only, then fill antisymmetric part
            coeffs = structure_constants(basis, i, j)
            for k in range(dim):
                c[k][i][j] = coeffs[k]
                c[k][j][i] = -coeffs[k]  # antisymmetry enforced here

    return c

# end commutator algebra

def pairing_forms_vectors(forms, vectors):
    for i in range(len(forms)):  # Iterate over indices of forms
        for j in range(len(vectors)):  # Iterate over indices of vectors
            print((i, j), forms[i](vectors[j]).display())  # Apply forms[i] to vectors[j] and display the result

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
def lower_connection_coefficient(i, j, metric = None, basis = None, forms = None, gauge = ""):
    
    if metric is None | basis is None | forms is None:
        exit("Wow! Some argument is not selected! Input all three geometrical quantities.")
    else:
        pass
    if gauge == "ortho-normal":
        term1 = 0
        term2 = 0
        term3 = 0
        omega = 0
    return 0

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
dim = 4
def lower_connection_coefficients(i, l, j, g, c, dim = dim):
    """
    Computes the sum over mu of:
    gamma_(i, l, j) = 1/2 * (g_{i mu} * c^mu_{l j} - g_{j mu} * c^mu_{l i} + g_{l mu} * c^mu_{i j})

    Parameters:
    - i, l, j: lower indices (integers)
    - g: a 2D list or matrix representing the tensor g[i][mu]
    - c: a 3D list representing the tensor c[mu][a][b]
    - dim: iterable of values mu runs over (e.g. range(4))

    Returns:
    - The symbolic or numeric result of summing gamma_(i, l, j) over mu
    """
    # orthonormal
    return sum(
        (1/2) * (g[i][mu] * c[mu][l][j] - g[j][mu] * c[mu][l][i] + g[l][mu] * c[mu][i][j])
        for mu in range(dim)
    )

def upper_connection_coefficients(i, l, j, g, c, dim = dim, g_inv = None):
    """
    gamma^{i}_(l, j) = sum of g^{i a} gamma_(a, l, j)
    """
    if g_inv == "eta":
        g_inv = g
    
    return sum(
        g_inv[i][a] * lower_connection_coefficients(a, l, j, g, c, dim = dim)
        for a in range(dim)
    )

def all_upper_coefficients(g, c, dim = dim, g_inv = None, manifold = None):
    # Initialize 3D list with zeros
    gamma = [[[0 for j in range(dim)] for l in range(dim)] for i in range(dim)]

    for i in range(dim):
        for l in range(dim):
            for j in range(dim):
                x = upper_connection_coefficients(i, l, j, g, c, dim = dim, g_inv = g_inv)
                gamma[i][l][j] = manifold.scalar_field(x)

    return gamma

def connection_one_forms(solution, gamma, dim=dim):
    omega = [[0 for j in range(dim)] for i in range(dim)]
    for i in range(dim):
        for j in range(dim):
            omega[i][j] = sum(gamma[i][l][j] * solution.omega_forms()[l] for l in range(dim))
            show(omega[i][j].display())
    return omega

def Riemannian_curvature_component(basis, gamma, c, i, j, k, l, dim=dim):
    # term1: e_k (gamma ^i_{lj}) + gamma ^μ_{lj} gamma ^i_{kμ}
    term1 = basis[k](gamma[i][l][j]) + sum(gamma[mu][l][j] * gamma[i][k][mu] for mu in range(dim))
    
    # term2: - (e_l (gamma ^i_{kj}) + gamma ^μ_{kj} gamma ^i_{lμ})
    term2 = -(basis[l](gamma[i][k][j]) + sum(gamma[mu][k][j] * gamma[i][l][mu] for mu in range(dim)))
    
    # term3: - c ^μ_{kl} gamma ^i_{μj}
    term3 = -sum(c[mu][k][l] * gamma[i][mu][j] for mu in range(dim))
    
    return term1 + term2 + term3

def Riemannian_curvature(basis, gamma, c, dim=dim):
    # Create a 4D zero array
    Riemann = [[[[0 for l in range(dim)] for k in range(dim)] for j in range(dim)] for i in range(dim)]
    for i in range(dim):
        for j in range(dim):
            for k in range(dim):
                for l in range(k+1, dim):  # k < l only, then fill antisymmetric part
                    x = Riemannian_curvature_component(basis, gamma, c, i, j, k, l, dim=dim)
                    Riemann[i][j][k][l] = x
                    Riemann[i][j][l][k] = -x
    return Riemann

def Ricci_tensor_component(Riem, i, j, dim=dim):
    Ric = 0
    for mu in range(dim):
        Ric += Riem[mu][i][mu][j]
    return Ric

def Ricci_tensor(Riem, dim=dim):
    Ric = [[0 for j in range(dim)] for i in range(dim)]
    for i in range(dim):
        for j in range(i, dim):  # i <= j only, then fill antisymmetric part
            Ric[i][j] = Ricci_tensor_component(Riem, i, j, dim=dim)
            Ric[j][i] = Ric[i][j]
    return Ric

def scalar_curvature(Ric, g_inv, dim=dim):
    R = 0
    for i in range(dim):
        for j in range(dim):
            R += g_inv[i][j] * Ric[i][j]
            # show(R.display())
    return R

# print 
def show_two_tensor(tensor, dim, array=True, Mathematica=False, result="result.txt"):
    for i in range(dim):
        for j in range(dim):
            # if i < j:
            if array:
                x = tensor[i][j]
            else:# matrix
                x = tensor[i,j]
            if x != 0:
                show(i, j)
                show(x.display())
                if Mathematica:
                    with open(result, "a") as f:
                        f.write(f"({i}, {j}) ")
                        f.write(str(x.display()))
                        f.write("\n")
                    # print(mathematica_code(x))
    return 0

def show_four_tensor(Riem, dim):
    for i in range(dim):
        for j in range(dim):
            for k in range(dim):
                for l in range(dim):
                    if k < l:
                        x = Riem[i][j][k][l]
                        if x != 0:
                            show(i, j, k, l)
                            show(x.display())