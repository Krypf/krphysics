class Schwarzschild:
    # Constants
    G = var('G')    # Gravitational constant
    Mass = var('M_S')    # Mass of the central object
    c = var('c')    # Speed of light
    dim = 4

    def __init__(self,
        args_name = 't r theta phi',
        manifold = Manifold(dim, 'Schwarzschild', structure='Lorentzian'),
        has_radius = True,
        neighborhood_name = 'U',
        ):
        if has_radius:
            self.radius = var('R_S')  # Schwarzschild radius
        else:
            self.radius = 2 * G * Mass / c**2  # Schwarzschild radius
        # Define the manifold
        self.manifold = manifold
        self.args_name = args_name
        self.args = var(args_name)# = t, r, theta, phi 
        self.neighborhood = manifold.open_subset(neighborhood_name) # = U

    def spherical_chart(self):
        # Define the coordinate chart
        C_spherical = self.neighborhood.chart(self.args_name)
        return (C_spherical)

    def coordinate_frame(self):
        return self.spherical_chart().frame()

    def dual_frame(self):
        return self.spherical_chart().coframe()
        
    def potential(self, number_radius = 1):
        r = self.args[number_radius]
        abs_g_00 = 1 - self.radius / r
        return abs_g_00

    def global_domain(self):
        U = self.spherical_chart().domain()# Open subset U of the 4-dimensional Lorentzian manifold Schwarzschild
        return U

    # Metric components
    def components(self, number_radius = 1, number_theta = 2):
        r = self.args[number_radius]
        theta = self.args[number_theta]
        F = self.potential()
        return r, theta, F

    def metric_tensor(self):
        r, theta, F = self.components()
        U = self.global_domain()
        g_S = U.metric('g_S')
        # Define the metric tensor
        g_S[0, 0] = -F
        # Diagonal components
        g_S[1, 1] = 1 / F
        g_S[2, 2] = r^2
        g_S[3, 3] = (r * sin(theta))^2
        # print(g_S.display())
        return g_S

    def ortho_normal_frame(self):
        r, theta, F = self.components()
        partials = self.coordinate_frame()
        e0 = partials[0] / sqrt(F)
        e1 = partials[1] * sqrt(F)
        e2 = partials[2] / r
        e3 = partials[3] / (r * sin(theta))
        es = [e0, e1, e2, e3]
        return es

    def ortho_normal_tetrads(self):
        U = self.global_domain()
        r, theta, F = self.components()

        f0 = U.one_form(sqrt(F), 0, 0, 0, name='f0')
        f1 = U.one_form(0, 1 / sqrt(F), 0, 0, name='f1')
        f2 = U.one_form(0, 0, r, 0, name='f2')
        f3 = U.one_form(0, 0, 0, r * sin(theta), name='f3')
        fs = [f0, f1, f2, f3]
        return fs

    def E_hat(self):
        r, theta, F = self.components()
        # Construct the matrix
        E = matrix([
            [sqrt(F), 0, 0, 0],
            [0, 1 / sqrt(F), 0, 0],
            [0, 0, r, 0],
            [0, 0, 0, r * sin(theta)]
        ])
        return E

    def structure_tetrad(self, i, j):
        basis = self.ortho_normal_frame()
        com = commutator_field(basis, i, j)
        E = self.E_hat()
        n = len(basis)
        cs = [0 for _ in range(n)]
        for k in range(n):
            cs[k] = sum(com[a] * E[k][a] for a in range(n))
        return cs

    def compute_structure_coefficients(self, dim=dim):
        """
        Computes the full structure constants tensor c[k][i][j],
        where c^k_{ij} satisfies c^k_{ji} = -c^k_{ij} (antisymmetry).

        Parameters:
        - dim: dimension of the space
        - basis: basis to be used for commutators

        Returns:
        - c: a 3D list [k][i][j] representing c^k_{ij}
        """
        # Initialize 3D list with zeros
        c = [[[0 for j in range(dim)] for i in range(dim)] for k in range(dim)]

        for i in range(dim):
            for j in range(i+1, dim):  # i < j only, then fill antisymmetric part
                coeffs = self.structure_tetrad(i, j)
                for k in range(dim):
                    c[k][i][j] = coeffs[k]
                    c[k][j][i] = -coeffs[k]  # antisymmetry enforced here
        return c

load("~/krphysics/sagemath/geometry.sage")
import time
t0 = time.time()
spacetime = Schwarzschild()
# Display the metric tensor
metric = g_S = spacetime.metric_tensor()
basis = es = spacetime.ortho_normal_frame()
forms = fs = spacetime.ortho_normal_tetrads()

eta = metric_minkowski = diagonal_matrix([1, 1, 1, -1])
c = spacetime.compute_structure_coefficients()
gamma = all_upper_coefficients(eta, c, g_inv = "eta", manifold = spacetime.neighborhood)
Riem = Riemannian_curvature(basis, gamma, c, dim=dim)
Ric = Ricci_tensor(Riem)
R = scalar_curvature(Ric, eta)
show(R.display())
t1= time.time()
print(t1 - t0)
# Riemann_curvature = g_S.riemann()
# show_commutators(es)
# show_metric_values(g_S, es)
# pairing_forms_vectors(fs, es)
# show_first_christoffel_symbols(metric, basis)
# show_structure_coefficients(basis, forms)
def show_lower_christoffel():
    for i in range(4):  # Loop over all i, j, mu
        for L in range(4):
            for j in range(4):
                print(lower_connection_term(i, L, j, metric, basis, forms).display())