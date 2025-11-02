load("~/krphysics/sagemath/spacetime.sage")
class Schwarzschild(Spacetime):
    # Constants
    G = var('G')    # Gravitational constant
    Mass = var('M_S')    # Mass of the central object
    c = var('c')    # Speed of light

    def __init__(self, manifold,
        args_name = 't r theta phi',
        has_radius = True,
        neighborhood_name = 'U',
        ):
        if has_radius:
            self.radius = var('R_S')  # Schwarzschild radius
        else:
            self.radius = 2 * G * Mass / c**2  # Schwarzschild radius
        # Define the manifold
        super().__init__(manifold)
        self.args_name = args_name
        self.args = var(args_name)# = t, r, theta, phi 
        self.neighborhood = manifold.open_subset(neighborhood_name) # = U
        self.sign_matrix = (+1) * diagonal_matrix([- 1, 1, 1, 1]) # Minkowski metric # 0 component is negative.
    
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
        F_r = 1 - self.radius / r
        return F_r

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

    def E_inv(self):
        r, theta, F = self.components()
        # Construct the matrix
        E_inverse = matrix([
            [1 / sqrt(F), 0, 0, 0],
            [0, sqrt(F), 0, 0],
            [0, 0, 1 / r, 0],
            [0, 0, 0, 1 / (r * sin(theta))]
        ])
        return E_inverse

    def ortho_normal_tetrads(self):
        U = self.global_domain()
        omega = self.compute_orthonormal_tetrads(domain=U, _show = True)
        return omega

def set_const():
    _dim = 4
    M = Manifold(_dim, 'Schwarzschild', structure='Lorentzian')
    spacetime = Schwarzschild(M)
    return spacetime

def compute_curvatures():
    t0 = time.time()
    spacetime = set_const()
    U = spacetime.global_domain()
    # Display the metric tensor
    # metric = g_S = spacetime.metric_tensor()
    partials = spacetime.coordinate_frame()
    frame = spacetime.compute_orthonormal_frame(partials=partials)
    eta = spacetime.sign
    C = spacetime.compute_structure_coefficients(frame) # capital letter
    # gamma = all_upper_coefficients(eta, c, g_inv = "eta", domain = spacetime.neighborhood)

    Riem, Ric, R = spacetime.compute_curvatures(frame, eta, C, g_inv = eta, domain=U)
    t1 = time.time(); print(f"compute_curvatures in {t1 - t0} seconds");
    return Riem, Ric, R



spacetime = set_const()
# mode = "orthonormal"
mode = "coordinate"
K = compute_second_fundamental(spacetime, mode) # r = 3 R_S / 2
g = spacetime.metric_tensor()
_ans = print_second_fundamental_form(spacetime, mode, K)

# x = spacetime.spherical_chart()
# U = spacetime.global_domain()
# surface = U.scalar_field(x[1]) # r = const.

# Riemann_curvature = g_S.riemann()
# show_commutators(frame)
# forms = spacetime.ortho_normal_tetrads()
# show_metric_values(g_S, frame)
# pairing_forms_vectors(forms, frame)
# show_first_christoffel_symbols(metric, frame)
# show_structure_coefficients(frame, forms)
# main()
