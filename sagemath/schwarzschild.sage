class Schwarzschild:
    # Constants
    G = var('G')    # Gravitational constant
    M = var('M')    # Mass of the central object
    c = var('c')    # Speed of light
    def __init__(self,
        args_name = 't r theta phi',
        manifold = Manifold(4, 'Schwarzschild', structure='Lorentzian'),
        has_radius = True,
        neighborhood_name = 'U',
        ):
        if has_radius:
            self.radius = var('R_S')  # Schwarzschild radius
        else:
            self.radius = 2 * G * M / c**2  # Schwarzschild radius
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
        U = spacetime.spherical_chart().domain()# Open subset U of the 4-dimensional Lorentzian manifold Schwarzschild
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

spacetime = Schwarzschild()
# Display the metric tensor
g_S = spacetime.metric_tensor()
es = spacetime.ortho_normal_frame()
fs = spacetime.ortho_normal_tetrads()
dim = 4
def show_vectors(fs):
    print("Tetrads")
    for i in range(dim):
        print(fs[i].display())
    print("Exterior derivatives of tetrads")
    for i in range(dim):
        print(diff(fs[i]).display())
def show_commutators(es):
    print("Commutators of orthonormal basis")
    for i in range(dim):
        for j in range(dim):
            print((i, j), es[i].bracket(es[j]).display())

show_commutators(es)
# Riemann_curvature = g_S.riemann()
