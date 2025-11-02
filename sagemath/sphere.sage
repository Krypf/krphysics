# https://doc.sagemath.org/html/en/reference/manifolds/sage/manifolds/differentiable/vectorfield.html
load("spacetime.sage")

class Sphere(Spacetime):
    def __init__(self, manifold,
        radius_name = 'a',
        args_name = 'theta phi',
        neighborhood_name = 'U',
        ):
        # Define the manifold
        super().__init__(manifold)
        self.radius = var('a')
        self.args_name = args_name
        self.args = var(args_name)# = theta, phi 
        # the open set covered by spherical coord.
        self.neighborhood = manifold.open_subset(neighborhood_name) # = U
        self.sign_matrix = diagonal_matrix([1, 1])

    def spherical_chart(self):
        # Restrict the domain for the variables in the chart
        chart_expression = f"{self.args_name.split()[0]}:(0,pi):theta {self.args_name.split()[1]}:(-pi,pi):phi"
        C_spherical = self.neighborhood.chart(chart_expression)
        return C_spherical
        #XS.<theta, phi> = U.chart(r'θ:(0,pi):\theta ϕ:(0,2*pi):\phi')
    
    def coordinate_frame(self):
        return self.spherical_chart().frame()
    def dual_frame(self):
        return self.spherical_chart().coframe()
    def global_domain(self):
        U = self.spherical_chart().domain()
        return U
    
    def metric_tensor(self):
        U = self.global_domain()
        g_sphere = U.metric('g_sphere')
        # Define the metric tensor
        g_sphere[0, 0] = self.radius^2
        g_sphere[1, 1] = (self.radius * sin(theta))^2
        return g_sphere

    def E_hat(self):
        theta = self.spherical_chart()[0]
        # Construct the matrix
        E = matrix([
            [self.radius, 0],
            [0, self.radius * sin(theta)],
        ])
        return E

    def E_inv(self):
        theta = self.spherical_chart()[0]
        # Construct the matrix
        E = matrix([
            [1 / self.radius, 0],
            [0, 1 / (self.radius * sin(theta))],
        ])
        return E

    def ortho_normal_tetrads(self):
        U = self.global_domain()
        omega = self.compute_orthonormal_tetrads(domain=U, _show = True)
        return omega
    
    def show_fields(self):
        g_sphere = self.metric_tensor();
        print(g_sphere.display())
        print(g_sphere.connection().display())

def set_const():
    _dim = 2
    M = Manifold(_dim, 'Sphere', structure='Riemannian')
    S2 = Sphere(M)
    return S2

def main():
    S2 = set_const()
    g_sphere = S2.metric_tensor()
    partials = S2.coordinate_frame()
    frame = S2.compute_orthonormal_frame(partials = partials)
    """ frame
    eth, eph = S2.coordinate_frame()
    e1 = eth / self.radius
    e2 = eph / sin(theta) / self.radius # needs a 
    """

    c = S2.compute_structure_coefficients(frame)
    delta = S2.sign_matrix
    Riem, Ric, R = S2.compute_curvatures(frame, delta, c, g_inv = delta, domain=S2.global_domain())
    return 0

r = var('r')
S2_r = manifolds.Sphere(2, radius=r); S2_r
K = S2_r.second_fundamental_form()
K.display()
# main()