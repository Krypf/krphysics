# https://doc.sagemath.org/html/en/reference/manifolds/sage/manifolds/differentiable/vectorfield.html
class Sphere:
    dim = 2
    def __init__(self,
        args_name = 'theta phi',
        manifold = Manifold(dim, 'Sphere', structure='Riemannian'),
        radius_name = 'a',
        neighborhood_name = 'U',
        ):
        self.radius = var('a')
        # Define the manifold
        self.manifold = manifold
        self.args_name = args_name
        self.args = var(args_name)# = theta, phi 
        # the open set covered by spherical coord.
        self.neighborhood = manifold.open_subset(neighborhood_name) # = U

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
        # print(g_sphere.display())
        return g_sphere

    def ortho_normal_frame(self):
        eth, eph = S2.coordinate_frame()
        e1 = eth
        e2 = eph / sin(theta)
        return [e1, e2]
    def ortho_normal_tetrads(self):
        fth = S2.one_form(1, 0, name='wth')
        fph = S2.one_form(0, sin(theta), name='wth')
        return [fth, fph]

S2 = Sphere()
g = S2.metric_tensor(); print(g_sphere.display());
eth, eph = S2.coordinate_frame()
e1, e2 = S2.ortho_normal_frame()

print(g.connection().display())

commutator_e1_e2 = e1.bracket(e2)
print(commutator_e1_e2.display())
