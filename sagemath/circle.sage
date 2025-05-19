class CircleCoordinates:
    # Constants
    dim = 2
    def __init__(self,
        args_name = 'r theta',
        manifold = Manifold(dim, 'Euclid', structure='Riemannian'),
        neighborhood_name = 'U',
        ):
        # Define the manifold
        self.manifold = manifold
        self.args_name = args_name
        self.args = var(args_name)
        self.neighborhood = manifold.open_subset(neighborhood_name) # = U

    def spherical_chart(self):
        # Restrict the domain for the variables in the chart
        chart_expression = f"{self.args_name.split()[0]}:(0,oo):r {self.args_name.split()[1]}:(0,2*pi):theta"
        C_spherical = self.neighborhood.chart(chart_expression)
        return C_spherical

    def coordinate_frame(self):
        return self.spherical_chart().frame()

    def dual_frame(self):
        return self.spherical_chart().coframe()
    def global_domain(self):
        U = self.spherical_chart().domain()
        return U
    def metric_tensor(self):
        U = self.global_domain()
        g_E = U.metric('g_E')
        # Define the metric tensor
        g_E[0, 0] = 1
        g_E[1, 1] = r^2
        print(g_E.display())
        return g_E


E = CircleCoordinates()
g_Euclid_spherical = g_E = E.metric_tensor()
