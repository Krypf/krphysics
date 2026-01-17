###########################################################
#  VaidyaMetric.sage
#  Class for the Vaidya spacetime (outgoing radiation)
#  Written in a modular style similar to SphericalMetric
###########################################################

load("~/krphysics/sagemath/spacetime.sage")

class VaidyaMetric(Spacetime):
    def __init__(self, manifold, args_name='u r th ph', neighborhood_name='U'):
        # Step 1: Initialize the spacetime
        super().__init__(manifold)
        self.args_name = args_name
        self.args = var(args_name)  # (u, r, th, ph)
        self.neighborhood = manifold.open_subset(neighborhood_name)
        self.chart = self.neighborhood.chart(r'u r:(0,oo) th:(0,pi):\theta ph:(0,2*pi):\phi')

        # Step 2: Define symbolic mass function m(u)
        self.u = self.args[0]
        self.r = self.args[1]
        self.th = self.args[2]
        self.ph = self.args[3]
        self.m = function('m')(self.u)

        # Step 3: Define the Lorentzian metric
        self.g = self.neighborhood.metric('g')
        self._set_metric_components()

    def _set_metric_components(self, null = +1):
        u, r, th = self.u, self.r, self.th

        g = self.g
        m = self.m

        # Metric components for ingoing Vaidya metric
        g[0,0] = -(1 - 2*m/r)    # g_uu
        g[0,1] = null            # g_ur = null
        g[2,2] = r^2             # g_thth
        g[3,3] = r^2 * sin(th)^2 # g_phph

    def metric_tensor(self):
        return self.g

    def mass_function(self):
        return self.m

    def compute_curvature(self):
        Ric = self.g.ricci()
        R_scalar = self.g.ricci_scalar()
        return Ric, R_scalar

    def display_metric(self):
        print("\n--- Metric g ---")
        print(self.g.display())

    def delta(self):
        g = self.g
        Delta = g[0,0] * g[1,1] - g[0,1]^2
        return Delta

def compute_second_fundamental(spacetime, mode):
    x = spacetime.chart
    U = spacetime.neighborhood
    surface = U.scalar_field(x[1] - 3 * spacetime.m) # U.scalar_field(x[1]) # r = const.
    partials = spacetime.chart.frame()
    _dim = spacetime.dimension
    if mode == "orthonormal":
        frame = spacetime.compute_orthonormal_frame(partials=partials)
        normal_form = [frame[mu](surface) for mu in range(_dim)]
        eta = spacetime.sign_matrix
        C = spacetime.compute_structure_coefficients(frame)
        gamma = all_upper_coefficients(eta, C, g_inv = eta, domain = U)
    if mode == "coordinate": 
        g = spacetime.metric_tensor()
        gamma = g.christoffel_symbols()
        gamma = gamma[:]
        frame = partials
        normal_form = [frame[mu](surface) for mu in range(_dim)]
    K = second_fundamental_form(gamma, normal_form)
    return K

def print_second_fundamental_form(spacetime, mode, K):
    g = spacetime.metric_tensor()
    print(g.display())
    _dim = spacetime.dimension
    _ans = [[0 for i in range(_dim)] for j in range(_dim)]
    for i in range(_dim): 
        for j in range(_dim): 
            if K[i][j] != 0:
                print(i, j)
                print("second_fundamental_form K:")
                print(K[i][j].display())
                # if mode == "orthonormal":
                #     _ans[i][j] = K[i][j] / K[3][3]
                #     print("K_[3][3] is set to be one", _ans[i][j].display())
                    
                if mode == "coordinate": 
                    _ans[i][j] = K[i][j] / g[i, j]
                    print("K[i][j] / g[i, j] = ", _ans[i][j].display())
                    print("metric g:")
                    print(g[i, j].display())
                    print(latex(_ans[i][j].display()))
    return _ans

_dim = 4
M = Manifold(_dim, 'M', structure='Lorentzian')
v = spacetime = VaidyaMetric(M)

g = v.metric_tensor()
Ric, R = v.compute_curvature()
G = v.einstein11type()
mode = "coordinate"

# K = compute_second_fundamental(spacetime, mode)
# _ans = print_second_fundamental_form(spacetime, mode, K)
