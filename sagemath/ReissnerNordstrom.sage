###########################################################
#  ReissnerNordstrom.sage
#  Reissner–Nordström spacetime metric class
#  (inherits from SphericalMetric)
###########################################################

load("~/krphysics/sagemath/spacetime.sage")
load("~/krphysics/sagemath/spherical_metric.sage")

class ReissnerNordstrom(SphericalMetric):
    # Physical constants
    G = var('G')      # Gravitational constant
    M = var('M')      # Mass of the central body
    Q = var('Q')      # Electric charge
    c = var('c')      # Speed of light
    k_e = var('k_e')

    def __init__(self, manifold,
                 args_name='t r theta phi',
                 neighborhood_name='U'):
        # Define the radial coordinate
        t, r, th, ph = var(args_name)

        # Define f(r) and g(r)
        f_r = 1 - (2 * self.G * self.M) / (self.c^2 * r) + (self.k_e * self.G * self.Q^2) / (self.c^4 * r^2)
        g_r = 1 / f_r

        # Initialize parent (SphericalMetric)
        super().__init__(manifold, f_r=f_r, g_r=g_r,
                         args_name=args_name,
                         neighborhood_name=neighborhood_name)
    def global_domain(self):
        U = self.spherical_chart().domain()# Open subset U of the 4-dimensional Lorentzian manifold 
        return U

    # Optionally override for clarity
    def potential(self, number_radius=1):
        r = self.args[number_radius]
        F_r = 1 - (2 * self.G * self.M) / (self.c^2 * r) + (self.G * self.Q^2) / (self.c^4 * r^2)
        return F_r
    def elemag_potential(self):
        # Define Coulomb's constant as a new variable
        U = self.spherical_chart().domain()# not just self.neighborhood
        # 1. Define the 4-potential A (a 1-form)
        A = U.tensor_field(0, 1, name='potential')
        # A = A.set_comp(self.chart.frame())
        A[0] = -k_e * Q / r  # Set the A_t component
        return A
    def field_strength_coord(self):
        A = self.elemag_potential()
        F = A.exterior_derivative() #.components()
        return F
    def energy_momentum_tensor(self):
        print("Returns the energy momentum tensor...")
        F = matrix(self.field_strength_coord().components())
        g = matrix(self.metric_tensor().components())
        g_inv = matrix(self.metric_tensor().inverse().components())
        trace_term = 1 / 4 * (F.transpose() * g_inv * F * g_inv).trace() * g # not / 4
        four_pi_T = F * g_inv * F.transpose() - trace_term
        # multiplied by 4 \pi
        return four_pi_T
    def energy_momentum_tensor_11(self):
        print("Returns the energy momentum tensor...")
        F = matrix(self.field_strength_coord().components())
        g = matrix(self.metric_tensor().components())
        g_inv = matrix(self.metric_tensor().inverse().components())
        ID = identity_matrix(4)
        trace_term = 1 / 4 * (F.transpose() * g_inv * F * g_inv).trace() * ID # not / 4
        four_pi_T = g_inv * F * g_inv * F.transpose() - trace_term
        # multiplied by 4 \pi
        return four_pi_T

###########################################################
# Example usage:
def set_const():
    _dim = 4
    M = Manifold(_dim, 'ReissnerNordstrom', structure='Lorentzian')
    RN = ReissnerNordstrom(M)
    
    return RN
###########################################################
def main():
    t0 = time.time()
    spacetime = set_const()
    # U = spacetime.global_domain()
    metric = g_RN = spacetime.metric_tensor()
    # Display the metric tensor
    print("metric tensor:")
    print(g_RN.display())
    # g_mat = matrix(RN.metric_tensor().components())
    partials = spacetime.coordinate_frame()
    # show(partials)d
    G_type11 = g_RN.ricci().up(g_RN, 0)
    print("Ricci ^i_j tensor")
    print(G_type11.display())

    four_pi_T = spacetime.energy_momentum_tensor_11()
    for i in range(4): 
        for j in range(4): 
            print(four_pi_T[i, j].display())
    t1 = time.time()
    
    # frame = spacetime.compute_orthonormal_frame(partials=partials)
    # eta = spacetime.sign
    # C = spacetime.compute_structure_coefficients(frame) # capital letter
    # gamma = all_upper_coefficients(eta, c, g_inv = "eta", domain = spacetime.neighborhood)
    # Riem, Ric, R = spacetime.compute_curvatures(frame, eta, C, g_inv = eta, domain=U)


RN = spacetime = set_const()
# mode = "orthonormal"
mode = "coordinate"
K = compute_second_fundamental(spacetime, mode) # r = 3 R_S / 2
g = spacetime.metric_tensor()
_ans = print_second_fundamental_form(spacetime, mode, K)
# main()
