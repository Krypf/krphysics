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

    def __init__(self, manifold,
                 args_name='t r theta phi',
                 neighborhood_name='U'):
        # Define the radial coordinate
        t, r, th, ph = var(args_name)

        # Define f(r) and g(r)
        f_r = 1 - (2 * self.G * self.M) / (self.c^2 * r) + (self.G * self.Q^2) / (self.c^4 * r^2)
        g_r = 1 / f_r

        # Initialize parent (SphericalMetric)
        super().__init__(manifold, f_r=f_r, g_r=g_r,
                         args_name=args_name,
                         neighborhood_name=neighborhood_name)

    # Optionally override for clarity
    def potential(self, number_radius=1):
        r = self.args[number_radius]
        F_r = 1 - (2 * self.G * self.M) / (self.c^2 * r) + (self.G * self.Q^2) / (self.c^4 * r^2)
        return F_r

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
    print(g_RN.display())
    partials = spacetime.coordinate_frame()
    show(partials)
    G_type11 = g_RN.ricci().up(g_RN, 0)
    print(G_type11.display())
    # frame = spacetime.compute_orthonormal_frame(partials=partials)
    # eta = spacetime.sign
    # C = spacetime.compute_structure_coefficients(frame) # capital letter
    # gamma = all_upper_coefficients(eta, c, g_inv = "eta", domain = spacetime.neighborhood)
    # Riem, Ric, R = spacetime.compute_curvatures(frame, eta, C, g_inv = eta, domain=U)

main()