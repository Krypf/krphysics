# https://gemini.google.com/share/987375768a6b

load("~/krphysics/sagemath/spacetime.sage")

class Kerr(Spacetime):
    # Constants
    G = var('G')       # Gravitational constant
    Mass = var('M')    # Mass of the black hole
    a = var('a')       # Spin parameter (angular momentum per unit mass, J/M)
    c = var('c')       # Speed of light

    def __init__(self, manifold, 
                 args_name='t r theta phi', 
                 neighborhood_name='U'):
        
        # Define the manifold via the parent Spacetime class
        super().__init__(manifold)
        
        # Set up coordinate arguments
        self.args_name = args_name
        self.args = var(args_name)  # t, r, theta, phi 
        
        # Define the open subset (global domain for Boyer-Lindquist coordinates)
        self.neighborhood = manifold.open_subset(neighborhood_name) 
        
        # Minkowski sign convention (- + + +)
        self.sign_matrix = (+1) * diagonal_matrix([-1, 1, 1, 1]) 

    def spherical_chart(self):
        # Define the coordinate chart for Boyer-Lindquist coordinates
        C_spherical = self.neighborhood.chart(self.args_name)
        return (C_spherical)

    def coordinate_frame(self):
        return self.spherical_chart().frame()

    def dual_frame(self):
        return self.spherical_chart().coframe()

    def global_domain(self):
        U = self.spherical_chart().domain() # Open subset U of the 4-dimensional Lorentzian manifold
        return U
    def components(self):
        # Extract coordinates
        r = self.args[1]
        theta = self.args[2]
        
        # Extract parameters (assuming self.Mass and self.a are defined in __init__)
        M = self.Mass
        a = self.a
        return r, theta, M, a

    # Metric components
    def metric_tensor(self):
        r, theta, M, a = self.components()
        
        # Define Kerr metric functions
        Sigma = r^2 + a^2 * cos(theta)^2
        Delta = r^2 - 2 * M * r + a^2
        
        U = self.global_domain()
        g_K = U.metric('g_K')
        
        # Define the metric tensor components
        g_K[0, 0] = -(1 - 2 * M * r / Sigma)
        g_K[0, 3] = -2 * M * a * r * sin(theta)^2 / Sigma
        g_K[3, 0] = g_K[0, 3]
        g_K[1, 1] = Sigma / Delta
        g_K[2, 2] = Sigma
        g_K[3, 3] = (r^2 + a^2 + 2 * M * a^2 * r * sin(theta)^2 / Sigma) * sin(theta)^2
        
        return g_K

def set_const():
    """
    Initializes the 4-dimensional Lorentzian manifold and returns 
    the instantiated Kerr spacetime object.
    """
    _dim = 4
    # Create the Lorentzian manifold for Kerr spacetime
    M = Manifold(_dim, 'Kerr_Manifold', structure='Lorentzian')
    
    # Instantiate the Kerr class
    spacetime = Kerr(M)
    
    return spacetime


def compute_curvatures():
    t0 = time.time()
    spacetime = set_const()
    U = spacetime.global_domain()
    # Display the metric tensor
    metric = g_K = spacetime.metric_tensor()
    partials = spacetime.coordinate_frame()
    frame = partials
    # frame = spacetime.compute_orthonormal_frame(partials=partials)
    eta = spacetime.sign_matrix
    C = spacetime.compute_structure_coefficients(frame) # capital letter
    # gamma = all_upper_coefficients(eta, c, g_inv = "eta", domain = spacetime.neighborhood)

    Riem, Ric, R = spacetime.compute_curvatures(frame, eta, C, g_inv = eta, domain=U)
    t1 = time.time(); print(f"compute_curvatures in {t1 - t0} seconds");
    return Riem, Ric, R

# main()
K = set_const()
g = K.metric_tensor()
g_inv = g.inverse()
print(g.display())
partials = K.coordinate_frame()
# g.ricci().display() # zero
# g.riemann().display() # too heavy to show in the terminal.

# Assuming 'U' is the global domain and 'r', 'a', 'Delta' are defined
U = K.global_domain()
r, theta, M, a = K.components()
Sigma = r^2 + a^2 * cos(theta)^2
Delta = r^2 - 2 * M * r + a^2

# Define the null vectors l and n in the coordinate basis
l = U.vector_field([r^2 + a^2, -Delta, 0, a], frame=partials, name='l')
n = U.vector_field([r^2 + a^2, Delta, 0, a], frame=partials, name='n')
# g(l,l) == 0 # True
# g(n,n) == 0 # True
ln_sym = (l * n + n * l) / 2
K_tensor = (1 / (Delta)) * ln_sym # + r^2 * g_inv
p_t, p_r, p_th, p_ph = var('p_t p_r p_th p_ph')
p_form = U.diff_form(1, 'p'); p_form.set_comp(partials)[:] = [p_t, p_r, p_th, p_ph]
# carter_constant = K_tensor.contract(0, p_form, 0).contract(0, p_form, 0)
# print(simplify(carter_constant).display()) # too complicated to show in the terminal.