load("~/krphysics/sagemath/spacetime.sage")

class EuclidSpace(Spacetime):
    def __init__(self, manifold, args_name='x', neighborhood_name='U'):
        # Initialize the parent Spacetime class
        super().__init__(manifold)
        
        # Define the coordinate arguments based on dimension
        self.args_name = ' '.join([f'{args_name}{i}' for i in range(self.dimension)])
        self.args = var(self.args_name)
        self.chart = self.manifold.chart(self.args_name)
        
        # Define the subset domain
        self.neighborhood = manifold.open_subset(neighborhood_name)
        
        # Isolated setup for charts and transitions
        self.cartesian_chart()
        self.spherical_chart()
        self.set_transition_map()
        
        # Define the default sign matrix for Euclidean metric
        self.sign_matrix = diagonal_matrix([1] * self.dimension)

    def cartesian_chart(self):
        """Isolated definition of the Cartesian chart."""
        # Define the coordinate chart
        C_cartesian = self.neighborhood.chart(self.args_name)
        return C_cartesian

    def spherical_chart(self):
        """Isolated definition of the Spherical chart."""
        C_spherical = self.neighborhood.chart(r'r:(0,+oo) th:(0,pi):periodic ph:(0,2*pi):periodic')
        self.r, self.th, self.ph = C_spherical[:]
        return (C_spherical)

    def set_transition_map(self):
        """Isolated definition of the Spherical -> Cartesian transition."""
        S = self.spherical_chart()
        C = self.cartesian_chart()
        return S.transition_map(C, [
            self.r * sin(self.th) * cos(self.ph), 
            self.r * sin(self.th) * sin(self.ph), 
            self.r * cos(self.th)
        ])

    def metric_tensor(self):
        """
        Defines the Euclidean metric delta_ij.
        """
        _frame = self.chart.frame()
        n = self.dimension
        # Use signature (n, 0, 0) for a positive-definite Euclidean metric
        delta = self.manifold.metric('delta', signature=(n, 0, 0))
        for i in range(n):
            delta[_frame, i, i] = self.sign_matrix[i, i]
        return delta

    def E_hat(self):
    # Construct the matrix
        E = matrix([
            [1, 0,                      0],
            [0,      self.r,  0],
            [0,      0,                      self.r * sin(self.th)],
        ])
        return E

    def E_inv(self):
        # Construct the matrix
        E = matrix([
            [1, 0,                            0],
            [0,          1 / (self.r),  0],
            [0,          0,                            1 / (self.r * sin(self.th))],
        ])
        return E
        
    def diffeomorphisms(self):
        # Define the diffeomorphisms for Cartesian and Spherical charts
        C_cartesian = self.cartesian_chart()
        C_spherical = self.spherical_chart()
        
        # Transition map from Spherical to Cartesian
        transition_map = self.set_transition_map()
        
        return C_cartesian, C_spherical, transition_map

    def get_killing_vectors(self):
        X, S, phi = self.diffeomorphisms()
        x, y, z = X
        e_x, e_y, e_z = X.frame()
        # better written explicitly:
        Lz = x*e_y - y*e_x
        Lx = y*e_z - z*e_y
        Ly = z*e_x - x*e_z
        show(Lx.display(S))
        show(Ly.display(S))
        show(Lz.display(S))
        return Lx, Ly, Lz

def initialization(n = 3):
    # Define a 3D manifold for Euclidean space
    Eucdim = Manifold(n, 'Euclidean3', structure='Riemannian')
    Euc = EuclidSpace(Eucdim)
    
    return Euc

Euc = initialization()
Lx, Ly, Lz = Euc.get_killing_vectors()
# Now you can use existing methods like:
S = Euc.spherical_chart()
orthonormal_frame = Euc.compute_orthonormal_frame(chart=S, _show=False)

# delta = Euc.metric_tensor()
# main()