# Parent class Spacetime defined
# Subclasses are Sphere, Schwarzschild
load("~/krphysics/sagemath/geometry.sage")

class Spacetime():
    def __init__(self, manifold):
        # Define the self
        self.manifold = manifold
        self.dimension = manifold.dimension()

    def metric_tensor(self):
        raise NotImplementedError("Spacetime must implement metric")
    
    def E_hat(self):
        raise NotImplementedError("E_hat must implement metric")
    
    def E_inv(self):
        raise NotImplementedError("E_inv must implement metric")
    
    def compute_orthonormal_frame(self, chart = None, partials = None, _show = False):
        """
        Computes and displays the transformed vector fields forming an orthonormal frame.
        
        Parameters:
            self.chart: The chart or manifold coordinate frame (with .frame() method)
            self = solution: An object with an .E_inv() method 
        Returns:
            List of vector fields forming an orthonormal (moving) frame.
        """
        orthonormal_frame = []
        if chart is None and partials is None:
            raise NotImplementedError("Select some chart or moving frame")
        if partials is None:
            partials = chart.frame()

        yo = self.E_inv() # katakana letter ﾖ
        n = self.dimension
        for j in range(n):
            vec = sum(partials[i] * yo[i, j] for i in range(n))
            orthonormal_frame.append(vec)
            
            if _show:
                show(vec.display())
        return orthonormal_frame

    def compute_orthonormal_tetrads(self, chart = None, domain = None, _show = False, name=None):
    # Get tetrads as 1-forms: f^a = E^a_μ dx^μ
        if chart is None and domain is None:
            raise NotImplementedError("Select some chart or domain")
        if name is None:
            name = "Tetrad"
        if domain is None:
            domain = chart.domain()
        
        tetrads = []
        E = self.E_hat()
        n = self.dimension
        for a in range(n):
            dx = domain.one_form(E.row(a), name=f"{name}^{a}")
            if _show:
                show(dx.display())
            tetrads.append(dx)
        return tetrads

    def connection_one_forms(self, gamma, **kwargs):
        N = self.dimension
        omega = [[0 for j in range(N)] for i in range(N)]
        for i in range(N):
            for j in range(N):
                omega[i][j] = sum(gamma[i][l][j] * self.compute_orthonormal_tetrads(**kwargs)[l] for l in range(N)) # change omega_forms
                show(omega[i][j].display())
        return omega

    def show_derivative_tetrads(self, **kwargs):
        # Display the tetrads and their exterior derivatives
        tetrads = self.compute_orthonormal_tetrads(**kwargs)
        for x in tetrads:
            print(x.display())
            show(x.exterior_derivative().display())

    def structure_coefficient_tetrad(self, frame, i, j):
        com = commutator_field(frame, i, j)
        E = self.E_hat()
        n = self.dimension
        cs = [0 for _ in range(n)]
        for k in range(n):
            cs[k] = sum(com[a] * E[k][a] for a in range(n))
        return cs

    def compute_structure_coefficients(self, frame):
        """
        Computes the full structure constants tensor c[k][i][j],
        where c^k_{ij} satisfies c^k_{ji} = -c^k_{ij} (antisymmetry).

        Returns:
        - c: a 3D list [k][i][j] representing c^k_{ij}
        """
        # Initialize 3D list with zeros
        n = self.dimension
        c = [[[0 for j in range(n)] for i in range(n)] for k in range(n)]

        for i in range(n):
            for j in range(i+1, n):  # i < j only, then fill antisymmetric part
                coeffs = self.structure_coefficient_tetrad(frame, i, j)
                for k in range(n):
                    c[k][i][j] = coeffs[k]
                    c[k][j][i] = -coeffs[k]  # antisymmetry enforced here
        return c

    def compute_curvatures(self, frame, metric, structure, g_inv = "eta", domain = None):
        t0 = time.time()
        _dim = self.dimension
        g = metric
        c = structure
        gamma = all_upper_coefficients(g, c, g_inv = g_inv, domain = domain)

        Riem = Riemannian_curvature(frame, gamma, c)
        show_four_tensor(Riem, _dim)

        Ric = Ricci_tensor(Riem)
        show(show_two_tensor(Ric, _dim))
        R = scalar_curvature(Ric, metric)
        show(R.display())
        t1 = time.time(); print(f"compute_curvatures in {t1 - t0} seconds");

        return Riem, Ric, R

class Minkowski(Spacetime):
    var('c')
    def __init__(self, manifold, args_name = 'x', neighborhood_name = 'U', eta00 = -1):
        # Define the manifold
        super().__init__(manifold)
        self.args_name = ' '.join([f'{args_name}{i}' for i in range(self.dimension)])
        self.args = var(self.args_name)
        self.chart = self.manifold.chart(self.args_name)

        self.neighborhood = manifold.open_subset(neighborhood_name) # = U
        self.sign_matrix = (- eta00) * diagonal_matrix([- 1, 1, 1, 1]) # Minkowski metric # 0 component is negative.

    def metric_tensor(self):
        _frame = self.chart.frame()
        n = self.dimension
        eta = M.metric('eta', signature=(n-1, 1, 0))
        for i in range(n):
            eta[_frame, i, i] = self.sign_matrix[i, i]
        return eta

    def electromagnetic(self):
        dom = self.chart.domain()
        F = dom.diff_form(2)
        x = self.chart
        E = [function('E1')(*x), function('E2')(*x), function('E3')(*x)]
        for i in range(1, 3+1):
            F[i,0] = E[i-1] / c
        # B = [function('B1')(*x), function('B2')(*x), function('B3')(*x)]
        F[1,2] = function('B12')(*x) # B[3-1]  # dx∧dy
        F[2,3] = function('B23')(*x) # B[1-1]  # dy∧dz
        F[3,1] = function('B31')(*x) # B[2-1]  # dz∧dx
        return F     
    def excitation(self):
        dom = self.chart.domain()
        exc = dom.diff_form(2)
        x = self.chart
        H = [function('H1')(*x), function('H2')(*x), function('H3')(*x)]
        for i in range(1, 3+1):
            exc[i,0] = - H[i-1] / c # magnetic field
            # with the minus sign
        # D = [function('E1')(*x), function('E2')(*x), function('E3')(*x)]
        exc[1,2] = function('D12')(*x) # D[3-1] # dx∧dy
        exc[2,3] = function('D23')(*x) # D[1-1] # dy∧dz
        exc[3,1] = function('D31')(*x) # D[2-1] # dz∧dx
        return exc     

n = 4
M = Manifold(n, 'Minkowski', structure='Lorentzian')
Min = Minkowski(M)
eta = Min.metric_tensor()
N = Manifold(3, 'de Sitter', ambient=M)
# class EuclidSpace(Spacetime):