# Parent class Spacetime defined
# Subclasses are Sphere, Schwarzschild

class Spacetime():
    def __init__(self, manifold):
        # Define the self
        self.manifold = manifold
        self.dimension = manifold.dimension()

    def metric_tensor(self):
        raise NotImplementedError("Spacetime must implement metric")
    
    def compute_orthonormal_frame(self, partials = None, _show = False):
        """
        Computes and displays the transformed vector fields forming an orthonormal frame.
        
        Parameters:
            self.chart: The chart or manifold coordinate frame (with .frame() method)
            self = solution: An object with an .E_inv() method 
        Returns:
            List of vector fields forming an orthonormal (moving) frame.
        """
        orthonormal_frame = []
        if partials is None:
            partials = self.chart.frame()

        yo = self.E_inv() # katakana letter ﾖ
        n = self.dimension
        for j in range(n):
            vec = sum(partials[i] * yo[i, j] for i in range(n))
            orthonormal_frame.append(vec)
            
            if _show:
                show(vec.display())
        return orthonormal_frame

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

        print("Calculating the Riemannian_curvature... ")
        Riem = Riemannian_curvature(frame, gamma, c)
        t1 = time.time(); print(f"in {t1 - t0} seconds");
        show_four_tensor(Riem, _dim)

        print("Calculating the Ricci_tensor and the scalar_curvature... ")
        Ric = Ricci_tensor(Riem)
        
        R = scalar_curvature(Ric, metric)
        t2 = time.time(); print(f"in {t2 - t1} seconds");
        show(show_two_tensor(Ric, _dim))
        show(R.display())

        return Riem, Ric, R
