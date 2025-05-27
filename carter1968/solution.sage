import time
load("~/krphysics/sagemath/geometry.sage")

var('Lambda e h')
m1, m2 = var('m1 m2')
n1, n2 = var('n1 n2')
e1, e2 = var('e1 e2')
n = var('n')

class CarterSolution(Spacetime):
    def __init__(self, manifold):
        super().__init__(manifold)
        self.chart = self.manifold.chart(r'lambda_ mu psi chi')
    
    def Delta_lambda(self):
        raise NotImplementedError("Subclasses must implement Delta_lambda")
    
    def Delta_mu(self):
        raise NotImplementedError("Subclasses must implement Delta_mu")

    def P_lambda(self):
        raise NotImplementedError("Subclasses must implement P_lambda")

    def P_mu(self):
        raise NotImplementedError("Subclasses must implement P_mu")
    
    def Q_lambda(self):
        raise NotImplementedError("Subclasses must implement Q_lambda")
    
    def Q_mu(self):
        raise NotImplementedError("Subclasses must implement Q_mu")

    def Z(self):
        return self.P_lambda() * self.Q_mu() - self.P_mu() * self.Q_lambda()

    def get_parameters(self):
        """ 
        Retrieve all required parameters and return them as a tuple.
        This allows easy access to multiple values at once.
        """
        return (
            self.Z(),              # Get Z
            self.Delta_lambda(),   # Get Delta_lambda
            self.Delta_mu(),       # Get Delta_mu
            self.P_lambda(),       # Get P_lambda
            self.P_mu(),           # Get P_mu
            self.Q_lambda(),       # Get Q_lambda
            self.Q_mu()            # Get Q_mu
        )

    def omega_forms(self):
        Z, Delta_lambda, Delta_mu, P_lambda, P_mu, Q_lambda, Q_mu = self.get_parameters() 

        # Define the 1-forms
        omega_p1 = M.diff_form(1, name='omega^{+1}')
        omega_p1.set_comp(chart.frame())[0] = sqrt(Z / Delta_lambda)  # lambda_ component

        omega_m1 = M.diff_form(1, name='omega^{-1}')
        omega_m1.set_comp(chart.frame())[1] = sqrt(Z / Delta_mu)  # mu component

        omega_p2 = M.diff_form(1, name='omega^{+2}')
        omega_p2.set_comp(chart.frame())[2] = (sqrt(Z * Delta_mu) / Z) * P_lambda  # psi component
        omega_p2.set_comp(chart.frame())[3] = -(sqrt(Z * Delta_mu) / Z) * Q_lambda  # chi component

        omega_m2 = M.diff_form(1, name='omega^{-2}')
        omega_m2.set_comp(chart.frame())[2] = (sqrt(Z * Delta_lambda) / Z) * P_mu  # psi component
        omega_m2.set_comp(chart.frame())[3] = -(sqrt(Z * Delta_lambda) / Z) * Q_mu  # chi component

        return [omega_p1, omega_m1, omega_p2, omega_m2]
    
    def metric_tensor(self):
        """ Calculate the metric ds^2 """
        g = self.manifold.metric(name="g")  # metric tensor g_ij
        Z, Delta_lambda, Delta_mu, P_lambda, P_mu, Q_lambda, Q_mu = self.get_parameters() 

        g[0, 0] = Z / Delta_lambda  # coefficient of dλ² 
        g[1, 1] = Z / Delta_mu      # coefficient of dμ² 
        g[2, 2] = Delta_mu / Z * P_lambda^2  # dψ²
        g[3, 3] = Delta_mu / Z * Q_lambda^2  # dχ²
        g[2, 3] = -Delta_mu / Z * P_lambda * Q_lambda  # dψ dχ (cross-term)

        g[2, 2] -= Delta_lambda / Z * P_mu^2  # add to the dψ²-term
        g[3, 3] -= Delta_lambda / Z * Q_mu^2  # add to the dχ²-term
        g[2, 3] += Delta_lambda / Z * P_mu * Q_mu  # add to the dψ dχ-term

        # symmetric
        g[3, 2] = g[2, 3]

        return g
    def E_hat(self):
        Z, Delta_lambda, Delta_mu, P_lambda, P_mu, Q_lambda, Q_mu = self.get_parameters() 

        E = matrix([
        [sqrt(Z/Delta_lambda), 0, 0, 0],
        [0, sqrt(Z/Delta_mu), 0, 0],
        [0, 0, (sqrt(Z*Delta_mu)/Z) * P_lambda, -(sqrt(Z*Delta_mu)/Z) * Q_lambda],
        [0, 0, (sqrt(Z*Delta_lambda)/Z) * P_mu, -(sqrt(Z*Delta_lambda)/Z) * Q_mu]
        ])
        return E

    def E_inv(self):
        Z, Delta_lambda, Delta_mu, P_lambda, P_mu, Q_lambda, Q_mu = self.get_parameters() 
        
        return matrix(SR, [
            [1/sqrt(Z / Delta_lambda), 0, 0, 0],
            [0, 1/sqrt(Z / Delta_mu), 0, 0],
            [0, 0, Q_mu/sqrt(Z * Delta_mu), -Q_lambda/sqrt(Z * Delta_lambda)],
            [0, 0, P_mu/sqrt(Z * Delta_mu), -P_lambda/sqrt(Z * Delta_lambda)]
        ])

    def elemag_potential(self):
        Z, Delta_lambda, Delta_mu, P_lambda, P_mu, Q_lambda, Q_mu = self.get_parameters() 
        X_lambda = e1 * self.chart[0]
        X_mu = e2 * self.chart[1]

        A = self.manifold.diff_form(1, name='potential')
        A.set_comp(self.chart.frame())[2] = (P_lambda * X_mu + P_mu * X_lambda) / Z # d_psi
        A.set_comp(self.chart.frame())[3] = - (Q_lambda * X_mu + Q_mu * X_lambda) / Z # d_chi
        return A

    def field_strength_coord(self):
        A = self.elemag_potential()
        F = A.exterior_derivative().components()
        return F
    def field_strength_tetrad(self):
        Einv = self.E_inv() # coordinates to tetrads
        N = self.dimension
        F = self.field_strength_coord()
        field_strength = [[0 for j in range(N)] for i in range(N)]
        for i in range(N):
            for j in range(N):
                field_strength[i][j] = sum(F[a, b] * Einv[a][i] * Einv[b][j] for a in range(N) for b in range(N))
        return field_strength

    def energy_momentum_tensor(self, g, g_inv):
        field_strength = self.field_strength_tetrad()
        F = matrix(field_strength)
        trace_term = 1 / 4 * (F.transpose() * g_inv * F * g_inv).trace() * g # not / 4
        four_pi_T = F * g_inv * F.transpose() - trace_term
        # multiplied by 4 \pi
        return four_pi_T

    def show(self):
        # Print specific computed expressions
        print(f"Show the methods of {type(self)}")
        print(f"Delta_lambda = {self.Delta_lambda()}")
        print(f"Delta_mu = {self.Delta_mu()}")
        print(f"P_lambda = {self.P_lambda()}")
        print(f"P_mu = {self.P_mu()}")
        print(f"Q_lambda = {self.Q_lambda()}")
        print(f"Q_mu = {self.Q_mu()}")

    def show_latex(self):
        # Print computed expressions in LaTeX format
        print(f"Show the latex representations of {type(self)}")

        print(f"Delta_lambda = {latex(self.Delta_lambda())}")
        print(f"Delta_mu = {latex(self.Delta_mu())}")
        print(f"P_lambda = {latex(self.P_lambda())}")
        print(f"P_mu = {latex(self.P_mu())}")
        print(f"Q_lambda = {latex(self.Q_lambda())}")
        print(f"Q_mu = {latex(self.Q_mu())}")

    # Display results
    def print_solution(self):
        self.show()
        self.show_latex()

class TypeA(CarterSolution):
    def Delta_lambda(self):
        lamda = self.chart[0]
        x = Lambda * lamda**4 / 3 + h * lamda**2 - 2 * m1 * lamda + n + e**2
        return x
    
    def Delta_mu(self):
        mu = self.chart[1]
        x = Lambda * mu**4 / 3 - (h * mu**2 - 2 * m2 * mu - n)
        return x

    def P_lambda(self):
        lamda = self.chart[0]
        return lamda**2

    def P_mu(self):
        mu = self.chart[1]
        return - mu**2  

    def Q_lambda(self):
        return 1
    
    def Q_mu(self):
        return 1

class TypeB_plus(CarterSolution):
    def Delta_lambda(self):
        lamda = self.chart[0]
        return Lambda * (lamda**4 / 3 + 2 * lamda**2 - 1) + (h * (lamda**2 - 1) - 2 * m1 * lamda + e**2)
    
    def Delta_mu(self):
        mu = self.chart[1]
        return - (h * (mu**2) - 2 * m2 * mu + n)
    
    def P_lambda(self):
        lamda = self.chart[0]
        return lamda**2 + 1
    
    def P_mu(self):
        mu = self.chart[1]
        return 0 # + 2 * mu

    def Q_lambda(self):
        return 0

    def Q_mu(self):
        return 1

class TypeC_plus(CarterSolution):
    def Delta_lambda(self):
        lamda = self.chart[0]
        return Lambda * (lamda**4) / 3 + (h * (lamda**2) - 2 * m1 * lamda + e**2)# 20250514 one third added
    
    def Delta_mu(self):
        mu = self.chart[1]
        return - (h * (mu**2) - 2 * m2 * mu + n)
    
    def P_lambda(self):
        lamda = self.chart[0]
        return lamda**2  # P_lambda = lambda^2
    
    def P_mu(self):
        return SR(0)  # P_mu = 0
    
    def Q_lambda(self):
        return SR(0)  # Q_lambda = 0
    
    def Q_mu(self):
        return SR(1)  # Q_mu = 1
    
class TypeD(CarterSolution):
        
    def Delta_lambda(self):
        lamda = self.chart[0]
        return (Lambda * (lamda^2) + 
                (+1) * (e^2 * (lamda^2) - 2 * m1 * lamda + n1))
    def Delta_mu(self):
        mu = self.chart[1]
        return (Lambda * (mu^2) + 
                (-1) * (e^2 * (mu^2) - 2 * m2 * mu + n2))
    
    def P_lambda(self):
        return SR(1/2)  # Convert to Sage symbolic representation
    
    def P_mu(self):
        # negative
        return SR(-1/2)  # Convert to Sage symbolic representation
    
    def Q_lambda(self):
        return 1
    
    def Q_mu(self):
        return 1

def show_derivative(solution):
    # Display the tetrads and their exterior derivatives
    tetrads = solution.omega_forms()
    for x in tetrads:
        print(x.display())
        show(x.exterior_derivative().display())

def calculate_energy_momentum(solution):
    potential = solution.elemag_potential()
    # print(potential.display())
    F = solution.field_strength_coord()
    # field_strength = solution.field_strength_tetrad()
    eta = diagonal_matrix([1, 1, 1, -1])
    four_pi_T = solution.energy_momentum_tensor(eta, eta)
    show(2 * four_pi_T)
    return four_pi_T

def calculate_einstein_tensor(solution, cosmological = True):
    t0 = time.time()
    frame = solution.compute_orthonormal_frame(chart=solution.chart)# show_commutator(frame)
    # g = solution.metric_tensor() # show(g.display())
    # omega = connection_one_forms(solution, gamma)
    eta = diagonal_matrix([1, 1, 1, -1]) # Minkowski metric
    g = eta
    c = solution.compute_structure_coefficients(frame)

    gamma = all_upper_coefficients(g, c, g_inv = "eta", domain = M)
    Riem = Riemannian_curvature(frame, gamma, c)
    t1= time.time(); print(t1 - t0);
    Ric = Ricci_tensor(Riem)
    # show(show_two_tensor(Ric, _dim))
    R = scalar_curvature(Ric, eta)
    t2= time.time(); print(t2 - t1);
    Einstein_tensor = matrix(Ric) - 1 / 2 * R * g
    # show_two_tensor(Einstein_tensor, _dim, array=False)
    if cosmological:
        LHS = Einstein_tensor - (Lambda * g)
        return LHS
    return Einstein_tensor


# symbolic = 'f'
# f_sym = M.scalar_field(function(symbolic)(*chart), name=symbolic)# symbolic function

# Create a row vector of vector fields
# part_lambda, part_mu, part_psi, part_chi = chart.frame()  # coordinate vector fields

# Define the 4-dimensional manifold with coordinates
_dim = 4
M = Manifold(_dim, 'M', structure='differentiable')
# Instantiate and use the class
# solution = TypeA(chart=chart)
# solution = TypeB_plus(chart=chart)
# solution = TypeC_plus(chart=chart)
solution = TypeD(M)

def main(solution):
    LHS = calculate_einstein_tensor(solution)
    RHS = 2 * calculate_energy_momentum(solution)
    eq = LHS - RHS
    # show_two_tensor(Einstein_tensor, _dim, array=False)
    # show_two_tensor(four_pi_T, _dim, array=False)
    show_two_tensor(eq, _dim, array=False, Mathematica=False)
    return 0

main(solution)
print(0)