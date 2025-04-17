# Define the 4-dimensional manifold with coordinates
M = Manifold(4, 'M', structure='differentiable')
chart = M.chart(r'lambda_ mu psi chi')
var('Lambda e h')
m1, m2 = var('m1 m2')
n1, n2 = var('n1 n2')

class CarterSolution:
    def __init__(self, chart):
        self.chart = chart
    
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
    
    def metric(self):
        """ Calculate the metric ds^2 """
        M = self.chart.domain()
        g = M.metric(name="g")  # metric tensor g_ij

        # 座標変数
        # lambda_, mu, psi, chi = self.chart._first, self.chart._second, self.chart._third, self.chart._fourth

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
    def hat_E(self):
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

    def compute_orthonormal_basis(self):
        """
        Computes and displays the transformed vector fields forming an orthonormal basis.
        
        Parameters:
            self.chart: The chart or manifold coordinate frame (with .frame() method)
            self = solution: An object with an .E_inv() method 
        Returns:
            List of vector fields forming an orthonormal basis.
        """
        orthonormal_basis = []
        for j in range(4):
            vec = sum(self.chart.frame()[i] * self.E_inv()[i, j] for i in range(4))
            orthonormal_basis.append(vec)
            show(vec.display())
        return orthonormal_basis

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
        return 2 * mu

    def Q_lambda(self):
        return 0

    def Q_mu(self):
        return 1

class TypeC_plus(CarterSolution):
    n = var('n')
        
    def Delta_lambda(self):
        lamda = self.chart[0]
        return Lambda * (lamda**4) + (h * (lamda**2) - 2 * m1 * lamda + e**2)
    
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
        return SR(1/2)  # Convert to Sage symbolic representation
    
    def Q_lambda(self):
        return 1
    
    def Q_mu(self):
        return 1

symbolic = 'f'
f_sym = M.scalar_field(function(symbolic)(*chart), name=symbolic)# symbolic function

def commutator_fields(basis, i, j, dummy_sym = f_sym):
    a = basis[i](basis[j](f_sym)) - basis[j](basis[i](f_sym))
    return a

def show_commutator(basis):
    n = len(basis)
    for i in range(n):
        for j in range(n):
            if i < j:
                show(i, j)
                x = commutator_fields(basis, i, j).display()
                show(x)
                show(x.coefficient(diff(f, mu)))

# Instantiate and use the class
# solution = TypeB_plus(chart=chart)
solution = TypeC_plus(chart=chart)
# solution = TypeD(chart=chart)

# Create a row vector of vector fields
part_lambda, part_mu, part_psi, part_chi = chart.frame()  # coordinate vector fields

orthonormal_basis = solution.compute_orthonormal_basis()
# show_commutator(orthonormal_basis)
x = commutator_fields(orthonormal_basis, 0, 1)
x_expr = x.expr()
coeff = x_expr.coefficient(diff(f_sym.expr(), chart[1]))
show(coeff)

# g = solution.metric()
# show(g.display())


# tetrads = solution.omega_forms()
# # Display the tetrads and their exterior derivatives
# for x in tetrads:
#     print(x.display())
#     show(x.exterior_derivative().display())
#     print(latex(diff(x).display()))