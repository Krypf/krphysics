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

    def omega_forms(self):
        Z = self.Z()
        Delta_lambda = self.Delta_lambda()
        Delta_mu = self.Delta_mu()
        P_lambda = self.P_lambda()
        P_mu = self.P_mu()
        Q_lambda = self.Q_lambda()
        Q_mu = self.Q_mu()

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

    def show(self):
        # Print specific computed expressions
        print(f"Delta_lambda = {self.Delta_lambda()}")
        print(f"Delta_mu = {self.Delta_mu()}")
        print(f"P_lambda = {self.P_lambda()}")
        print(f"P_mu = {self.P_mu()}")
        print(f"Q_lambda = {self.Q_lambda()}")
        print(f"Q_mu = {self.Q_mu()}")

    def show_latex(self):
        # Print computed expressions in LaTeX format
        # from sage.misc.latex import latex
        print(f"Delta_lambda = {latex(self.Delta_lambda())}")
        print(f"Delta_mu = {latex(self.Delta_mu())}")
        print(f"P_lambda = {latex(self.P_lambda())}")
        print(f"P_mu = {latex(self.P_mu())}")
        print(f"Q_lambda = {latex(self.Q_lambda())}")
        print(f"Q_mu = {latex(self.Q_mu())}")

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

    
# Instantiate class
solution_TypeB_plus = TypeB_plus(chart=chart)
solution_TypeC_plus = TypeC_plus(chart=chart)
solution_TypeD = TypeD(chart=chart)

# Display results
# print("TypeB_plus")
# solution_TypeB_plus.show()
# solution_TypeB_plus.show_latex()
print("TypeC_plus")
solution_TypeC_plus.show()
solution_TypeC_plus.show_latex()
# print("TypeD")
# solution_TypeD.show()
# solution_TypeD.show_latex()

# Instantiate and use the class
solution = solution_TypeC_plus
tetrads = solution.omega_forms()

# Display the tetrads and their exterior derivatives
for x in tetrads:
    print(x.display())
    show(x.exterior_derivative().display())
    print(latex(diff(x).display()))