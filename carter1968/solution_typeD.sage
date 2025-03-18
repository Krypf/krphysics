# Define the 4-dimensional manifold with coordinates
M = Manifold(4, 'M', structure='differentiable')
chart = M.chart(r'lambda_ mu psi chi')
var('Lambda e h')
m1, m2 = var('m1 m2')
n1, n2 = var('n1 n2')

class typeD:
    def __init__(self, chart):
        self.chart = chart
        
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

    def show(self):
        # Print specific computed expressions
        print(f"Delta_lambda = {self.Delta_lambda()}")
        print(f"Delta_mu = {self.Delta_mu()}")
        print(f"P_lambda = {self.P_lambda()}")
        print(f"P_mu = {self.P_mu()}")
        print(f"Q_lambda = {self.Q_lambda()}")
        print(f"Q_mu = {self.Q_mu()}")

# Instantiate class
solution_typeD = typeD(chart=chart)

# Display results
solution_typeD.show()