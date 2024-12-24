# Function to define the metric tensor
def define_metric(M, Z, Delta_lamda, Delta_mu, l, mu):
    # Define the metric tensor
    g = M.metric()

    # Define the metric components
    g[0, 0] = Z / Delta_lamda                # g_lambda_lambda
    g[1, 1] = Z / Delta_mu                   # g_mu_mu
    g[2, 2] = -Delta_lamda / Z               # g_chi_chi
    g[3, 3] = (Z - 4 * l^2 * mu^2 * Delta_lamda) * Delta_mu  # g_psi_psi
    g[2, 3] = g[3, 2] = -2 * l * mu * Delta_lamda / Z  # g_chi_psi = g_psi_chi

    return g

# Example: Call the function with appropriate arguments

M = Manifold(4, 'M', structure='Lorentzian')
# Define coordinates
C.<lamda, mu, chi, psi> = M.chart()

Z, Delta_lamda, Delta_mu, l, mu = var('Z Delta_lamda Delta_mu l mu')
g = define_metric(M, Z, Delta_lamda, Delta_mu, l, mu)
print(g.display())

time ginv = g.inverse()
print(ginv.display())

print("========")
# 変数の定義
var('lamda h m e q p Lamda')

# 数式の定義
l = 0
Lamda = 0
Z = lamda^2 + l^2
Delta_lamda = Lamda * (1/3 * lamda^4 + 2 * l^2 * lamda^2 - l^4) \
              + h * (lamda^2 - l^2) - 2 * m * lamda + e^2
Delta_mu = -h * mu^2 + 2 * q * mu + p

# g_at_Lamda_zero
g = define_metric(M, Z, Delta_lamda, Delta_mu, l, mu)
print(g.display())

time ric = g.ricci()
print(ric.display())