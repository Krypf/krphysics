# Define the 4-dimensional manifold with coordinates
M = Manifold(4, 'M', structure='differentiable')
chart = M.chart(r'lambda_ mu psi chi')

# Define the scalar functions
# Z, Delta_lambda, Delta_mu, P_lambda, Q_lambda, P_mu, Q_mu = var('Z Delta_lambda Delta_mu P_lambda Q_lambda P_mu Q_mu')
Delta_lambda = M.scalar_field(function('Delta_lambda')(chart[0]))
Delta_mu = M.scalar_field(function('Delta_mu')(chart[1]))
P_lambda = M.scalar_field(function('P_lambda')(chart[0]))
Q_lambda = M.scalar_field(function('Q_lambda')(chart[0]))
P_mu = M.scalar_field(function('P_mu')(chart[1]))
Q_mu = M.scalar_field(function('Q_mu')(chart[1]))
# Define Z
Z = P_lambda * Q_mu - P_mu * Q_lambda

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

# Display the forms
tetrads = [omega_p1, omega_m1, omega_p2, omega_m2]
for x in tetrads:
    print(x.display())
    show(x.exterior_derivative().display())
    print(latex(diff(x).display()))