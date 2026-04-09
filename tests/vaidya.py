import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def solve_vaidya_photon_surface(u1_list, lam):
    # Equation 103 setup
    def dydt(u, y):
        r, drdu = y
        if r < 1e-6: return [drdu, 0]
        m = lam * u
        dmdu = lam
        term1 = 1 - 3 * m / r
        term2 = 1 - 2 * m / r - 3 * drdu
        inside_braces = term1 * term2 - dmdu + 2 * (drdu**2)
        d2rdu2 = (1/r) * inside_braces
        return [drdu, d2rdu2]

    plt.figure(figsize=(10, 6))
    
    # Plotting loop
    for u1 in u1_list:
        m1 = lam * u1
        y0 = [3 * m1, 0] # Boundary condition Eq 105
        t_span = (u1, 1e-4) # Backward evolution
        
        sol = solve_ivp(dydt, t_span, y0, max_step=0.01, rtol=1e-6, atol=1e-8)
        
        u_vals = sol.t
        r_vals = sol.y[0]
        
        # Plot r/m1
        plt.plot(u_vals, r_vals / m1, label=f'$u_1 = {u1}$')

    # Add reference line
    plt.axhline(y=3, color='gray', linestyle='--', alpha=0.5, label=r'$r = 3m_1$')
    
    # Formatting
    plt.xlabel('$u$')
    plt.ylabel('$r/m_1$')
    plt.title(f'Backwards evolution of photon surface ($\lambda = {lam}$)')
    plt.legend()
    plt.grid(True)
    plt.xlim(0, max(u1_list) + 1)
    
    # Set Y-axis limit as requested
    plt.ylim(2.5, 4.0) # Lower limit set to 2.5 to capture the dip, upper to 4 as requested
    
    plt.savefig('vaidya_photon_surface_clipped.png')

# Parameters
lambda_val = 1.0 / 16.0
u1_values = [1, 5, 10]

solve_vaidya_photon_surface(u1_values, lambda_val)