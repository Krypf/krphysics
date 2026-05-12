import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

rs = 1.0

def f(r): return 1 - rs/r
def df(r): return rs / r**2

# ハミルトン方程式に基づくODE系
def hamiltonian_odes(lam, y):
    r, phi, pr, pphi, pt = y
    
    r_dot = f(r) * pr
    phi_dot = pphi / r**2
    pr_dot = -0.5 * ( (pt**2 * df(r)) / (f(r)**2) + df(r) * pr**2 - 2 * pphi**2 / r**3 )
    
    pphi_dot = 0
    pt_dot = 0
    
    return [r_dot, phi_dot, pr_dot, pphi_dot, pt_dot]

def event_horizon(lam, y):
    return y[0] - 1.05 * rs
event_horizon.terminal = True
event_horizon.direction = -1

def event_escape(lam, y):
    return y[0] - 45.0
event_escape.terminal = True
event_escape.direction = 1

bc = 3 * np.sqrt(3) / 2 * rs
b_values = [2.0, bc, 2.8, 3.5, 5.2]
colors = ['purple', 'red', 'orange', 'green', 'blue']

plt.figure(figsize=(10, 10))

for b, color in zip(b_values, colors):
    E = 1.0
    L = b * E
    
    # 初期条件 (無限遠に近い場所からスタート)
    r0 = 40.0
    phi0 = 0.0
    pt0 = -E
    pphi0 = L
    
    # ハミルトニアン H = 0 から初期の p_r を決定
    val = (pt0**2) / (f(r0)**2) - (pphi0**2) / (r0**2 * f(r0))
    pr0 = -np.sqrt(max(val, 0)) # マイナスは内向き(ブラックホールに向かう)を意味する
    
    y0 = [r0, phi0, pr0, pphi0, pt0]
    
    # 数値積分
    sol = solve_ivp(hamiltonian_odes, [0, 150], y0, 
                    events=[event_horizon, event_escape], 
                    max_step=0.1, rtol=1e-8, atol=1e-10)
    
    r_sol = sol.y[0]
    phi_sol = sol.y[1]
    
    # 描画用の座標変換 (右から左へ入射するように回転)
    x = r_sol * np.cos(np.pi - phi_sol)
    y = r_sol * np.sin(np.pi - phi_sol)
    
    label = f'b = $b_c$' if b == bc else f'b = {b}'
    plt.plot(x, y, color=color, label=label)

# ブラックホールと光子球の描画
circle_bh = plt.Circle((0, 0), rs, color='black')
plt.gca().add_patch(circle_bh)
circle_ph = plt.Circle((0, 0), 1.5*rs, color='red', fill=False, linestyle='--', label='Photon Sphere ($1.5r_s$)')
plt.gca().add_patch(circle_ph)

plt.xlim(-10, 10)
plt.ylim(-10, 10)
plt.axhline(0, color='gray', linestyle='-', linewidth=0.5)
plt.axvline(0, color='gray', linestyle='-', linewidth=0.5)
plt.xlabel('x ($r_s$)')
plt.ylabel('y ($r_s$)')
plt.title('Null Geodesics computed via Hamiltonian Formalism')
plt.legend()
plt.grid(alpha=0.3)
plt.gca().set_aspect('equal')
plt.savefig('tests/hamiltonian_lensing.png', dpi=300)
print("Image saved")