"""
https://claude.ai/chat/bddcbd78-06c6-4d81-8939-c947d8e4020fhttps://claude.ai/chat/bddcbd78-06c6-4d81-8939-c947d8e4020f
Kerr時空での自由落下プラズマ中の光線軌道（ray tracing）
変数分離不可能なケース → ハミルトン方程式を直接数値積分

プラズマ密度（あなたの論文式45）:
  n(r,θ) = C(θ)/sqrt(2Mr(r^2+a^2))
  ω_p^2 = A_p * ω_0^2 / sqrt(r*(r^2+a^2))  (C(θ)=const, M=1 単位系)
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import warnings
warnings.filterwarnings('ignore')

m  = 1.0
a  = 1.0 #0.9 * m   # 高スピン Kerr
A_p = 1.0      # プラズマ係数

r_h = m + np.sqrt(m**2 - a**2)   # 外側ホライズン

# ============================================================
# Kerr逆計量
# ============================================================
def kerr_inv(r, th):
    Sigma  = r**2 + a**2*np.cos(th)**2 #OK
    Delta  = r**2 + a**2 - 2*m*r #OK
    A_   = (r**2+a**2)**2 - a**2*Delta*np.sin(th)**2 #OK
    s2   = np.sin(th)**2 # sin theta squared
    gtt   = -A_/(Sigma*Delta) #OK
    gtphi = -2*a*m*r/(Sigma*Delta) #OK
    grr   =  Delta/Sigma #OK
    gthth =  1.0/Sigma #OK
    gphph = (Delta - a**2*s2)/(Sigma*Delta*s2) if abs(s2) > 1e-10 else 1e10
    return gtt, gtphi, grr, gthth, gphph, Sigma, Delta

# ============================================================
# 自由落下プラズマ
# ============================================================
def wp2(r, th):
    d = 2*m*r*(r**2 + a**2)
    return A_p / np.sqrt(d) if d > 0 else 0.0

# ============================================================
# ハミルトン方程式の右辺
# ============================================================
def make_rhs(use_plasma):
    def rhs(s, y):
        r, th, phi, pr, pth, pph = y
        if r < r_h * 1.005:
            return [0.0]*6
        gtt, gtph, grr, gthth, gphph, Sigma, Delta = kerr_inv(r, th)
        w2  = wp2(r, th) if use_plasma else 0.0
        pt  = -1.0  # ω0=1

        dr_ds   = grr   * pr
        dth_ds  = gthth * pth
        dphi_ds = gtph  * pt + gphph * pph

        # ∂H/∂r, ∂H/∂θ を数値微分
        eps_r  = max(r * 1e-5, 1e-8)
        eps_th = 1e-5

        def H_val(rv, thv):
            g1,g2,g3,g4,g5,_,_ = kerr_inv(rv, thv)
            w = wp2(rv, thv) if use_plasma else 0.0
            return 0.5*(g1*pt**2 + 2*g2*pt*pph +
                        g3*pr**2 + g4*pth**2 + g5*pph**2 + w)

        dHdr  = (H_val(r+eps_r, th) - H_val(r-eps_r, th)) / (2*eps_r)
        dHdth = (H_val(r, th+eps_th) - H_val(r, th-eps_th)) / (2*eps_th)

        return [dr_ds, dth_ds, dphi_ds, -dHdr, -dHdth, 0.0]
    return rhs

def ev_hor(s, y): return y[0] - r_h*1.005
ev_hor.terminal = True;  ev_hor.direction = -1

def ev_far(s, y): return y[0] - 55.0*m
ev_far.terminal = True;  ev_far.direction = +1

# ============================================================
# 初期条件生成（赤道面 θ=π/2 入射）
# ============================================================
r0 = 30.0*m
# the initial angle = pi
# てん回てん = R 
def make_ic(pph, use_plasma):
    gtt,gtph,grr,_,gphph,_,_ = kerr_inv(r0, np.pi/2)
    w2  = wp2(r0, np.pi/2) if use_plasma else 0.0
    pt  = -1.0
    rhs = -(gtt*pt**2 + 2*gtph*pt*pph + gphph*pph**2 + w2) / grr
    if rhs <= 0:
        return None
    pr0 = -np.sqrt(rhs)
    return [r0, np.pi/2, np.pi, pr0, 0.0, pph]

def trace_ray(pph, use_plasma):
    y0 = make_ic(pph, use_plasma)
    if y0 is None:
        return None, None
    rhs = make_rhs(use_plasma)
    sol = solve_ivp(rhs, (0, 10000), y0,
                    method='DOP853',
                    events=[ev_hor, ev_far],
                    rtol=1e-9, atol=1e-11, max_step=0.4)
    x = sol.y[0] * np.cos(sol.y[2])
    y = sol.y[0] * np.sin(sol.y[2])
    return x, y

# ============================================================
# impact parameters
# ============================================================
pph_values = [2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 12.0, 14.0, 16.0]

cmap   = plt.cm.plasma
colors = cmap(np.linspace(0.1, 0.9, len(pph_values)))

lim = 22

fig, axes = plt.subplots(1, 2, figsize=(14, 7), facecolor='white')
plt.subplots_adjust(wspace=0.35)

for ax, use_plasma, title in zip(
    axes,
    [False, True],
    [r'Vacuum  ($\omega_p=0$)',
     r'Free-fall plasma  ($\omega_p^2 \propto [r(r^2+a^2)]^{-1/2}$)']
):
    ax.set_facecolor('white')
    ax.grid(True, color='lightgray', lw=0.5, zorder=0)
    ax.set_axisbelow(True)

    for pph, col in zip(pph_values, colors):
        x_arr, y_arr = trace_ray(pph, use_plasma)
        if x_arr is None:
            continue
        ax.plot(x_arr, y_arr, color=col, lw=1.0, alpha=0.9, zorder=3)

    # ブラックホール（黒塗り）
    bh = plt.Circle((0,0), r_h, color='black', zorder=6)
    ax.add_patch(bh)

    # エルゴ領域境界（赤道面: r_ergo=2m）
    r_ergo = 2*m
    ergo = plt.Circle((0,0), r_ergo, color='dimgray', fill=False,
                       lw=0.9, ls='--', alpha=0.8, zorder=5)
    ax.add_patch(ergo)

    # 真空光子球（参考）
    r_ph_vac = m*(1 + 2*np.cos(np.arccos(-a/m)/3))  # Kerr光子球（赤道面）
    ps = plt.Circle((0,0), r_ph_vac, color='steelblue', fill=False,
                     lw=0.9, ls=':', alpha=0.8, zorder=5)
    ax.add_patch(ps)

    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect('equal')
    ax.set_xlabel(r'$x\,/\,M$', fontsize=12)
    ax.set_ylabel(r'$y\,/\,M$', fontsize=12)
    ax.set_title(title, fontsize=11, pad=8)
    ax.set_xticks(np.arange(-20, 21, 5))
    ax.set_yticks(np.arange(-20, 21, 5))
    for sp in ax.spines.values(): sp.set_color('black')
    ax.tick_params(colors='black', labelsize=9)

    # 凡例
    ax.text(r_ergo*np.cos(np.radians(50)),
            r_ergo*np.sin(np.radians(50))+0.5,
            r'$r_{\rm ergo}$', color='dimgray', fontsize=8, ha='center')
    ax.text(r_ph_vac*np.cos(np.radians(130)),
            r_ph_vac*np.sin(np.radians(130))+0.5,
            r'$r_{\rm ph}$', color='steelblue', fontsize=8, ha='center')

# カラーバー
sm = plt.cm.ScalarMappable(cmap=cmap,
     norm=plt.Normalize(vmin=min(pph_values), vmax=max(pph_values)))
sm.set_array([])
cbar = fig.colorbar(sm, ax=axes, orientation='vertical',
                    fraction=0.022, pad=0.02, shrink=0.85)
cbar.set_label(r'$p_\varphi\,/\,(\omega_0 M)$', fontsize=10)
cbar.ax.tick_params(labelsize=9)

fig.suptitle(
    f'Ray Tracing in Kerr Spacetime ($a={a/m:g}M$)  —  '
    r'Free-fall plasma: $n(r,\theta)\propto[r(r^2+a^2)]^{-1/2}$',
    fontsize=12, y=1.01, color='black')

plt.savefig('docs/kerr_freefall_raytracing.png',
            dpi=160, bbox_inches='tight', facecolor='white')
plt.close()
print("Done.")
