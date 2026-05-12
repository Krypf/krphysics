"""
弱偏向角の解析式 vs 厳密数値積分（修正版）

修正点:
  - 積分の上限∞を u=R/r 置換で有限区間 [0,1] に変換
  - 弱偏向式は b ではなく R の関数として正しく評価（式117）
"""

import numpy as np
from scipy import integrate, special
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

M = 1.0
k = 1.5
C = 7.0

# ============================================================
# 基本関数
# ============================================================
def f(r):
    return 1.0 - 2*M/r

def h2(r, C_val=C, k_val=k):
    """Tsupko (2021) 式(22): h^2 = r^2/f * (1 - f*C*(M/r)^k)"""
    fr = f(r)
    return r**2/fr * (1.0 - fr*C_val*(M/r)**k_val)

def b_from_R(R, C_val=C, k_val=k):
    """b = h(R) (Tsupko eq.23)"""
    h2R = h2(R, C_val, k_val)
    return np.sqrt(h2R) if h2R > 0 else np.nan

# ============================================================
# 厳密積分（u=R/r 置換で安定化）
# ============================================================
def integrand_u(u, R, C_val=C, k_val=k):
    """u = R/r 置換後の被積分関数"""
    if u <= 0 or u > 1: return 0.0
    r    = R / u
    h2r  = h2(r, C_val, k_val)
    h2R  = h2(R, C_val, k_val)
    fr   = f(r)
    ratio = h2r/h2R - 1.0
    if ratio <= 0 or fr <= 0 or h2R <= 0: return 0.0
    # |dr/du| = R/u^2
    return R / (r * np.sqrt(fr) * np.sqrt(ratio) * u**2)

def alpha_exact(R, C_val=C, k_val=k):
    """Tsupko 式(26): alpha = 2*int_R^inf ... dr - pi"""
    h2R = h2(R, C_val, k_val)
    if h2R <= 0: return np.nan
    try:
        res, err = integrate.quad(
            integrand_u, 0, 1,
            args=(R, C_val, k_val),
            limit=2000, epsabs=1e-10, epsrel=1e-10
        )
        return 2.0*res - np.pi
    except:
        return np.nan

# ============================================================
# 弱偏向式（Perlick & Tsupko 2024 式116・117）
# ============================================================
def integrand_refr(r, R, C_val=C, k_val=k):
    """式(116)の被積分関数: dN/dr / sqrt(r^2-R^2)"""
    # N(r) = C*(M/r)^k  →  dN/dr = -k*C*M^k / r^{k+1}
    dNdr = -k_val * C_val * M**k_val / r**(k_val+1)
    denom = np.sqrt(r**2 - R**2)
    if denom <= 0: return 0.0
    return dNdr / denom

def alpha_weak_R(R, C_val=C, k_val=k):
    """
    式(117): alpha(R) = 4M/R + alpha_refr(R)
    式(116): alpha_refr(R) = R * int_R^inf dN/dr/sqrt(r^2-R^2) dr
    """
    alpha_grav = 4*M / R
    # 同じく u=R/r 置換で安定化
    def integrand_refr_u(u, R_, C_val_, k_val_):
        if u <= 0 or u >= 1: return 0.0
        r = R_/u
        dNdr = -k_val_*C_val_*M**k_val_ / r**(k_val_+1)
        denom = np.sqrt(r**2 - R_**2)
        if denom <= 0: return 0.0
        return dNdr / denom * R_/u**2
    res, _ = integrate.quad(
        integrand_refr_u, 0, 1,
        args=(R, C_val, k_val),
        limit=1000, epsabs=1e-10, epsrel=1e-10
    )
    alpha_refr = R * res
    return alpha_grav, alpha_refr, alpha_grav + alpha_refr

# ============================================================
# 光子球の位置
# ============================================================
from scipy.optimize import brentq

def dh2_dr(r, C_val=C, k_val=k):
    eps = r*1e-6
    return (h2(r+eps, C_val, k_val) - h2(r-eps, C_val, k_val))/(2*eps)

rph_vac = 3.0*M
b_ph_vac = b_from_R(rph_vac, 0.0)

try:
    rph_pl = brentq(dh2_dr, 3.01*M, 6.0*M)
    b_ph_pl = b_from_R(rph_pl)
except:
    rph_pl = 3.5*M
    b_ph_pl = b_from_R(rph_pl)

print(f"真空光子球:    r_ph={rph_vac:.4f}M,  b_crit={b_ph_vac:.4f}M")
print(f"プラズマ光子球: r_ph={rph_pl:.4f}M,  b_crit={b_ph_pl:.4f}M")

# ============================================================
# 数値計算
# ============================================================
print("\n=== R を引数とした比較 ===")
print("R/M    b/M      exact        weak(R)      grav         refr         diff%")

R_test = [10., 15., 20., 30., 50., 80., 120.]
for R_val in R_test:
    ae = alpha_exact(R_val)
    ag, ar, aw = alpha_weak_R(R_val)
    b_val = b_from_R(R_val)
    diff = abs(aw-ae)/abs(ae)*100 if abs(ae)>1e-8 else float('nan')
    print(f"R={R_val:5.1f}  b={b_val:7.4f}  "
          f"exact={ae:.6f}  weak={aw:.6f}  "
          f"grav={ag:.6f}  refr={ar:.6f}  diff={diff:.2f}%")

# ============================================================
# プロット
# ============================================================
R_arr = np.concatenate([
    np.linspace(rph_pl*1.01, 15*M, 150),
    np.linspace(15*M, 120*M, 100)
])

res_exact, res_weak, res_grav, res_refr, res_b = [], [], [], [], []

for R in R_arr:
    ae = alpha_exact(R)
    ag, ar, aw = alpha_weak_R(R)
    b  = b_from_R(R)
    if np.isnan(ae) or np.isnan(aw): continue
    res_exact.append(ae); res_weak.append(aw)
    res_grav.append(ag);  res_refr.append(ar)
    res_b.append(b)

res_exact = np.array(res_exact); res_weak = np.array(res_weak)
res_grav  = np.array(res_grav);  res_refr = np.array(res_refr)
res_b     = np.array(res_b)

# 真空の厳密解
R_vac = np.concatenate([np.linspace(3.01*M,20*M,150), np.linspace(20*M,120*M,80)])
ae_vac_list, b_vac_list = [], []
for Rv in R_vac:
    av = alpha_exact(Rv, 0.0)
    bv = b_from_R(Rv, 0.0)
    if not np.isnan(av):
        ae_vac_list.append(av); b_vac_list.append(bv)

fig, axes = plt.subplots(2, 1, figsize=(9, 10), facecolor='white')

# --- 上段: 偏向角 ---
ax = axes[0]
ax.set_facecolor('white')
ax.grid(True, color='lightgray', lw=0.5)
ax.plot(res_b, np.degrees(res_grav),  'b--', lw=1.2, label=r'gravity: $4M/R$')
ax.plot(res_b, np.degrees(res_refr),  'r--', lw=1.2, label=r'refraction: $\alpha_{\rm refr}(R)$')
ax.plot(res_b, np.degrees(res_weak),  'g-',  lw=2.0,
        label=r'weak approx: $\alpha_{\rm grav}+\alpha_{\rm refr}$ [eq.117]')
ax.plot(res_b, np.degrees(res_exact), 'k-',  lw=2.0,
        label=r'exact [Tsupko 2021 eq.26]')
ax.plot(b_vac_list, np.degrees(ae_vac_list), color='gray', lw=1.2, ls=':',
        label='vacuum exact')
ax.axhline(0, color='black', lw=0.8)
ax.set_xlabel(r'Impact parameter $b\,/\,M$', fontsize=12)
ax.set_ylabel(r'Deflection angle $\hat{\alpha}$ [deg]', fontsize=12)
ax.set_title(r'Free-fall plasma: $\omega_p^2=7\omega_0^2(M/r)^{3/2}$', fontsize=11)
ax.legend(fontsize=9, loc='upper right')
ax.set_xlim(0, 50); ax.set_ylim(-2, 20)
for sp_ in ax.spines.values(): sp_.set_color('black')
ax.tick_params(colors='black')

# --- 下段: 相対誤差 ---
ax2 = axes[1]
ax2.set_facecolor('white')
ax2.grid(True, color='lightgray', lw=0.5)
mask = np.abs(res_exact) > 1e-6
rel_err = np.where(mask, np.abs(res_weak-res_exact)/np.abs(res_exact)*100, np.nan)
ax2.semilogy(res_b[mask], rel_err[mask], 'purple', lw=1.5,
             label='|weak - exact| / |exact| [%]')
ax2.axhline(1.0,  color='gray', lw=0.8, ls='--', label='1%')
ax2.axhline(10.0, color='gray', lw=0.8, ls=':',  label='10%')
ax2.set_xlabel(r'Impact parameter $b\,/\,M$', fontsize=12)
ax2.set_ylabel('Relative error [%]', fontsize=12)
ax2.set_title('Accuracy of weak deflection approximation', fontsize=11)
ax2.legend(fontsize=9)
ax2.set_xlim(0, 50)
for sp_ in ax2.spines.values(): sp_.set_color('black')
ax2.tick_params(colors='black')

plt.tight_layout()
plt.savefig('tests/weak_deflection_corrected.png',
            dpi=160, bbox_inches='tight', facecolor='white')
plt.close()
print("\nDone.")
