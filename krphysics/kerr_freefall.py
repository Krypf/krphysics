"""
Kerr時空での自由落下プラズマ中の光線軌道（ray tracing）
========================================================

変数分離不可能なケース → ハミルトン方程式を直接数値積分する。

プラズマ密度（論文式45）::

    n(r,θ) = C(θ) / sqrt(2 M r (r^2 + a^2))
    ω_p^2 = A_p * ω_0^2 / sqrt(r (r^2 + a^2))     (C(θ)=const, M=1 単位系)
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import Callable, Sequence

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.patches import Circle
from scipy.integrate import solve_ivp

warnings.filterwarnings("ignore")


# =============================================================================
# 設定
# =============================================================================
@dataclass(frozen=True)
class KerrConfig:
    """ブラックホール・プラズマ・初期条件の物理パラメータをまとめる。"""

    # 質量とスピン（M=1 単位系）
    M: float = 1.0
    a: float = 1.0

    # 自由落下プラズマ係数
    A_p: float = 1.0

    # 入射初期半径（赤道面 θ=π/2 入射、進入方向 φ=π）
    r0: float = 30.0

    # 数値積分のホライズン安全係数（r_h * factor を打ち切り半径とする）
    horizon_safety: float = 1.005

    # 遠方の打ち切り半径（M単位）
    r_far: float = 55.0

    # 角運動量サンプル（impact parameter, p_φ / (ω_0 M)）
    pph_values: tuple[float, ...] = (
        2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 12.0, 14.0, 16.0,
    )

    @property
    def r_horizon(self) -> float:
        """外側ホライズン半径 r_+."""
        return self.M + np.sqrt(self.M**2 - self.a**2)

    @property
    def r_ergo_eq(self) -> float:
        """赤道面のエルゴ面（r_ergo = 2M）。"""
        return 2.0 * self.M

    @property
    def r_photon_eq(self) -> float:
        """赤道面 Kerr 光子球（解析表現の一つ）。"""
        return self.M * (1.0 + 2.0 * np.cos(np.arccos(-self.a / self.M) / 3.0))


# =============================================================================
# Kerr 計量（Boyer–Lindquist の逆計量）
# =============================================================================
@dataclass(frozen=True)
class KerrMetric:
    cfg: KerrConfig

    def inverse(self, r: float, th: float):
        """逆計量成分 (g^tt, g^tφ, g^rr, g^θθ, g^φφ) と Σ, Δ を返す。"""
        M, a = self.cfg.M, self.cfg.a
        s2 = np.sin(th) ** 2
        Sigma = r**2 + a**2 * np.cos(th) ** 2
        Delta = r**2 + a**2 - 2.0 * M * r
        A_ = (r**2 + a**2) ** 2 - a**2 * Delta * s2

        gtt = -A_ / (Sigma * Delta)
        gtph = -2.0 * a * M * r / (Sigma * Delta)
        grr = Delta / Sigma
        gthth = 1.0 / Sigma
        if abs(s2) > 1e-10:
            gphph = (Delta - a**2 * s2) / (Sigma * Delta * s2)
        else:
            gphph = 1e10  # 極軸近傍の発散回避用ガード
        return gtt, gtph, grr, gthth, gphph, Sigma, Delta


# =============================================================================
# 自由落下プラズマ
# =============================================================================
@dataclass(frozen=True)
class FreeFallPlasma:
    cfg: KerrConfig

    def omega_p_squared(self, r: float, th: float) -> float:
        """ω_p^2(r, θ)。式 (45) を M=1 単位系で評価。"""
        M, a, A_p = self.cfg.M, self.cfg.a, self.cfg.A_p
        d = 2.0 * M * r * (r**2 + a**2)
        return A_p / np.sqrt(d) if d > 0 else 0.0


# =============================================================================
# ハミルトン方程式の右辺
# =============================================================================
def make_hamiltonian_rhs(
    metric: KerrMetric,
    plasma: FreeFallPlasma,
    use_plasma: bool,
) -> Callable[[float, Sequence[float]], list[float]]:
    """光子用ハミルトニアン H = ½(g^μν p_μ p_ν + ω_p^2) の右辺を生成する。

    保存量 (p_t = -ω_0 = -1, p_φ = const) を用い、
    ∂H/∂r と ∂H/∂θ は中心差分で評価する。
    """
    cfg = metric.cfg
    r_cut = cfg.r_horizon * cfg.horizon_safety
    pt = -1.0  # ω_0 = 1

    def hamiltonian(r: float, th: float, pr: float, pth: float, pph: float) -> float:
        gtt, gtph, grr, gthth, gphph, _, _ = metric.inverse(r, th)
        w2 = plasma.omega_p_squared(r, th) if use_plasma else 0.0
        return 0.5 * (
            gtt * pt**2
            + 2.0 * gtph * pt * pph
            + grr * pr**2
            + gthth * pth**2
            + gphph * pph**2
            + w2
        )

    def rhs(_s: float, y: Sequence[float]) -> list[float]:
        r, th, _phi, pr, pth, pph = y
        if r < r_cut:
            return [0.0] * 6

        gtt, gtph, grr, gthth, gphph, _, _ = metric.inverse(r, th)

        # 配位変数の発展
        dr_ds = grr * pr
        dth_ds = gthth * pth
        dphi_ds = gtph * pt + gphph * pph

        # 共役運動量の発展（∂H/∂r, ∂H/∂θ を中心差分で評価）
        eps_r = max(r * 1e-5, 1e-8)
        eps_th = 1e-5
        dHdr = (
            hamiltonian(r + eps_r, th, pr, pth, pph)
            - hamiltonian(r - eps_r, th, pr, pth, pph)
        ) / (2.0 * eps_r)
        dHdth = (
            hamiltonian(r, th + eps_th, pr, pth, pph)
            - hamiltonian(r, th - eps_th, pr, pth, pph)
        ) / (2.0 * eps_th)

        return [dr_ds, dth_ds, dphi_ds, -dHdr, -dHdth, 0.0]

    return rhs


# =============================================================================
# レイトレーサ
# =============================================================================
class RayTracer:
    """赤道面（θ=π/2）入射の光線追跡。"""

    def __init__(self, cfg: KerrConfig):
        self.cfg = cfg
        self.metric = KerrMetric(cfg)
        self.plasma = FreeFallPlasma(cfg)

    # --- 終端イベント ---------------------------------------------------------
    def _event_horizon(self):
        cfg = self.cfg

        def _ev(_s, y):
            return y[0] - cfg.r_horizon * cfg.horizon_safety

        _ev.terminal = True
        _ev.direction = -1
        return _ev

    def _event_far(self):
        cfg = self.cfg

        def _ev(_s, y):
            return y[0] - cfg.r_far * cfg.M

        _ev.terminal = True
        _ev.direction = +1
        return _ev

    # --- 初期条件 -------------------------------------------------------------
    def initial_state(self, pph: float, use_plasma: bool):
        """ヌル条件 H = 0 を解いて p_r(初期) を決定する。"""
        cfg = self.cfg
        gtt, gtph, grr, _, gphph, _, _ = self.metric.inverse(cfg.r0, np.pi / 2)
        w2 = self.plasma.omega_p_squared(cfg.r0, np.pi / 2) if use_plasma else 0.0
        pt = -1.0
        rhs = -(gtt * pt**2 + 2.0 * gtph * pt * pph + gphph * pph**2 + w2) / grr
        if rhs <= 0:
            return None
        pr0 = -np.sqrt(rhs)  # 内向き
        return [cfg.r0, np.pi / 2, np.pi, pr0, 0.0, pph]

    # --- 軌道計算 -------------------------------------------------------------
    def trace(self, pph: float, use_plasma: bool):
        y0 = self.initial_state(pph, use_plasma)
        if y0 is None:
            return None, None
        rhs = make_hamiltonian_rhs(self.metric, self.plasma, use_plasma)
        sol = solve_ivp(
            rhs,
            (0.0, 10000.0),
            y0,
            method="DOP853",
            events=[self._event_horizon(), self._event_far()],
            rtol=1e-9,
            atol=1e-11,
            max_step=0.4,
        )
        r_arr, phi_arr = sol.y[0], sol.y[2]
        x = r_arr * np.cos(phi_arr)
        y = r_arr * np.sin(phi_arr)
        return x, y


# =============================================================================
# プロッタ
# =============================================================================
@dataclass
class PlotConfig:
    lim: float = 22.0
    figsize: tuple[float, float] = (14.0, 7.0)
    cmap_name: str = "plasma"
    dpi: int = 160


class RayPlotter:
    """真空 / 自由落下プラズマの2パネル比較プロット。"""

    def __init__(self, cfg: KerrConfig, plot_cfg: PlotConfig | None = None):
        self.cfg = cfg
        self.plot_cfg = plot_cfg or PlotConfig()
        self.tracer = RayTracer(cfg)

    def _draw_panel(self, ax, use_plasma: bool, title: str, colors):
        cfg = self.cfg
        plot_cfg = self.plot_cfg

        ax.set_facecolor("white")
        ax.grid(True, color="lightgray", lw=0.5, zorder=0)
        ax.set_axisbelow(True)

        # 光線
        for pph, col in zip(cfg.pph_values, colors):
            x, y = self.tracer.trace(pph, use_plasma)
            if x is None:
                continue
            ax.plot(x, y, color=col, lw=1.0, alpha=0.9, zorder=3)

        # 構造物（ホライズン・エルゴ面・光子球）
        ax.add_patch(Circle((0, 0), cfg.r_horizon, color="black", zorder=6))
        ax.add_patch(
            Circle(
                (0, 0), cfg.r_ergo_eq,
                color="dimgray", fill=False, lw=0.9, ls="--", alpha=0.8, zorder=5,
            )
        )
        ax.add_patch(
            Circle(
                (0, 0), cfg.r_photon_eq,
                color="steelblue", fill=False, lw=0.9, ls=":", alpha=0.8, zorder=5,
            )
        )

        # 軸装飾
        ax.set_xlim(-plot_cfg.lim, plot_cfg.lim)
        ax.set_ylim(-plot_cfg.lim, plot_cfg.lim)
        ax.set_aspect("equal")
        ax.set_xlabel(r"$x\,/\,M$", fontsize=12)
        ax.set_ylabel(r"$y\,/\,M$", fontsize=12)
        ax.set_title(title, fontsize=11, pad=8)
        ax.set_xticks(np.arange(-20, 21, 5))
        ax.set_yticks(np.arange(-20, 21, 5))
        for sp in ax.spines.values():
            sp.set_color("black")
        ax.tick_params(colors="black", labelsize=9)

        # ラベル
        ax.text(
            cfg.r_ergo_eq * np.cos(np.radians(50)),
            cfg.r_ergo_eq * np.sin(np.radians(50)) + 0.5,
            r"$r_{\rm ergo}$",
            color="dimgray", fontsize=8, ha="center",
        )
        ax.text(
            cfg.r_photon_eq * np.cos(np.radians(130)),
            cfg.r_photon_eq * np.sin(np.radians(130)) + 0.5,
            r"$r_{\rm ph}$",
            color="steelblue", fontsize=8, ha="center",
        )

    def render(self, out_path: str) -> None:
        cfg = self.cfg
        plot_cfg = self.plot_cfg

        cmap = plt.get_cmap(plot_cfg.cmap_name)
        colors = cmap(np.linspace(0.1, 0.9, len(cfg.pph_values)))

        fig, axes = plt.subplots(1, 2, figsize=plot_cfg.figsize, facecolor="white")
        plt.subplots_adjust(wspace=0.35)

        panels = [
            (False, r"Vacuum  ($\omega_p=0$)"),
            (True, r"Free-fall plasma  ($\omega_p^2 \propto [r(r^2+a^2)]^{-1/2}$)"),
        ]
        for ax, (use_plasma, title) in zip(axes, panels):
            self._draw_panel(ax, use_plasma, title, colors)

        # カラーバー
        sm = ScalarMappable(
            cmap=cmap,
            norm=Normalize(vmin=min(cfg.pph_values), vmax=max(cfg.pph_values)),
        )
        sm.set_array([])
        cbar = fig.colorbar(
            sm, ax=axes, orientation="vertical",
            fraction=0.022, pad=0.02, shrink=0.85,
        )
        cbar.set_label(r"$p_\varphi\,/\,(\omega_0 M)$", fontsize=10)
        cbar.ax.tick_params(labelsize=9)

        fig.suptitle(
            f"Ray Tracing in Kerr Spacetime ($a={cfg.a / cfg.M:g}M$)  —  "
            r"Free-fall plasma: $n(r,\theta)\propto[r(r^2+a^2)]^{-1/2}$",
            fontsize=12, y=1.01, color="black",
        )

        plt.savefig(out_path, dpi=plot_cfg.dpi, bbox_inches="tight", facecolor="white")
        plt.close(fig)


# =============================================================================
# エントリポイント
# =============================================================================
def main() -> None:
    cfg = KerrConfig()
    plotter = RayPlotter(cfg)
    out_path = "kerr_freefall_raytracing.png"
    plotter.render(out_path)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
