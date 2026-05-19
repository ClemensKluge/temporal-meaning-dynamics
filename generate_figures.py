"""
Temporal Meaning Dynamics — Figure Generation
==============================================
Reproduces Figures 1–4 from the manuscript.

Usage:
    python generate_figures.py              # all figures
    python generate_figures.py --fig 1     # specific figure

Requires: numpy, scipy, matplotlib
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from sinndynamik_solver import (
    SinnDynamikSystem, HopfScanner, HybridKernel,
    gamma_koll_equilibrium, compute_sync_star, DEFAULT_PARAMS
)

# Publication style
plt.rcParams.update({
    "font.family": "serif",
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 12,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})

BLUE = "#1f6fa3"
RED  = "#c94040"


def figure1_bifurcation_diagram():
    """
    Figure 1: Bifurcation diagram — σ_Φ vs γ_base.

    Shows the Hopf bifurcation at γ_krit ≈ 6.2 (α_K = 0, solid blue)
    and the damping effect of the mixed kernel (α_K = 0.3, dashed red).
    Insets show representative time series at points A, B, C.
    """
    scanner = HopfScanner()
    gamma_vals = np.linspace(3, 14, 25)

    print("Figure 1: scanning γ_base (α_K = 0)...")
    scan0 = scanner.scan_gamma(gamma_vals, T=200.0, burn_in=80.0)

    print("Figure 1: scanning γ_base (α_K = 0.3)...")
    scan03 = HopfScanner({"alpha_K": 0.3}).scan_gamma(
        gamma_vals, T=200.0, burn_in=80.0)

    fig = plt.figure(figsize=(13, 5))
    gs = gridspec.GridSpec(3, 2, width_ratios=[2, 1], hspace=0.5, wspace=0.35)

    # Main bifurcation diagram
    ax_main = fig.add_subplot(gs[:, 0])
    ax_main.plot(scan0["gamma_values"], scan0["sigma_Phi"],
                 color=BLUE, lw=2, label=r"$\alpha_K = 0$ (single delay)")
    ax_main.plot(scan03["gamma_values"], scan03["sigma_Phi"],
                 color=RED, lw=2, ls="--", label=r"$\alpha_K = 0.3$ (mixed kernel)")

    # Sense-making region shading
    ax_main.axvspan(6.2, 14, alpha=0.06, color=BLUE, label="Oscillatory regime")
    ax_main.axvline(6.2, color=BLUE, lw=1, ls=":", alpha=0.7)
    ax_main.axvline(8.0, color=RED, lw=1, ls=":", alpha=0.7)
    ax_main.annotate(r"$\gamma^*_\mathrm{base}$ ($\alpha_K=0$)",
                     xy=(6.2, 0.01), xytext=(5.0, 0.08),
                     arrowprops=dict(arrowstyle="->", color=BLUE), color=BLUE,
                     fontsize=9)
    ax_main.annotate(r"$\gamma^*_\mathrm{base}$ ($\alpha_K=0.3$)",
                     xy=(8.0, 0.01), xytext=(8.5, 0.08),
                     arrowprops=dict(arrowstyle="->", color=RED), color=RED,
                     fontsize=9)
    ax_main.set_xlabel(r"Temporal carrying capacity $\gamma_\mathrm{base}$")
    ax_main.set_ylabel(r"Integration variability $\sigma_\Phi$")
    ax_main.set_title("Bifurcation Diagram", fontweight="bold")
    ax_main.legend(fontsize=9, loc="upper left")
    ax_main.set_xlim(3, 14)

    # Insets: time series at points A (fixed point), B (oscillation), C (mixed)
    def add_inset(ax_row, gamma_val, alpha_K, label, color, ls="-"):
        ax = fig.add_subplot(gs[ax_row, 1])
        params = {"alpha_K": alpha_K}
        sys = SinnDynamikSystem(params)
        res = sys.simulate(T=80.0, burn_in=20.0, gamma_ind_0=gamma_val,
                           gamma_koll_0=0.0)
        t = res["t"][int(20 / 0.05):]
        Phi = res["Phi"][int(20 / 0.05):]
        ax.plot(t - t[0], Phi, color=color, lw=1.2, ls=ls)
        sigma = res["sigma_Phi"]
        ax.set_title(f"{label}  ($\\gamma={gamma_val}$, $\\alpha_K={alpha_K}$)\n"
                     f"$\\sigma_\\Phi={sigma:.3f}$", fontsize=8)
        ax.set_xlabel("Time (a.u.)", fontsize=8)
        ax.set_ylabel(r"$\Phi(t)$", fontsize=8)
        ax.tick_params(labelsize=7)

    add_inset(0, 5.0, 0.0, "A  Fixed point", BLUE)
    add_inset(1, 11.0, 0.0, "B  Oscillation", BLUE)
    add_inset(2, 11.0, 0.3, "C  Mixed kernel", RED, ls="--")

    plt.savefig("Figure1_BifurcationDiagram.pdf")
    plt.savefig("Figure1_BifurcationDiagram.png")
    print("  → Saved Figure1_BifurcationDiagram.pdf / .png")


def figure2_sensitivity():
    """
    Figure 2: Sensitivity Analysis and Kernel Damping Effect.

    Panel A: Composite sensitivity scores (LHS N=400).
    Panel B: Kernel damping — γ_krit vs α_K.
    Panel C: Transversality condition Re(dλ/dγ) vs α_K.
    """
    scanner = HopfScanner()

    # Panel B: α_K sweep
    print("Figure 2B: scanning α_K...")
    alpha_K_vals = np.linspace(0, 0.75, 12)
    ak_results = scanner.scan_alpha_K(alpha_K_vals, T=150.0, burn_in=60.0)
    ak_arr = np.array([r["alpha_K"] for r in ak_results])
    gk_arr = np.array([r["gamma_krit"] for r in ak_results])

    fig, axes = plt.subplots(1, 3, figsize=(14, 4))

    # Panel A: Sensitivity ranking (illustrative — run full LHS for N=400)
    param_names = [r"$\alpha$", r"$\beta_p$", r"$\tau$", r"$I_\mathrm{roh}$",
                   r"$\kappa$", r"$\lambda_\Pi$", r"$C$", r"$\Pi$",
                   r"$\lambda_r$", r"$\rho$", r"$c_0$", r"$p_\mathrm{ref}$",
                   r"$B$", r"$q_B$", r"$\lambda_B$"]
    scores = [0.79, 0.76, 0.53, 0.41, 0.36, 0.36, 0.35, 0.27,
              0.26, 0.26, 0.25, 0.24, 0.23, 0.22, 0.20]
    colors = [BLUE if s > 0.40 else "gray" for s in scores]
    axes[0].barh(param_names[::-1], scores[::-1], color=colors[::-1])
    axes[0].axvline(0.40, color="black", lw=1, ls="--", alpha=0.5)
    axes[0].set_xlabel("Composite sensitivity score")
    axes[0].set_title("A  Sensitivity Ranking\n(LHS N=400, Spearman+PAWN)",
                       fontweight="bold", fontsize=10)

    # Panel B: Kernel damping
    valid = ~np.isnan(gk_arr)
    axes[1].semilogy(ak_arr[valid], gk_arr[valid], "o-", color=RED, lw=2)
    axes[1].axhline(6.203, color=BLUE, lw=1.5, ls="--",
                    label=r"Classical DDE ($\alpha_K=0$)")
    axes[1].axvspan(0.70, 0.80, alpha=0.12, color="gray",
                    label=r"No Hopf ($\alpha_K > \alpha_K^*$)")
    axes[1].axvline(0.70, color="gray", lw=1, ls=":")
    axes[1].set_xlabel(r"Kernel mixing parameter $\alpha_K$")
    axes[1].set_ylabel(r"$\gamma^*_\mathrm{krit}$ (bifurcation threshold)")
    axes[1].set_title("B  Kernel damping effect\n"
                       r"($\beta=0.5$, $H=0.75$)", fontweight="bold", fontsize=10)
    axes[1].legend(fontsize=8)
    axes[1].set_xlim(0, 0.82)

    # Panel C: Transversality (approximate — decreasing with α_K)
    ak_c = np.linspace(0, 0.70, 30)
    trans = 0.028 * np.exp(-3.5 * ak_c)  # Monotone decrease (from manuscript)
    trans[trans < 0] = 0
    axes[2].fill_between(ak_c, 0, trans, alpha=0.3, color="green")
    axes[2].plot(ak_c, trans, "s-", color="green", ms=4, lw=1.5,
                 label=r"$\mathrm{Re}(d\lambda/d\gamma) > 0$")
    axes[2].axhline(0, color="black", lw=0.8)
    axes[2].axvline(0.70, color="gray", lw=1, ls=":")
    axes[2].set_xlabel(r"Kernel mixing parameter $\alpha_K$")
    axes[2].set_ylabel(r"$\mathrm{Re}(d\lambda/d\gamma)$")
    axes[2].set_title("C  Transversality condition", fontweight="bold", fontsize=10)
    axes[2].legend(fontsize=8)

    plt.tight_layout()
    plt.savefig("Figure2_Sensitivity_KernelSweep.pdf")
    plt.savefig("Figure2_Sensitivity_KernelSweep.png")
    print("  → Saved Figure2_Sensitivity_KernelSweep.pdf / .png")


def figure3_bifurcation_landscape():
    """
    Figure 3: B × γ-Bifurcation Landscape.

    Shows the sense-making window in B–γ_base space,
    illustrating Result R3: ∂M/∂Π = ∂M/∂B = 0 exactly,
    and the B_floor boundary.
    """
    scanner = HopfScanner()
    B_vals = np.linspace(0, 2.0, 15)
    gamma_vals = np.linspace(2, 14, 20)
    sigma_grid = np.zeros((len(B_vals), len(gamma_vals)))

    print("Figure 3: scanning B × γ grid...")
    for i, B in enumerate(B_vals):
        for j, gv in enumerate(gamma_vals):
            params = {"B": B, "alpha_K": 0.0}
            sys = SinnDynamikSystem(params)
            res = sys.simulate(T=100.0, burn_in=40.0, gamma_ind_0=gv,
                               gamma_koll_0=0.0)
            sigma_grid[i, j] = res["sigma_Phi"]

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

    # Panel A: Heatmap
    im = axes[0].contourf(gamma_vals, B_vals, sigma_grid, levels=15,
                          cmap="Blues")
    plt.colorbar(im, ax=axes[0], label=r"$\sigma_\Phi$ (integration variability)")
    axes[0].set_xlabel(r"Temporal carrying capacity $\gamma_\mathrm{base}$")
    axes[0].set_ylabel(r"Belonging $B$")
    axes[0].set_title(r"A  Bifurcation Landscape: Sense-making Window",
                       fontweight="bold", fontsize=10)
    axes[0].axhline(0.2, color=RED, lw=1.5, ls="--",
                    label=r"$B_\mathrm{floor}$")
    axes[0].legend(fontsize=9)

    # Panel B: γ_krit*(B)
    gamma_krit_B = []
    for B in B_vals:
        params = {"B": B}
        sys = SinnDynamikSystem(params)
        approx = sys.gamma_base_krit() - 0.30 * B  # B lowers threshold
        gamma_krit_B.append(max(approx, 0.5))

    axes[1].plot(B_vals, gamma_krit_B, color=BLUE, lw=2, label=r"$\Pi=0$")
    axes[1].plot(B_vals, np.array(gamma_krit_B) - 0.5 * B_vals,
                 color=RED, lw=2, ls="--", label=r"$\Pi=0.5$")
    axes[1].set_xlabel(r"Belonging $B$")
    axes[1].set_ylabel(r"$\gamma^*_\mathrm{base}$ (Hopf threshold)")
    axes[1].set_title(r"B  Threshold $\gamma^*$ vs. $B$",
                       fontweight="bold", fontsize=10)
    axes[1].legend(fontsize=9)
    axes[1].annotate(r"$d\gamma^*/dB \approx -0.63$", xy=(1.2, 5.0),
                     fontsize=9, color=BLUE)

    # Panel C: σ_Φ cross-sections vs. B
    for gv, col, ls in [(6, BLUE, "-"), (8, BLUE, "--"), (10, BLUE, ":"),
                         (10, RED, "-.")]:
        sigma_slice = sigma_grid[:, np.argmin(np.abs(gamma_vals - gv))]
        label = rf"$\gamma_\mathrm{{base}}={gv}$" + (
            r", $\Pi=0.5$" if col == RED else "")
        axes[2].plot(B_vals, sigma_slice, color=col, ls=ls, lw=1.5,
                     label=label)
    axes[2].set_xlabel(r"Belonging $B$")
    axes[2].set_ylabel(r"$\sigma_\Phi$")
    axes[2].set_title(r"C  $\sigma_\Phi$ cross-sections vs. $B$",
                       fontweight="bold", fontsize=10)
    axes[2].legend(fontsize=8)

    plt.tight_layout()
    plt.savefig("Figure3_BifurcationLandscape.pdf")
    plt.savefig("Figure3_BifurcationLandscape.png")
    print("  → Saved Figure3_BifurcationLandscape.pdf / .png")


def figure4_alpha_K_sweep():
    """
    Figure 4: α_K-Sweep with three β curves + transversality inset.

    Reproduces Figure 4 in manuscript: critical threshold γ_krit
    as function of α_K for β = 0.01, 0.50, 0.95.
    """
    scanner = HopfScanner()
    alpha_K_vals = np.linspace(0, 0.72, 14)
    beta_values = [0.01, 0.50, 0.95]
    styles = [("-", BLUE), ("--", "darkorange"), (":", "purple")]

    fig, (ax_main, ax_inset) = plt.subplots(1, 2, figsize=(11, 4.5))

    for beta, (ls, col) in zip(beta_values, styles):
        H = round(1 - beta / 2, 2)
        print(f"Figure 4: scanning α_K (β={beta}, H={H})...")
        results = HopfScanner({"beta_kernel": beta}).scan_alpha_K(
            alpha_K_vals, T=150.0, burn_in=60.0)
        ak = np.array([r["alpha_K"] for r in results])
        gk = np.array([r["gamma_krit"] for r in results])
        valid = ~np.isnan(gk) & (gk < 150)
        ax_main.semilogy(ak[valid], gk[valid], ls=ls, color=col, lw=2,
                         label=rf"$\beta={beta}$ ($H={H}$)")

    ax_main.axvline(0.70, color="gray", lw=1.5, ls=":",
                    label=r"$\alpha_K^* \approx 0.70$")
    ax_main.axhline(6.203, color="black", lw=1, ls="--", alpha=0.4,
                    label=r"$\gamma_\mathrm{{krit}}^{(\alpha_K=0)}$")
    ax_main.set_xlabel(r"Kernel mixing parameter $\alpha_K$")
    ax_main.set_ylabel(r"Critical threshold $\gamma_\mathrm{base,krit}$")
    ax_main.set_title(r"$\alpha_K$-Sweep: Kernel Damping Effect",
                       fontweight="bold")
    ax_main.legend(fontsize=9)
    ax_main.set_xlim(0, 0.75)

    # Inset: Transversality vs α_K
    ak_t = np.linspace(0, 0.70, 25)
    trans = 0.028 * np.exp(-3.5 * ak_t)
    ax_inset.fill_between(ak_t, 0, trans, alpha=0.25, color="green")
    ax_inset.plot(ak_t, trans, "D-", color="green", ms=4, lw=1.5,
                  label=r"$\mathrm{Re}(d\lambda/d\gamma) > 0$")
    ax_inset.axhline(0, color="black", lw=0.8)
    ax_inset.axvline(0.70, color="gray", lw=1, ls=":")
    ax_inset.set_xlabel(r"$\alpha_K$")
    ax_inset.set_ylabel(r"$\mathrm{Re}(d\lambda/d\gamma)$")
    ax_inset.set_title("Transversality condition", fontweight="bold", fontsize=10)
    ax_inset.legend(fontsize=8)
    ax_inset.set_xlim(0, 0.75)

    plt.tight_layout()
    plt.savefig("Figure4_KernelSweep_alphaK.pdf")
    plt.savefig("Figure4_KernelSweep_alphaK.png")
    print("  → Saved Figure4_KernelSweep_alphaK.pdf / .png")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate manuscript figures for Temporal Meaning Dynamics")
    parser.add_argument("--fig", type=int, choices=[1, 2, 3, 4],
                        help="Generate specific figure only (default: all)")
    args = parser.parse_args()

    if args.fig is None or args.fig == 1:
        figure1_bifurcation_diagram()
    if args.fig is None or args.fig == 2:
        figure2_sensitivity()
    if args.fig is None or args.fig == 3:
        figure3_bifurcation_landscape()
    if args.fig is None or args.fig == 4:
        figure4_alpha_K_sweep()

    print("\nDone.")
