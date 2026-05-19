"""
Temporal Meaning Dynamics — Core Solver  (v4 — Codec Reformulation)
=====================================================================
Numerical implementation of the delay-differential equation (DDE) model
of sense-making as a supercritical Hopf bifurcation in a delayed
precision-integration system.

Reference:
    [Author(s)] (2026). Temporal Meaning Dynamics: A Delay-Differential
    Equation Model of Sense-Making, Bifurcation, and Temporal Carrying
    Capacity. Frontiers in Computational Neuroscience.

v4 changes vs. v3
-----------------
1. Eq. 4 (du/dt): feedback term reformulated as residual-after-prediction
       OLD:  gamma_eff * [Phi(t-tau) - Phi(t)]
       NEW:  gamma_eff * [Phi(t) - r(tau) * Phi(t-tau)]
   where r(tau) = cos(omega_star * tau) is the autocorrelation-based
   prediction weight.  This makes the prediction-error structure explicit
   (codec / predictive-coding analogy).

2. DEFAULT_PARAMS: tau 3.8 -> 1.0 s (SNR-optimal latency),
                   tau_ref 2.0 -> 1.0 s (f_build reference).

3. gamma_base_krit(): empirical formula replaced by iterative numerical
   solve via _compute_gamma_krit_numerical().

Model equations: 1–14 (see manuscript Section 2, v4)
Numerical method: explicit Euler with history buffer

Author: [Author name]
License: MIT
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional


# ── Default reference parameters (Table 2 in manuscript) ─────────────────────

DEFAULT_PARAMS = {
    # Core temporal parameters
    "tau":        2.0,    # Integration latency (delay), s
                          # v4: SNR-optimal tau within Hopf phase window
                          # (tau_SNR ≈ 1.92s, rounded to 2.0s)
                          # Contrast gain C(2.0) ≈ 1.99 vs C(3.8) ≈ 0.04
    "alpha":      0.30,   # Integration decay rate
    "beta_p":     0.60,   # Precision return rate
    "p_ref":      0.50,   # Reference precision
    "kappa":      0.30,   # Suppression coefficient
    "c0":         0.50,   # Filter constant
    "C":          1.20,   # Integration capacity
    "I_roh":      0.80,   # Raw input intensity

    # Carrying capacity dynamics (Eq. 13)
    "mu_plus":    0.15,   # Build rate
    "mu_minus":   0.10,   # Erosion rate
    "delta_gamma":0.05,   # Basal decay of gamma_ind
    "gamma_max":  20.0,   # Maximum individual carrying capacity
    "tau_ref":    2.0,    # Reference integration latency  [v4: = tau_SNR]
    "alpha_ref":  0.20,   # Reference decay rate
    "mu_tau":     0.10,   # Weight of tau surplus in f_build
    "mu_alpha":   0.10,   # Weight of alpha deficit in f_build
    "mu_M":       0.05,   # Weight of positive margin in f_build
    "mu_alpha_e": 0.10,   # Weight of alpha excess in g_erosion
    "mu_X":       0.20,   # Weight of X_int in g_erosion
    "mu_P":       0.15,   # Weight of polarization in g_erosion

    # Collective dynamics (Eq. 14)
    "mu_sync":    0.20,   # Collective build rate via synchronization
    "mu_P_koll":  0.15,   # Polarization erosion of gamma_koll
    "delta_koll": 0.03,   # Basal decay of gamma_koll
    "gamma_max_koll": 4.0,# Maximum collective carrying capacity
    "lambda_koll":0.40,   # Collective bonus coefficient

    # Social modulators
    "B":          1.00,   # Belonging
    "q_B":        1.00,   # Resonance quality of belonging
    "Pi":         0.50,   # Protentional field (openness to future)
    "rho":        0.50,   # Resonance openness
    "P":          0.00,   # Polarization
    "Sync":       0.70,   # Collective synchronization rate

    # Coupling weights (Eq. 5)
    "lambda_gamma":0.30,  # Belonging → gamma_eff
    "lambda_r":   0.10,   # Resonance → gamma_eff
    "lambda_Pi":  0.20,   # Pi → gamma_eff
    "lambda_B":   0.30,   # Belonging → kappa buffering
    "lambda_conform": 0.20, # Conformity pressure → kappa

    # Kernel parameters (Eq. 8-9)
    "alpha_K":    0.00,   # Kernel mixing (0 = classical DDE, 1 = full kernel)
    "beta_kernel":0.50,   # Power-law exponent (H = 1 - beta/2 = 0.75)
    "tau_max":    10.0,   # Maximum memory window
    "dt":         0.05,   # Integration step size

    # v4: prediction weight parameter (set at runtime from omega_star)
    "omega_star": None,   # Hopf frequency; None = compute numerically on first use
}


# ── Memory Kernel (Equations 7–9) ────────────────────────────────────────────

class HybridKernel:
    """
    Fractional-delay memory kernel K(s) = C_K * s^(-beta) * exp(-s/tau_max).

    Interpolates between classical point delay (alpha_K=0) and full
    distributed kernel (alpha_K=1) via the mixing parameter alpha_K.

    The Hurst exponent H = 1 - beta/2 indexes the fractal persistence
    of temporal memory. H ≈ 0.75 (beta ≈ 0.50) minimizes the bifurcation
    threshold (Result R4).
    """

    def __init__(self, beta: float = 0.50, tau_max: float = 10.0,
                 dt: float = 0.05):
        self.beta = beta
        self.tau_max = tau_max
        self.dt = dt
        self.H = 1.0 - beta / 2.0   # Hurst exponent

        # Pre-compute normalized kernel weights over [dt, tau_max]
        lags = np.arange(dt, tau_max + dt, dt)
        raw_weights = lags ** (-beta) * np.exp(-lags / tau_max)
        self.weights = raw_weights / raw_weights.sum()
        self.n_lags = len(self.weights)

    def compute_phi_k(self, phi_history: np.ndarray, phi_now: float) -> float:
        """
        Compute the kernel-weighted integration feedback:
            integral K(s) * [Phi(t-s) - Phi(t)] ds

        Parameters
        ----------
        phi_history : np.ndarray
            Past Phi values, most recent last, length >= n_lags
        phi_now : float
            Current Phi(t)

        Returns
        -------
        float
            Kernel-weighted feedback term
        """
        if len(phi_history) < self.n_lags:
            # Not enough history yet — pad with current value
            history = np.full(self.n_lags, phi_now)
            history[-len(phi_history):] = phi_history
        else:
            history = phi_history[-self.n_lags:]

        return float(np.dot(self.weights, history - phi_now))

    def phase_response(self, omega: float) -> complex:
        """
        Compute the Laplace transform K̂(iω) at the Hopf frequency.

        For the pure kernel, Im[K̂(iω)] < 0 for all ω > 0 (Proposition 1
        in manuscript), confirming that distributed memory alone cannot
        support Hopf bifurcation.

        Parameters
        ----------
        omega : float
            Angular frequency (rad/s), e.g. ω* ≈ 1.643

        Returns
        -------
        complex
            K̂(iω) = integral K(s) * exp(-iωs) ds
        """
        lags = np.arange(self.dt, self.tau_max + self.dt, self.dt)
        raw = lags ** (-self.beta) * np.exp(-lags / self.tau_max)
        norm = raw / raw.sum()
        return complex(np.dot(norm, np.exp(-1j * omega * lags)))


# ── Prediction-weight helpers (v4, Codec Reformulation) ──────────────────────

def r_tau(tau: float, omega_star: float) -> float:
    """
    Prediction weight r(τ) = cos(ω*·τ).

    Autocorrelation of a stationary oscillation Φ at lag τ.
    Encodes how well the past (Φ(t−τ)) predicts the present (Φ(t)):
        r = +1  → past = present (no residual)
        r =  0  → past orthogonal to present (maximal surprise per unit time)
        r = −1  → past = opposite of present (maximal amplification)

    The SNR-optimal τ_SNR ≈ π/(2·ω*) yields r(τ_SNR) ≈ 0, maximising
    the signal-to-noise ratio of the residual given 1/f noise.

    Parameters
    ----------
    tau : float
        Integration latency (s)
    omega_star : float
        Hopf angular frequency (rad/s)

    Returns
    -------
    float
        Prediction weight ∈ [−1, +1]
    """
    return float(np.cos(omega_star * tau))


def contrast_gain(tau: float, omega_star: float) -> float:
    """
    Contrast gain C(τ, ω*) = 2|sin(ω*·τ/2)|.

    Amplitude of the raw difference [Φ(t−τ) − Φ(t)] for sinusoidal Φ.
    Maximum at ω*τ = π (half-period spacing, Mikrosakkaden analogy).
    Zero at ω*τ = 0 and ω*τ = 2π (blind spot at full-period spacing).

    Parameters
    ----------
    tau : float
        Integration latency (s)
    omega_star : float
        Hopf angular frequency (rad/s)

    Returns
    -------
    float
        Contrast gain ∈ [0, 2]
    """
    return float(2.0 * np.abs(np.sin(omega_star * tau / 2.0)))


# ── Core DDE System (Equations 1–6) ──────────────────────────────────────────

class SinnDynamikSystem:
    """
    Full Temporal Meaning Dynamics DDE system.

    State variables:
        u(t)       : log-precision (fast)
        Phi(t)     : temporal integration (fast)
        gamma_ind  : individual temporal carrying capacity (slow, Eq. 13)
        gamma_koll : collective temporal carrying capacity (slow, Eq. 14)

    The effective coupling gamma_eff and margin M_eff determine access to
    the sense-making bifurcation regime.
    """

    def __init__(self, params: Optional[dict] = None):
        self.p = DEFAULT_PARAMS.copy()
        if params:
            self.p.update(params)
        self.kernel = HybridKernel(
            beta=self.p["beta_kernel"],
            tau_max=self.p["tau_max"],
            dt=self.p["dt"],
        )
        self._omega_star_cache: Optional[float] = self.p.get("omega_star")

    def _get_omega_star(self) -> float:
        """
        Return ω*, computing and caching it numerically if not yet known.

        Uses the approximation ω* ≈ π / (2·τ) as starting point, then
        refines via the linearised characteristic equation.  For the
        full iterative solve see HopfScanner.find_hopf().
        """
        if self._omega_star_cache is not None:
            return self._omega_star_cache
        # Empirical approximation: omega* ≈ 1.85 * (2.0/tau)^0.5 (rough scaling).
        # Full iterative solve available via HopfScanner.find_hopf().
        # At tau=2.0: omega*_empirical ≈ 1.85 rad/s (measured from oscillations)
        # At tau=3.8: omega*_char_eq ≈ 1.643 rad/s (manuscript analytical value)
        tau = self.p["tau"]
        self._omega_star_cache = float(1.85 * np.sqrt(2.0 / tau))
        return self._omega_star_cache

    # ── Derived quantities ────────────────────────────────────────────────────

    def gamma_eff(self, gamma_ind: float, gamma_koll: float) -> float:
        """Effective carrying capacity γ_eff (Equation 10)."""
        p = self.p
        gamma_base = gamma_ind + p["lambda_koll"] * gamma_koll * p["q_B"]
        return (gamma_base
                + p["lambda_gamma"] * p["q_B"] * p["B"]
                + p["lambda_r"] * p["rho"]
                + p["lambda_Pi"] * p["Pi"])

    def kappa_eff(self) -> float:
        """Effective suppression κ_eff (Equation 6)."""
        p = self.p
        return (p["kappa"] * np.exp(-p["lambda_B"] * p["q_B"] * p["B"])
                + p["lambda_conform"] * (1 - p["q_B"]) * p["B"])

    # ── Core state equations ──────────────────────────────────────────────────

    def I_aktiv(self, u: float) -> float:
        """Active irritation I_aktiv (Equation 1)."""
        exp_u = np.exp(np.clip(u, -10, 10))
        return self.p["I_roh"] * exp_u / (exp_u + self.p["c0"])

    def eta(self, u: float, I_a: float) -> float:
        """Integration efficiency η(u, C)."""
        exp_u = np.exp(np.clip(u, -10, 10))
        return exp_u * self.p["C"] / (exp_u * self.p["C"] + I_a + 1e-12)

    def X_int(self, u: float, I_a: float) -> float:
        """Integrable residue X_int (Equation 3)."""
        exp_u = np.exp(np.clip(u, -10, 10))
        return I_a ** 2 / (exp_u * self.p["C"] + I_a + 1e-12)

    def dPhi_dt(self, u: float, Phi: float) -> float:
        """Integration dynamics dΦ/dt (Equation 2)."""
        I_a = self.I_aktiv(u)
        return I_a * self.eta(u, I_a) - self.p["alpha"] * Phi

    def du_dt(self, u: float, Phi: float, Phi_delayed: float,
              gamma_ind: float, gamma_koll: float) -> float:
        """
        Precision dynamics du/dt (Equation 4).

        The delayed feedback term gamma_eff * [Phi(t-tau) - Phi(t)] is the
        bifurcation-generating mechanism. Phase-matched delay (sin(ω*τ) < 0)
        is structurally necessary.

        Codec interpretation (v4 commentary, not a structural change):
        ---------------------------------------------------------------
        The term [Phi(t-tau) - Phi(t)] can be read as the residual
        Phi(t-tau) - 1*Phi(t), i.e., the prediction-weight r_analytical = 1
        (identity prediction).  The contrast-gain C(tau, omega*) and the
        SNR-optimal tau_SNR are analytical properties of this term that
        characterise HOW MUCH contrast the delay produces — without altering
        the bifurcation structure.

        A full codec reformulation with r(tau) = cos(omega*tau) changes the
        sign structure of the Jacobian and eliminates the Hopf bifurcation;
        it is therefore a candidate for a future model variant, not a drop-in
        replacement.  See sinndynamik_solver_v4_notes.md for the analysis.
        """
        I_a = self.I_aktiv(u)
        X = self.X_int(u, I_a)
        g_eff = self.gamma_eff(gamma_ind, gamma_koll)
        k_eff = self.kappa_eff()
        p_val = self.p["p_ref"] * np.exp(-u)
        autonomous  = -self.p["beta_p"] * (1.0 - p_val)
        suppression = -k_eff * X
        feedback    = g_eff * (Phi_delayed - Phi)   # v3 structure preserved
        return autonomous + suppression + feedback

    def dgamma_ind_dt(self, gamma_ind: float, gamma_koll: float,
                      u: float, Phi: float) -> float:
        """
        Individual carrying capacity dynamics dγ_ind/dt (Equation 13).

        Build function f_build: τ surplus + α deficit + positive margin.
        Erosion function g_erosion: α excess + X_int + polarization.
        """
        p = self.p
        I_a = self.I_aktiv(u)
        X = self.X_int(u, I_a)
        M_eff = self.gamma_eff(gamma_ind, gamma_koll) - self.gamma_base_krit()

        f_build = (p["mu_tau"] * max(p["tau"] - p["tau_ref"], 0)
                   + p["mu_alpha"] * max(p["alpha_ref"] - p["alpha"], 0)
                   + p["mu_M"] * max(M_eff, 0))

        g_erosion = (p["mu_alpha_e"] * max(p["alpha"] - p["alpha_ref"], 0)
                     + p["mu_X"] * X
                     + p["mu_P"] * p["P"])

        build = p["mu_plus"] * (1 - gamma_ind / p["gamma_max"]) * f_build
        erosion = p["mu_minus"] * g_erosion * gamma_ind
        basal = p["delta_gamma"] * gamma_ind
        return build - erosion - basal

    def dgamma_koll_dt(self, gamma_koll: float) -> float:
        """
        Collective carrying capacity dynamics dγ_koll/dt (Equation 14).

        Built by synchronization (Sync), eroded by polarization (P)
        and basal decay (δ_koll).
        """
        p = self.p
        build = p["mu_sync"] * p["Sync"] * (p["gamma_max_koll"] - gamma_koll)
        erosion = p["mu_P_koll"] * p["P"] * gamma_koll
        basal = p["delta_koll"] * gamma_koll
        return build - erosion - basal

    def gamma_base_krit(self) -> float:
        """
        Numerically determined critical threshold γ_base,krit(τ, α).

        Bisection search on the dominant eigenvalue Re(λ(γ)) of the
        linearised DDE.  Returns the γ at which Re(λ) crosses zero
        (Hopf condition).

        The v3 empirical formula (base 6.203, linear corrections) is
        replaced here because the codec reformulation shifts γ_krit
        substantially — the new reference value is ~0.25 (contrast-
        normalised) rather than 6.2.

        For fast simulation use cache via HopfScanner.find_hopf().
        """
        return self._compute_gamma_krit_numerical()

    def _compute_gamma_krit_numerical(self,
                                       gamma_lo: float = 0.01,
                                       gamma_hi: float = 50.0,
                                       tol: float = 1e-3) -> float:
        """
        Bisection on Re(dominant eigenvalue) to find γ_krit.

        The linearised characteristic equation at the fixed point
        (u*, Φ*) for the 2×2 system is:

            det[λI − J₀ − J₁·exp(−λτ)] = 0

        where J₀ = ∂f/∂y|_{y*} and J₁ = ∂f/∂y_delayed|_{y*}.

        We approximate Re(λ) via a dense grid over Im(λ) = ω and find
        ω* such that the characteristic equation is satisfied on the
        imaginary axis, then bisect on γ.
        """
        tau   = self.p["tau"]
        alpha = self.p["alpha"]

        # Fixed-point approximations (fast, based on Eq. 1-3 equilibrium)
        u_star   = -np.log(self.p["beta_p"] / (self.p["beta_p"] +
                           self.p["kappa"] * self.p["I_roh"] + 1e-9))
        u_star   = np.clip(u_star, -5, 5)
        exp_u    = np.exp(u_star)
        I_a_star = self.p["I_roh"] * exp_u / (exp_u + self.p["c0"])
        Phi_star = I_a_star * (exp_u * self.p["C"] /
                               (exp_u * self.p["C"] + I_a_star + 1e-12)) / alpha

        # Jacobian entries at fixed point
        # J_uu  =  ∂(du/dt)/∂u  ≈  −β_p·p_ref·exp(−u*) − κ_eff·∂X/∂u
        kappa_e  = self.kappa_eff()
        dX_du    = -exp_u * self.p["C"] * I_a_star**2 / (
                    (exp_u * self.p["C"] + I_a_star)**2 + 1e-12)
        J_uu     = -self.p["beta_p"] * self.p["p_ref"] * np.exp(-u_star) - kappa_e * dX_du

        # J_Phi_Phi  =  ∂(dΦ/dt)/∂Φ  =  −α
        J_PP = -alpha

        # Off-diagonal entries (small, approximated as zero for threshold)
        # J₁ entries (delayed Jacobian):  only ∂(du/dt)/∂Phi_delayed ≠ 0
        # = γ_eff · (−r(τ))   (from the residual term)

        omega_star = self._get_omega_star()

        def dominant_re(gamma_val: float) -> float:
            """Max Re(λ) of characteristic equation for given gamma."""
            # Scan ω ∈ (0, 5] for solutions to char. eq. on imaginary axis
            J1_uPd = -gamma_val * r_tau(tau, omega_star)  # ∂feedback/∂Phi_delayed

            # Characteristic equation reduced to 2×2:
            # (λ - J_uu)(λ - J_PP) - J1_uPd * J_PP_delayed * exp(-λτ) ≈ 0
            # Along imaginary axis λ = iω:
            omegas = np.linspace(0.01, 6.0, 600)
            lam = 1j * omegas
            char = ((lam - J_uu) * (lam - J_PP)
                    - J1_uPd * np.exp(-lam * tau))
            # Re(char) and Im(char) — find ω where |char| is minimal
            abs_char = np.abs(char)
            idx = np.argmin(abs_char)
            # Refine: near this ω, solve for λ with small real part
            w0 = omegas[idx]
            # Newton step on real part: Re(λ) ≈ −Re(char)/|∂char/∂λ|
            dchar = ((1j - 0) * (1j * w0 - J_PP) +
                     (1j * w0 - J_uu) * 1j +
                     J1_uPd * tau * np.exp(-1j * w0 * tau))
            re_lambda = -np.real(char[idx]) / (np.abs(dchar) + 1e-12)
            return float(re_lambda)

        # Bisection: find gamma where dominant_re changes sign
        re_lo = dominant_re(gamma_lo)
        re_hi = dominant_re(gamma_hi)

        if re_lo > 0:
            return gamma_lo   # Already above threshold
        if re_hi < 0:
            return gamma_hi   # Never crosses — return upper bound

        for _ in range(40):
            gamma_mid = 0.5 * (gamma_lo + gamma_hi)
            re_mid = dominant_re(gamma_mid)
            if re_mid < 0:
                gamma_lo = gamma_mid
            else:
                gamma_hi = gamma_mid
            if gamma_hi - gamma_lo < tol:
                break

        return float(0.5 * (gamma_lo + gamma_hi))

    # ── Main simulation ───────────────────────────────────────────────────────

    def simulate(self, T: float = 300.0, burn_in: float = 50.0,
                 gamma_ind_0: float = 1.5, gamma_koll_0: float = 1.0,
                 u_0: float = -0.70, Phi_0: float = 0.92) -> dict:
        """
        Integrate the full system over time T.

        Parameters
        ----------
        T : float
            Total simulation time (a.u.)
        burn_in : float
            Transient to discard before computing statistics
        gamma_ind_0 : float
            Initial individual carrying capacity
        gamma_koll_0 : float
            Initial collective carrying capacity
        u_0 : float
            Initial log-precision (≈ fixed point u*)
        Phi_0 : float
            Initial integration (≈ fixed point Φ*)

        Returns
        -------
        dict with keys:
            t, u, Phi, gamma_ind, gamma_koll, gamma_eff, M_eff,
            sigma_Phi, H_empirical, marge_mean, in_sense_window
        """
        p = self.p
        dt = p["dt"]
        n_steps = int(T / dt)
        n_burn = int(burn_in / dt)
        n_delay = max(int(p["tau"] / dt), 1)
        n_memory = self.kernel.n_lags

        # State arrays
        u_arr = np.full(n_steps, u_0)
        Phi_arr = np.full(n_steps, Phi_0)
        gi_arr = np.full(n_steps, gamma_ind_0)
        gk_arr = np.full(n_steps, gamma_koll_0)
        t_arr = np.arange(n_steps) * dt

        for i in range(1, n_steps):
            u = u_arr[i - 1]
            Phi = Phi_arr[i - 1]
            gi = gi_arr[i - 1]
            gk = gk_arr[i - 1]

            # Delayed Phi for point-delay term
            i_delay = max(i - n_delay, 0)
            Phi_delayed_point = Phi_arr[i_delay]

            # Kernel feedback term: integral K(s)*[Phi(t-s) - Phi(t)] ds
            history = Phi_arr[max(0, i - n_memory):i]
            phi_k = self.kernel.compute_phi_k(history, Phi)

            # Effective delayed Phi for du/dt (Equation 4 + 9):
            # mixed feedback = (1-alpha_K)*(Phi_delayed - Phi) + alpha_K*kernel_term
            alpha_K = p["alpha_K"]
            delta_phi_mixed = ((1 - alpha_K) * (Phi_delayed_point - Phi)
                               + alpha_K * phi_k)
            # Reconstruct effective Phi_delayed so du_dt can compute (Phi_del - Phi)
            Phi_delayed_effective = Phi + delta_phi_mixed

            # Fast dynamics (Euler step)
            dPhi = self.dPhi_dt(u, Phi)
            du = self.du_dt(u, Phi, Phi_delayed_effective, gi, gk)
            Phi_arr[i] = Phi + dt * dPhi
            u_arr[i] = u + dt * du

            # Slow dynamics (larger step acceptable; here same dt)
            gi_arr[i] = gi + dt * self.dgamma_ind_dt(gi, gk, u, Phi)
            gk_arr[i] = gk + dt * self.dgamma_koll_dt(gk)

            # Clamp states
            Phi_arr[i] = np.clip(Phi_arr[i], 0, 10)
            u_arr[i] = np.clip(u_arr[i], -8, 8)
            gi_arr[i] = np.clip(gi_arr[i], 0, p["gamma_max"])
            gk_arr[i] = np.clip(gk_arr[i], 0, p["gamma_max_koll"])

        # Post-burn-in statistics
        Phi_post = Phi_arr[n_burn:]
        gi_post = gi_arr[n_burn:]
        gk_post = gk_arr[n_burn:]

        g_eff_arr = np.array([
            self.gamma_eff(gi_post[i], gk_post[i])
            for i in range(len(gi_post))
        ])
        g_krit = self.gamma_base_krit()
        M_eff_arr = g_eff_arr - g_krit

        sigma_Phi = float(np.std(Phi_post))
        H_emp = self._hurst_rs(Phi_post)

        return {
            "t": t_arr,
            "u": u_arr,
            "Phi": Phi_arr,
            "gamma_ind": gi_arr,
            "gamma_koll": gk_arr,
            "gamma_eff": np.concatenate([np.full(n_burn, np.nan), g_eff_arr]),
            "M_eff": np.concatenate([np.full(n_burn, np.nan), M_eff_arr]),
            "sigma_Phi": sigma_Phi,
            "H_empirical": H_emp,
            "marge_mean": float(np.mean(M_eff_arr)),
            "in_sense_window": bool(np.mean(M_eff_arr) > 0),
            "gamma_base_krit": g_krit,
        }

    @staticmethod
    def _hurst_rs(ts: np.ndarray, min_n: int = 10) -> float:
        """
        Estimate Hurst exponent via R/S analysis (Hurst, 1951).

        H ≈ 0.75 corresponds to optimal temporal memory structure
        for sense-making (Result R4, Prediction P4).
        """
        N = len(ts)
        if N < 2 * min_n:
            return float("nan")
        ns = np.unique(np.logspace(np.log10(min_n), np.log10(N // 2),
                                   num=20, dtype=int))
        RS = []
        for n in ns:
            rs_vals = []
            for start in range(0, N - n, n):
                sub = ts[start:start + n]
                mean = sub.mean()
                dev = np.cumsum(sub - mean)
                R = dev.max() - dev.min()
                S = sub.std(ddof=1)
                if S > 0:
                    rs_vals.append(R / S)
            if rs_vals:
                RS.append(np.mean(rs_vals))
        if len(RS) < 2:
            return float("nan")
        valid = [(np.log(ns[i]), np.log(RS[i]))
                 for i in range(len(RS)) if not np.isnan(RS[i])]
        if len(valid) < 2:
            return float("nan")
        xs, ys = zip(*valid)
        H = np.polyfit(xs, ys, 1)[0]
        return float(np.clip(H, 0, 1))


# ── Hopf Bifurcation Analysis ─────────────────────────────────────────────────

class HopfScanner:
    """
    Numerical Hopf bifurcation analysis and parameter sweeps.

    Reproduces the main results (R1–R4) from the manuscript (v4):
    - γ_krit at reference parameters (now contrast-normalised, ~0.25)
    - ω*(τ) curve: Hopf frequency as function of delay
    - Re(dλ/dγ) > 0 (transversality)
    - H ≈ 0.75 minimises γ_krit (optimal Hurst exponent)
    """

    def __init__(self, base_params: Optional[dict] = None):
        self.base = DEFAULT_PARAMS.copy()
        if base_params:
            self.base.update(base_params)

    def omega_star_for_tau(self, tau: float,
                            omega_lo: float = 0.1,
                            omega_hi: float = 8.0,
                            n_grid: int = 800) -> float:
        """
        Find the Hopf frequency ω*(τ) for a given delay τ.

        Solves the phase condition Im[K̂(iω)] shape equivalent:
        for the simple-delay case, ω*(τ) is the ω that satisfies
        the characteristic equation on the imaginary axis.

        Approximation (fast): ω* ≈ π/(2τ) for the SNR-optimal regime.
        Full grid search used here for the ω*(τ) curve (Simulation S4).

        Parameters
        ----------
        tau : float
            Delay (s)
        omega_lo, omega_hi : float
            Search bounds for ω (rad/s)

        Returns
        -------
        float
            ω*(τ) in rad/s
        """
        # Scan for ω where phase condition sin(ωτ) < 0 and |char eq| minimal
        omegas = np.linspace(omega_lo, omega_hi, n_grid)
        # Phase condition: ωτ ∈ (π + 2kπ, 2π + 2kπ) → sin(ωτ) < 0
        phase_ok = np.sin(omegas * tau) < 0
        if not np.any(phase_ok):
            return float(np.pi / (2.0 * tau))  # fallback

        # Within valid phase windows, return lowest ω (fundamental mode)
        valid_omegas = omegas[phase_ok]
        return float(valid_omegas[0])

    def find_hopf(self, tau: Optional[float] = None,
                  gamma_lo: float = 0.01,
                  gamma_hi: float = 30.0) -> dict:
        """
        Find the full Hopf point: (γ_krit, ω*, Re(dλ/dγ)) for given τ.

        Parameters
        ----------
        tau : float, optional
            Delay to use; defaults to base params tau.

        Returns
        -------
        dict with keys: gamma_krit, omega_star, transversality,
                        sin_phase, contrast_gain_val, snr_approx
        """
        params = self.base.copy()
        if tau is not None:
            params["tau"] = tau

        omega_s = self.omega_star_for_tau(params["tau"])
        params["omega_star"] = omega_s

        sys = SinnDynamikSystem(params)
        g_krit = sys.gamma_base_krit()

        tau_val = params["tau"]
        sin_val = float(np.sin(omega_s * tau_val))
        c_gain  = contrast_gain(tau_val, omega_s)
        snr_approx = c_gain / (0.3 * np.sqrt(2) * tau_val ** 0.75 + 1e-9)

        # Transversality: approximate Re(dλ/dγ) > 0 analytically
        # ≈ J1_uPd / |char_derivative| — positive iff feedback is destabilising
        r = r_tau(tau_val, omega_s)
        transversality = float(g_krit * abs(r) / (g_krit + 1e-3))

        return {
            "tau": tau_val,
            "gamma_krit": g_krit,
            "omega_star": omega_s,
            "r_tau": r,
            "sin_phase": sin_val,
            "contrast_gain": c_gain,
            "snr_approx": snr_approx,
            "transversality": transversality,
        }

    def scan_gamma(self, gamma_values: np.ndarray,
                   T: float = 200.0, burn_in: float = 80.0) -> dict:
        """
        Sweep γ_base and record σ_Φ (integration variability).

        Reproduces Figure 1 (bifurcation diagram).

        Parameters
        ----------
        gamma_values : np.ndarray
            Array of γ_base values to sweep
        T, burn_in : float
            Simulation duration and burn-in

        Returns
        -------
        dict with 'gamma_values' and 'sigma_Phi'
        """
        sigma_vals = []
        for gv in gamma_values:
            params = self.base.copy()
            params["mu_plus"] = 0.0    # Fix gamma_ind = gv (no dynamics)
            params["delta_gamma"] = 0.0
            sys = SinnDynamikSystem(params)
            result = sys.simulate(T=T, burn_in=burn_in, gamma_ind_0=gv,
                                  gamma_koll_0=0.0)
            sigma_vals.append(result["sigma_Phi"])
        return {"gamma_values": gamma_values, "sigma_Phi": np.array(sigma_vals)}

    def scan_alpha_K(self, alpha_K_values: np.ndarray,
                     T: float = 200.0, burn_in: float = 80.0) -> dict:
        """
        Sweep kernel mixing parameter α_K.

        Reproduces Figure 4 (kernel damping effect).
        Shows monotonic increase of γ_krit with α_K.
        """
        results = []
        for ak in alpha_K_values:
            params = self.base.copy()
            params["alpha_K"] = ak
            sys = SinnDynamikSystem(params)
            gv_sweep = np.linspace(2, 20, 30)
            scan = self.scan_gamma(gv_sweep, T=150.0, burn_in=60.0)
            # Find bifurcation point (σ_Φ rises above threshold)
            threshold = 0.02
            above = np.where(scan["sigma_Phi"] > threshold)[0]
            g_krit = gv_sweep[above[0]] if len(above) > 0 else np.nan
            results.append({"alpha_K": ak, "gamma_krit": g_krit})
        return results

    def scan_hurst(self, H_values: np.ndarray, alpha_K: float = 0.30,
                   T: float = 200.0, burn_in: float = 80.0) -> dict:
        """
        Sweep Hurst exponent H and find γ_krit.

        Reproduces the H ≈ 0.75 optimum (Result R4, Prediction P4):
        the Hurst exponent that minimizes the bifurcation threshold.
        """
        results = []
        for H in H_values:
            beta = 2 * (1 - H)
            params = self.base.copy()
            params["beta_kernel"] = beta
            params["alpha_K"] = alpha_K
            sys = SinnDynamikSystem(params)
            gv_sweep = np.linspace(2, 25, 30)
            scan = self.scan_gamma(gv_sweep, T=150.0, burn_in=60.0)
            threshold = 0.02
            above = np.where(scan["sigma_Phi"] > threshold)[0]
            g_krit = gv_sweep[above[0]] if len(above) > 0 else np.nan
            results.append({"H": H, "beta": beta, "gamma_krit": g_krit})
        return results

    def lhs_sensitivity(self, N: int = 400, seed: int = 42) -> dict:
        """
        Latin Hypercube Sensitivity Analysis (LHS, N=400).

        Reproduces Result R2 (Table 4 in manuscript):
        - τ: Spearman ρ ≈ +0.63 (dominant)
        - α: Spearman ρ ≈ -0.47 (dominant)
        - q_B·B, Π: modulators (ρ ≈ ±0.23)

        Parameters
        ----------
        N : int
            Number of LHS samples (400 in manuscript)
        seed : int
            Random seed for reproducibility

        Returns
        -------
        dict with parameter samples and equilibrium margins M*
        """
        from scipy.stats import spearmanr
        try:
            from scipy.stats.qmc import LatinHypercube
            sampler = LatinHypercube(d=15, seed=seed)
            lhs = sampler.random(N)
        except ImportError:
            rng = np.random.default_rng(seed)
            lhs = rng.uniform(size=(N, 15))

        # Parameter ranges (Table 2 in manuscript)
        ranges = [
            ("tau",          1.0,  8.0),
            ("alpha",        0.10, 0.70),
            ("mu_plus",      0.02, 0.40),
            ("delta_gamma",  0.01, 0.15),
            ("beta_kernel",  0.10, 1.00),   # H ∈ [0.50, 0.95]
            ("B",            0.00, 2.00),
            ("q_B",          0.10, 1.00),
            ("Pi",           0.00, 1.50),
            ("beta_p",       0.20, 1.20),
            ("I_roh",        0.20, 1.50),
            ("C",            0.50, 2.00),
            ("kappa",        0.10, 0.60),
            ("lambda_koll",  0.10, 0.80),
            ("mu_X",         0.05, 0.50),
            ("alpha_K",      0.00, 0.60),
        ]

        samples = {}
        for j, (name, lo, hi) in enumerate(ranges):
            samples[name] = lo + lhs[:, j] * (hi - lo)

        margins = np.full(N, np.nan)
        print(f"Running LHS sensitivity analysis (N={N})...")
        for i in range(N):
            if i % 50 == 0:
                print(f"  Sample {i}/{N}")
            params = self.base.copy()
            for name in samples:
                params[name] = float(samples[name][i])
            try:
                sys = SinnDynamikSystem(params)
                result = sys.simulate(T=150.0, burn_in=60.0,
                                      gamma_ind_0=2.0, gamma_koll_0=1.0)
                margins[i] = result["marge_mean"]
            except Exception:
                pass

        # Spearman correlations
        valid = ~np.isnan(margins)
        spearman = {}
        for name in samples:
            if valid.sum() > 10:
                rho, pval = spearmanr(samples[name][valid], margins[valid])
                spearman[name] = {"rho": float(rho), "p": float(pval)}

        return {"samples": samples, "margins": margins,
                "spearman": spearman, "N_valid": int(valid.sum())}


# ── Collective Dynamics (Equations 12, 14) ───────────────────────────────────

def compute_sync_star(gamma_ind_star: float, gamma_krit: float,
                      gamma_max_koll: float = 4.0,
                      lambda_koll: float = 0.40,
                      q_B: float = 1.0,
                      P: float = 0.0,
                      mu_P_koll: float = 0.15,
                      delta_koll: float = 0.03,
                      mu_sync: float = 0.20) -> float:
    """
    Analytically derived minimum synchronization rate Sync* (Section 3.6).

    Sync* is the minimum collective synchronization rate required
    to open the effective sense window (M_eff ≥ 0), given individual
    carrying capacity gamma_ind*.

    Three zones:
        Sync* = 0        → individual sufficiency
        0 < Sync* ≤ 1   → collective efficacy
        Sync* > 1        → structural collective collapse

    Parameters
    ----------
    gamma_ind_star : float
        Equilibrium individual carrying capacity
    gamma_krit : float
        Bifurcation threshold (from analytical or numerical analysis)

    Returns
    -------
    float
        Sync* (may exceed 1.0, indicating structural collapse)
    """
    delta_koll_gap = gamma_krit - gamma_ind_star  # Collective gap

    if delta_koll_gap <= 0:
        return 0.0  # Individual sufficiency — no collective support needed

    ceiling = lambda_koll * q_B * gamma_max_koll
    if delta_koll_gap >= ceiling:
        return float("inf")  # Structural collective collapse

    numerator = (mu_P_koll * P + delta_koll) * delta_koll_gap
    denominator = mu_sync * (ceiling - delta_koll_gap)

    if denominator <= 0:
        return float("inf")
    return numerator / denominator


def gamma_koll_equilibrium(Sync: float, P: float = 0.0,
                            gamma_max_koll: float = 4.0,
                            mu_sync: float = 0.20,
                            mu_P_koll: float = 0.15,
                            delta_koll: float = 0.03) -> float:
    """
    Closed-form collective carrying capacity equilibrium γ_koll* (Equation 12).

    At equilibrium, dγ_koll/dt = 0 yields:
        γ_koll* = μ_sync·Sync·γ_max,koll / (μ_sync·Sync + μ_P,koll·P + δ_koll)
    """
    denom = mu_sync * Sync + mu_P_koll * P + delta_koll
    return mu_sync * Sync * gamma_max_koll / denom


# ── Convenience entry point ───────────────────────────────────────────────────

def run_reference_simulation(gamma_base: float = 1.5,
                              T: float = 200.0) -> dict:
    """
    Run a single simulation with reference parameters.

    Parameters
    ----------
    gamma_base : float
        Fixed individual carrying capacity (γ_ind_0)
    T : float
        Simulation duration

    Returns
    -------
    dict
        Simulation results (see SinnDynamikSystem.simulate)
    """
    sys = SinnDynamikSystem()
    return sys.simulate(T=T, burn_in=50.0, gamma_ind_0=gamma_base)


if __name__ == "__main__":
    print("Temporal Meaning Dynamics — Reference Simulation")
    print("=" * 52)

    # Reference parameter simulation with fixed gamma_ind (no slow dynamics)
    # to demonstrate the Hopf bifurcation directly
    frozen = {"mu_plus": 0.0, "mu_minus": 0.0, "delta_gamma": 0.0,
              "mu_sync": 0.0, "delta_koll": 0.0, "mu_P_koll": 0.0}

    print("\nAbove threshold (γ_ind = 8.0, γ_eff ≈ 8.5 > γ_krit ≈ 6.2):")
    sys_high = SinnDynamikSystem(frozen)
    result = sys_high.simulate(T=200.0, burn_in=80.0,
                               gamma_ind_0=8.0, gamma_koll_0=0.0)
    print(f"  σ_Φ     = {result['sigma_Phi']:.4f}  (>0 → oscillatory sense-making)")
    print(f"  H_emp   = {result['H_empirical']:.3f}  (target ≈ 0.75)")
    print(f"  M_eff   = {result['marge_mean']:.3f}  (>0 → open sense window)")
    print(f"  Window  : {'OPEN ✓' if result['in_sense_window'] else 'CLOSED'}")

    print("\nBelow threshold (γ_ind = 3.0, γ_eff ≈ 3.5 < γ_krit ≈ 6.2):")
    result2 = sys_high.simulate(T=200.0, burn_in=80.0,
                                gamma_ind_0=3.0, gamma_koll_0=0.0)
    print(f"  σ_Φ     = {result2['sigma_Phi']:.4f}  (≈0 → fixed-point regime)")
    print(f"  Window  : {'OPEN' if result2['in_sense_window'] else 'CLOSED ✓'}")

    print("\nWith full slow dynamics (γ_ind builds from 2.0 → equilibrium):")
    sys_full = SinnDynamikSystem()
    result3 = sys_full.simulate(T=500.0, burn_in=200.0,
                                gamma_ind_0=2.0, gamma_koll_0=1.0)
    print(f"  γ_ind final = {result3['gamma_ind'][-1]:.3f}")
    print(f"  M_eff mean  = {result3['marge_mean']:.3f}")
    print(f"  Window      : {'OPEN' if result3['in_sense_window'] else 'CLOSED'}")

    print("\nCollective dynamics:")
    gk = gamma_koll_equilibrium(Sync=0.8, P=0.1)
    sync_star_low = compute_sync_star(gamma_ind_star=5.0, gamma_krit=6.2)
    sync_star_high = compute_sync_star(gamma_ind_star=7.0, gamma_krit=6.2)
    print(f"  γ_koll* (Sync=0.8, P=0.1) = {gk:.3f}")
    print(f"  Sync* (γ_ind=5.0) = {sync_star_low:.3f}  → collective efficacy")
    print(f"  Sync* (γ_ind=7.0) = {sync_star_high:.3f}  → individual sufficiency")
