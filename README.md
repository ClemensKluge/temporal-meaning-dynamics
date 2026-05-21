# Temporal Meaning Dynamics

**Code repository for:**  
*Temporal Meaning Dynamics: A Delay-Differential Equation Model of Sense-Making, Bifurcation, and Temporal Carrying Capacity*  
Frontiers in Computational Neuroscience — Hypothesis and Theory (resubmitted 2026)

---

## Overview

This repository contains the numerical implementation of the Temporal Meaning Dynamics (TMD) model — a delay-differential equation (DDE) system in which sense-making is formally represented as a **subcritical** Hopf bifurcation governed by temporal carrying capacity and phase-matched feedback.

The model formalises the theoretical proposal that meaning-constitution requires a minimum capacity to sustain temporal openness (γ_eff ≥ γ_krit), and generates six empirically testable predictions (P1–P6) linking HRV spectral structure, EEG microstates, and narrative coherence to the bifurcation threshold.

> **v4.4 note:** Earlier versions of the README and manuscript described a *supercritical* Hopf bifurcation (amplitude ~ √(γ − γ_krit)). The first Lyapunov coefficient computed via the Kuznetsov–Faria–Magalhães formula gives **l₁ = +7.53 > 0**, establishing the bifurcation as **subcritical** (jump onset, bistability, hysteresis). All values below reflect the v4.4 revision.

---

## Quick Start

```bash
# Install dependencies
pip install -r requirements.txt

# Reproduce all manuscript values (10/10 assertions)
python reproduce_all.py

# Regenerate all figures
python generate_figures_v44.py

# Run a specific simulation
python sinndynamik_solver_v4.py
```

---

## Files

| File | Description |
|---|---|
| `sinndynamik_solver_v4.py` | Core DDE solver — all model equations (Eq. 1–14), Hopf scanner, sensitivity analysis |
| `reproduce_all.py` | **Reproducibility supplement** — recomputes all 10 key manuscript values with assertion checks (10/10 PASS) |
| `params_reference.json` | Machine-readable reference parameters and expected values for automated testing |
| `generate_figures_v44.py` | Reproduces Figures 1, 1c, 2, 3, 4, and S1 from the manuscript |
| `requirements.txt` | Python dependencies (`numpy`, `scipy`, `matplotlib`) |

---

## Reproducibility

Running `python reproduce_all.py` reproduces all quantitative claims in the manuscript. Expected output:

```
[PASS] Φ* (fixed point):            computed=0.7582, expected=0.758 ±0.003
[PASS] u* (fixed point):            computed=-0.7689, expected=-0.769 ±0.003
[PASS] J_uu (regulation):           computed=-0.658, expected=-0.658 ±0.025
[PASS] γ_krit (Hopf threshold):     computed=7.186, expected=7.186 ±0.15
[PASS] ω* (Hopf frequency):         computed=1.215, expected=1.215 ±0.015
[PASS] sin(ω*τ) (phase cond.):      computed=-0.996, expected=-0.996 ±0.005
[PASS] Re(dλ/dγ) (transversality):  computed=0.022, expected=0.021 ±0.015
[PASS] l₁ (Lyapunov, subcritical):  computed=7.53,  expected=7.09  ±0.60  → SUBCRITICAL ✓
[PASS] ∂M/∂Π (exact = λ_Π):        computed=0.200, expected=0.200 ±0.001
[PASS] ∂M/∂B (B=1, q_B=1):         computed=0.51,  expected=0.48  ±0.10

Result: 10 PASS / 0 FAIL
```

---

## Model Equations

The core system is a 2-variable DDE with delayed feedback:

```
dΦ/dt  =  I_aktiv(u) · η(u, C)  −  α · Φ                        (Eq. 2)

du/dt  =  −β_p · (1 − p_ref · exp(−u))
          −  κ_eff(B, q_B) · X_int(u)
          +  γ_eff · [Φ(t−τ) − Φ(t)]                            (Eq. 4)
```

where Φ is the integration state, u = log(p) is log-precision, τ is the integration latency, γ_eff is the effective temporal carrying capacity, and κ_eff encodes the belonging-modulated suppression coefficient.

Sense-making corresponds to a stable limit cycle (σ_Φ > 0), which arises via a **subcritical** Hopf bifurcation when γ_eff ≥ γ_krit(τ, α, κ_eff).

---

## Reference Parameters

All Hopf bifurcation values are derived from the characteristic equation  
`det(iω*·I − J₀ − J₁·exp(−iω*τ)) = 0`  
at the fixed point with κ = 0.30, B = 0 (Table 2, manuscript).

### v4.4 Reference (τ = 3.8 s, manuscript reference)

| Parameter | Value | Description |
|---|---|---|
| τ | 3.8 s | Integration latency |
| α | 0.30 | Integration decay rate |
| γ_krit | **7.19** | Hopf bifurcation threshold (char. eq.) |
| ω* | **1.215 rad/s** | Hopf frequency (char. eq.) |
| T* | 5.171 s | Period at onset |
| sin(ω*τ) | **−0.996** | Phase condition (< 0 required) ✓ |
| l₁ | **+7.53** | First Lyapunov coeff. → **subcritical** |
| Re(dλ/dγ) | +0.022 | Transversality coefficient > 0 ✓ |

### SNR-optimal parameters (τ = 2.0 s)

| Parameter | Value | Description |
|---|---|---|
| τ | 2.0 s | Integration latency (SNR-optimal) |
| α | 0.30 | Integration decay rate |
| γ_krit | **15.02** | Hopf bifurcation threshold (char. eq.) |
| ω* | **2.031 rad/s** | Hopf frequency (char. eq.) |
| T* | 3.093 s | Period at onset |
| sin(ω*τ) | **−0.908** | Phase condition (< 0 required) ✓ |
| l₁ | **+17.10** | First Lyapunov coeff. → **subcritical** |

> **Note on earlier values:** The README previously reported γ_krit = 6.50 and ω* = 1.131 rad/s (τ = 3.8 s), and γ_krit = 13.20 and ω* = 1.864 rad/s (τ = 2.0 s). These were computed with an earlier version of the Hopf scanner that contained a numerical error. The v4.4 values above are derived from the characteristic equation directly and are confirmed by `reproduce_all.py`.

### ω*(τ) Curve (v4.4, from characteristic equation)

| τ (s) | γ_krit | ω* (rad/s) | f* (Hz) | sin(ω*τ) | HRV band |
|---|---|---|---|---|---|
| 3.0 | 8.12 | 1.349 | 0.215 | −0.963 | HF (Atem) |
| 3.5 | 7.52 | 1.270 | 0.202 | −0.986 | HF (Atem) |
| **3.8** | **7.19** | **1.215** | **0.193** | **−0.996** | HF/LF boundary |
| 5.0 | 6.45 | 0.990 | 0.158 | −0.980 | LF ✓ |
| 6.0 | 5.88 | 0.851 | 0.135 | −0.967 | LF ✓ |
| 8.0 | 5.12 | 0.668 | 0.106 | −0.950 | LF ✓ |
| 10.0 | 4.71 | 0.558 | 0.089 | −0.920 | LF ✓ |

For prediction P1 (HRV): physiologically relevant τ ∈ [6–10 s] → f* ∈ [0.09–0.14 Hz] (LF band, consistent with baroreflex resonance ~0.1 Hz).

---

## Margin Sensitivities (Result R3, corrected)

An earlier version claimed ∂M/∂Π = ∂M/∂B = 0 (exact). This was incorrect.  
Correct results (Propositions R3a and R3b, §3.3):

| Modulator | Channel | Result |
|---|---|---|
| Π (protentional field) | γ_eff only (additive) | **∂M/∂Π = λ_Π = 0.20 (exact)** |
| B (belonging, q_B=1) | γ_eff and γ_krit via κ_eff | **∂M/∂B ≈ +0.39 > 0 (numerical, B=1)** |
| B (belonging, q_B→0) | conformity dominates | ∂M/∂B may be negative |

---

## Six Empirical Predictions

| # | Claim | Primary measure | Falsification |
|---|---|---|---|
| P1 | Higher SAH → LF-band HRV coherence peak (0.04–0.15 Hz), narrower bandwidth | HRV spectrum + SAH score | No LF coherence increase; r_S < 0.20 |
| P2 | Lower EEG microstate α → lower γ_krit (longer dwell times) | EEG persistence + PIL/MLQ | No correlation |
| P3 | B and Π substitutable along bifurcation boundary (B_min(Π) decreasing) | Loneliness + Openness + MLQ | No threshold effect |
| P4 | H ≈ 0.75 → highest SAH (inverted U) | R/S analysis HRV + SAH | No optimum at H ≈ 0.75 |
| P5 | Sense-loss has retrospective (subcritical) onset | Interview + retrospective archive | Gradual onset found |
| P6 | Rising SAH variance precedes crisis (critical slowing down) | Weekly EMA 6 months | No early-warning signal |

---

## Usage Example

```python
from sinndynamik_solver_v4 import SinnDynamikSystem, HopfScanner

# Reproduce Hopf threshold (reference: τ=3.8)
scanner = HopfScanner()
hopf = scanner.find_hopf(tau=3.8)
print(f"gamma_krit = {hopf['gamma_krit']:.3f}")   # expect ~7.19
print(f"omega*     = {hopf['omega_star']:.4f} rad/s")  # expect ~1.215
print(f"sin(omega*tau) = {hopf['sin_phase']:.4f}")     # expect ~-0.996

# Run simulation above threshold
sys = SinnDynamikSystem()
result = sys.simulate(T=300.0, burn_in=100.0,
                      gamma_ind_0=9.0, tau=3.8)
print(f"sigma_Phi = {result['sigma_Phi']:.4f}")
print(f"In sense window: {result['in_sense_window']}")
```

---

## Phenomenological Grounding

The model connects formal results to phenomenological concepts (see §4.3 of the manuscript):

- **τ** ↔ Husserlian retention window; Heidegger's *Zeitlichkeit*
- **γ_eff** ↔ temporal carrying capacity; Heidegger's *Sorge*
- **Subcritical Hopf bifurcation** ↔ tipping-point quality of sense-making onset
- **Bistability / hysteresis** ↔ robustness of established sense-making; resistance to perturbation
- **X_int** ↔ integrable residual; the processable remainder
- **γ_koll** ↔ collective carrying capacity; Heidegger's *Mitsein*
- **q_B = 1 vs q_B → 0** ↔ resonant belonging (*Mitsein*) vs conformity pressure (*das Man*)

---

## Changes in v4.4

| Area | Change |
|---|---|
| Bifurcation type | Supercritical → **subcritical** (l₁ = +7.53, computed via Kuznetsov formula) |
| γ_krit (τ=3.8) | 6.50 → **7.19** (from characteristic equation) |
| ω* (τ=3.8) | 1.131 → **1.215 rad/s** (from characteristic equation) |
| sin(ω*τ) (τ=3.8) | −0.915 → **−0.996** |
| γ_krit (τ=2.0) | 13.20 → **15.02** (from characteristic equation) |
| ω* (τ=2.0) | 1.864 → **2.031 rad/s** |
| R3 (∂M/∂Π) | ∂M/∂Π = 0 (wrong) → **∂M/∂Π = λ_Π = 0.20 (exact)** |
| R3 (∂M/∂B) | ∂M/∂B = 0 (wrong) → **∂M/∂B ≈ +0.39 (composite)** |
| Article type | Original Research → **Hypothesis and Theory** |
| Supplement | Added `reproduce_all.py` (10/10 PASS) and `params_reference.json` |

---

## Citation

```
[Author(s)] (2026). Temporal Meaning Dynamics: A Delay-Differential Equation
Model of Sense-Making, Bifurcation, and Temporal Carrying Capacity.
Frontiers in Computational Neuroscience (Hypothesis and Theory). doi: [pending]
```

---

## License

MIT License — see LICENSE file.

---

*For questions or correspondence: clemenskluge(at)gmx.de*
