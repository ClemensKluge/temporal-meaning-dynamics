# Temporal Meaning Dynamics

**Code repository for:**  
*Temporal Meaning Dynamics: A Delay-Differential Equation Model of Sense-Making, Bifurcation, and Temporal Carrying Capacity*  
Frontiers in Computational Neuroscience (submitted 2026)

---

## Overview

This repository contains the numerical implementation of the Temporal Meaning Dynamics (TMD) model — a delay-differential equation (DDE) system in which sense-making emerges as a supercritical Hopf bifurcation governed by temporal carrying capacity and phase-matched feedback.

The model formalises the phenomenological claim that meaning-constitution requires a minimum capacity to sustain temporal openness (γ_eff ≥ γ_krit), and makes six empirically testable predictions (P1–P7) linking HRV spectral structure, EEG microstates, and narrative coherence to the bifurcation threshold.

---

## Quick Start

```bash
# Install dependencies
pip install -r requirements.txt

# Run all figures (after simulations complete)
python generate_figures.py

# Run a specific figure
python generate_figures.py --fig 1
```

---

## Files

| File | Description |
|---|---|
| `sinndynamik_solver_v4.py` | Core DDE solver — all model equations (Eq. 1–14), Hopf scanner, sensitivity analysis |
| `generate_figures.py` | Reproduces Figures 1–4 from the manuscript |
| `requirements.txt` | Python dependencies |

---

## Model Equations

The core system is a 2-variable DDE with delayed feedback:

```
dΦ/dt  =  I_aktiv · η(u, C)  −  α · Φ                           (Eq. 2)

du/dt  =  −β_p · (1 − p_ref · exp(−u))
          −  κ_eff · X_int
          +  γ_eff · [Φ(t−τ) − Φ(t)]                           (Eq. 4)
```

where Φ is the integration state, u = log(p) is log-precision, τ is the integration latency, and γ_eff is the effective temporal carrying capacity.

Sense-making corresponds to a stable limit cycle (σ_Φ > 0), which arises via Hopf bifurcation when γ_eff ≥ γ_krit(τ, α, H).

---

## Reference Parameters

### v3 Reference (manuscript, τ = 3.8 s)

| Parameter | Value | Description |
|---|---|---|
| τ | 3.8 s | Integration latency |
| α | 0.30 | Integration decay rate |
| γ_krit (γ_eff) | 6.50 | Hopf bifurcation threshold |
| ω* | 1.131 rad/s | Hopf frequency (oscillation onset) |
| T* | 5.555 s | Period at onset |
| C(τ, ω*) | 1.675 | Contrast gain |
| sin(ω*τ) | −0.915 | Phase condition (< 0 required) ✓ |

### v4 Reference (SNR-optimal, τ = 2.0 s)

| Parameter | Value | Description |
|---|---|---|
| τ | 2.0 s | Integration latency (SNR-optimal) |
| α | 0.30 | Integration decay rate |
| γ_krit (γ_eff) | 13.20 | Hopf bifurcation threshold |
| ω* | 1.864 rad/s | Hopf frequency (oscillation onset) |
| T* | 3.371 s | Period at onset |
| C(τ, ω*) | 1.915 | Contrast gain |
| sin(ω*τ) | −0.553 | Phase condition (< 0 required) ✓ |

**Note on ω*:** ω* is τ-dependent. The value ω* = 1.643 rad/s reported in manuscript v3 is superseded by the empirically measured values above. See `ω*(τ) table` in Section 3.1.

### ω*(τ) Curve

| τ (s) | γ_krit (eff) | ω* (rad/s) | T* (s) | C(τ,ω*) | sin(ω*τ) |
|---|---|---|---|---|---|
| 1.5 | 19.95 | 2.388 | 2.631 | 1.952 | −0.426 |
| 1.8 | 15.45 | 2.061 | 3.049 | 1.920 | −0.538 |
| **2.0** | **13.20** | **1.864** | **3.371** | **1.915** | **−0.553** |
| 2.3 | 11.45 | 1.697 | 3.703 | 1.857 | −0.690 |
| 2.7 | 9.45 | 1.458 | 4.309 | 1.844 | −0.714 |
| 3.0 | 8.45 | 1.332 | 4.717 | 1.820 | −0.754 |
| 3.4 | 7.45 | 1.206 | 5.210 | 1.775 | −0.819 |
| **3.8** | **6.50** | **1.131** | **5.555** | **1.675** | **−0.915** |
| 4.5 | 6.45 | 0.980 | 6.411 | 1.611 | −0.955 |
| 5.0 | 6.45 | 0.905 | 6.943 | 1.540 | −0.983 |

---

## Six Empirical Predictions

| # | Claim | Primary measure | Falsification |
|---|---|---|---|
| P1 | Higher SAH → HRV peak near ω* ≈ 0.19–0.22 Hz | HRV spectrum + SAH score | Peak absent or uncorrelated |
| P2 | Lower EEG microstate α → lower γ_krit | EEG persistence + PIL/MLQ | No correlation |
| P3 | B and Π substitutable above B_floor | Loneliness + Openness + MLQ | No threshold effect |
| P4 | H ≈ 0.75 → highest SAH (inverted U) | R/S analysis HRV + SAH | No optimum |
| P5 | Sense-loss has retrospective onset | Interview + retrospective archive | Gradual onset |
| P6 | Rising SAH variance precedes crisis | Weekly EMA 6 months | No early-warning signal |
| P7 | HRV dominant period near τ_SNR ≈ 2.0 s | HRV spectrum + SAH | Peak not at 2.0 s |

---

## Usage Example

```python
from sinndynamik_solver_v4 import SinnDynamikSystem, HopfScanner

# Run a simulation at reference parameters (tau=2.0)
sys = SinnDynamikSystem()
result = sys.simulate(T=200.0, burn_in=80.0,
                      gamma_ind_0=15.0, gamma_koll_0=0.0)
print(f"sigma_Phi = {result['sigma_Phi']:.4f}")
print(f"In sense window: {result['in_sense_window']}")

# Scan bifurcation threshold across tau values
scanner = HopfScanner()
hopf = scanner.find_hopf(tau=2.0)
print(f"gamma_krit = {hopf['gamma_krit']:.3f}")
print(f"omega*     = {hopf['omega_star']:.4f} rad/s")
print(f"C(tau,w*)  = {hopf['contrast_gain']:.4f}")

# Contrast gain and prediction weight
from sinndynamik_solver_v4 import contrast_gain, r_tau
print(f"C(2.0, 1.864) = {contrast_gain(2.0, 1.864):.4f}")
print(f"r(2.0, 1.864) = {r_tau(2.0, 1.864):.4f}")
```

---

## Phenomenological Grounding

The model connects formal results to phenomenological concepts:

- **τ** ↔ Husserlian retention window; Heidegger's *Zeitlichkeit*
- **γ_eff** ↔ temporal carrying capacity; Heidegger's *Sorge*
- **Hopf bifurcation** ↔ threshold of metastable openness
- **X_int** ↔ integrable residual; the processable remainder
- **C(τ,ω*)** ↔ contrast gain; the Mikrosakkaden analogy — sense arises not from maximal difference but from meaningful residual after prediction
- **γ_koll** ↔ collective carrying capacity; Heidegger's *Mitsein*

---

## Citation

```
[Author(s)] (2026). Temporal Meaning Dynamics: A Delay-Differential Equation
Model of Sense-Making, Bifurcation, and Temporal Carrying Capacity.
Frontiers in Computational Neuroscience. doi: [pending]
```

---

## License

MIT License — see LICENSE file.

---

*For questions or correspondence: [author email]*
