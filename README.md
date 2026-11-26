# Constructive Lean 4 Kernel and Exact Zero-Float Certificates for OpenAI's Superpolynomial Simple Polytope Diameter

**Author:** Lluis Eriksson  
**Protocol:** `HIRSCH-DIAM-CERT-1.0`  
**Licenses:** MIT (code, Lean 4 proofs, certificates) / CC-BY-4.0 (manuscript & figures)

## Overview

Independent, explicitly scoped, zero-floating-point two-layer verification companion for the discrete convex-corner bend, prescribed dimension/facet product padding, Klee–Walkup wedge lifting to complementary $2d$-facet vertices, Schlegel projective prismatoid / Minkowski-sum middle-section lifting, and two-stage symmetric positive doubly-stochastic packing polytope mechanisms of OpenAI's November 24, 2026 preprint *"Superpolynomial graph diameter of simple polytopes"*.

- **Layer 1 (`lean/OpenAIHirschDiameter.lean`):** 18 standalone constructive Lean 4.35.0 sub-lemmas (`0` `sorry`, `0` custom `axiom`, `0` `Classical.choice`, `0` warnings; `1 none`, `1 [propext]`, `16 [propext, Quot.sound]`).
- **Layer 2 (`src/hirsch_exact.py`, `src/hirsch_certifier.py`):** Pure-stdlib Python `int`/`Fraction` zero-float auditor across 4 modes (`convex_corner_bend_and_crossing_cert`, `klee_walkup_wedge_and_complement_cert`, `schlegel_projective_and_prismatoid_cert`, `symmetric_doubly_stochastic_packing_cert`).
- **Canonical Atlas & Red-Team Audit:** 16 canonical certificates verified (`0.092–19.354 ms`), `82/82` pytest checks (`76/76` live adversarial forgeries rejected with exact branch regexes), and full 41-file SHA-256 integrity manifest.

## Reproducibility (5-Gate Verification)

```bash
python scripts/run_all_experiments.py
python scripts/make_figures.py
python scripts/verify_all.py
```
