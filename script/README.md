# DPPUv2 Computation Engine v5 — Scripts Directory

⇒ [日本語](README_ja.md)

**Paper**: "Unified Geometric Landau EFT of Homogeneous S³×S¹ Minisuperspace in Einstein-Cartan + Nieh-Yan Theory" (paper03)

Python packages and execution scripts for numerical and symbolic computation in homogeneous S³×S¹ minisuperspace with Einstein-Cartan + Nieh-Yan theory.

---

## Directory Structure

```
script/
├── dppu/                      # Main Python package (DPPUv2 Engine v5)
│   ├── geometry/              # Metric, volume form, structure constants
│   ├── connection/            # Levi-Civita connection, Contortion, EC connection
│   ├── curvature/             # Riemann, Ricci, Hodge dual, self-duality, Pontryagin, Weyl
│   ├── torsion/               # Torsion modes, Ansatz, Nieh-Yan density
│   ├── action/                # Lagrangian, effective potential, stability classification
│   ├── topology/              # Unified engine (S³×S¹, T³×S¹, Nil³×S¹)
│   ├── engine/                # Computation pipeline, logging, checkpointing
│   ├── scanning/              # Parameter scan, squashing scan, Phase 1 integration
│   └── utils/                 # Common utilities (Levi-Civita symbol, symbolic computation, visualization)
│
└── scripts/                   # Execution scripts
    ├── paper03/               # paper03-specific analysis scripts
    ├── proofs/                # Analytic and symbolic proof scripts
    └── visualize/             # Interactive viewer
```

### `docs/` - Documentation

Technical documentation and conventions:
- [DPPUv2 Engine CONVENTIONS](docs/CONVENTIONS.md) — Engine core conventions and specifications
- [DPPUv2 SymPy guideline](docs/SymPy_guideline.md) — SymPy usage guidelines and best practices

---

## Package Overview (dppu/)

| Module | Role | Key Classes / Functions |
|--------|------|------------------------|
| [`geometry`](dppu/geometry/README.md) | Metric and frame field definitions | `build_metric`, `frame_field` |
| [`connection`](dppu/connection/README.md) | EC connection construction | `levi_civita`, `contortion`, `ec_connection` |
| [`curvature`](dppu/curvature/README.md) | Curvature tensors, SD diagnostics, Pontryagin | `RiemannTensor`, `SDExtensionMixin`, `compute_pontryagin_inner_product` |
| [`torsion`](dppu/torsion/README.md) | Torsion structure | `Mode`, `NyVariant`, `build_torsion_tensor` |
| [`action`](dppu/action/README.md) | Action and stability analysis | `build_lagrangian`, `classify_stability` |
| [`topology`](dppu/topology/README.md) | Unified engine for all topologies and DOFs | `UnifiedEngine`, `DOFConfig`, `TopologyType`, `FiberMode` |
| [`engine`](dppu/engine/README.md) | 15-step computation pipeline | `BaseFrameEngine`, `ComputationLogger`, `CheckpointManager` |
| [`scanning`](dppu/scanning/README.md) | Parameter space scanning | `run_scan`, `Phase1ResultsLoader`, `SDDiagnosticsWithPhase1` |
| [`utils`](dppu/utils/README.md) | Common utilities | `LeviCivita`, `symbolic_simplify`, `visualize_potential` |

---

## Script Overview (scripts/)

### paper03/ — paper03-Specific Analysis

#### Pontryagin Protection

| Script | Description |
|--------|-------------|
| `pontryagin_scan_s3.py` | Pontryagin scan for S³ over torsion mode / parameter space |
| `pontryagin_scan_nil3.py` | Pontryagin scan for Nil³ |
| `pontryagin_anisotropic_twist.py` | Pontryagin check under anisotropic twist |
| `torsion_mode_analysis.py` | AX / VT / MX mode Pontryagin classification |

#### Spin-2 Sector

| Script | Description |
|--------|-------------|
| `spin2_schur_degeneracy.py` | Schur 5-fold degeneracy μ²_q = 48π² at isotropic point |
| `spin2_spectrum_splitting.py` | SO(3)→SO(2) splitting and phase transitions ε_c± |
| `spin2_diagonal_hessians.py` | H₁₁^q(ε), H₂₂^q(ε) via finite differences |
| `spin2_offdiag_hessian.py` | H₄₄^q (T₄, T₅ off-diagonal shear) |

#### Spin-1 Sector

| Script | Description |
|--------|-------------|
| `spin1_kinetic_metric.py` | Field-space metric G via velocity-mode engine (~300–600 s) |
| `spin1_null_mode_analysis.py` | Null mode classification (Case B: non-dynamical auxiliary) |
| `spin1_cross_hessian.py` | Cross-Hessian: spin-1 vs spin-2 / spin-0 decoupling |
| `spin1_auxiliary_exclusion.py` | Case C exclusion verification |

#### Spin-0 Sector

| Script | Description |
|--------|-------------|
| `spin0_rLetaV_hessian.py` | Spin-0 Hessian (r, η, V block) |
| `spin0_reta_block_eigenvalues.py` | Spin-0 (r, η) block eigenvalues |
| `spin0_torsion_cross_hessian.py` | Cross-Hessian: spin-0 torsion sector |

#### EFT

| Script | Description |
|--------|-------------|
| `eft_invariant_generation.py` | EFT invariant basis generation (Landau EFT selection rules) |
| `homogeneous_consistency.py` | Homogeneous geometry consistency check |

### proofs/ — Analytic and Symbolic Proofs

| Script | Proof Content |
|--------|---------------|
| `analytic_proof_P_equals_zero.py` | Algebraic proof of P = ⟨R, *R⟩ ≡ 0 for AX/VT modes |
| `symbolic_block_all_topologies.py` | Symbolic block-structure check for all topologies (S³, T³, Nil³) |
| `plucker_cancellation_proof.py` | **Pontryagin Protection (AX)**: Plücker-type cancellation proof |
| `pfaffian_decomposability.py` | **Pontryagin Protection (VT)**: Pfaffian decomposability proof |
| `pontryagin_block_decomposition.py` | Pontryagin block decomposition (P_SS, P_S3, P_33 vanish independently) |
| `riemann_block_classification.py` | Riemann block classification (R_{i3j3} ≠ 0 verification) |
| `vt_mode_cancellation.py` | VT mode cancellation mechanism |
| `analytical_proof.py` | Lemma 1: C²(ε=0)=0, ∂C²/∂ε|₀=0, ∂²C²/∂ε²|₀>0 (SymPy symbolic) |

### visualize/ — Interactive Viewer

| File | Contents |
|------|----------|
| `DPPUv2_interactive_viewer_p03.py` | Interactive potential viewer for paper03 (Jupyter, `%matplotlib widget` required) |

---

## Quick Start

### Install Dependencies

```bash
pip install sympy numpy scipy pandas tqdm matplotlib
```

### Run Engine with UnifiedEngine

```python
from dppu.topology.unified import UnifiedEngine, DOFConfig, TopologyType, FiberMode
from dppu.torsion.mode import Mode
from dppu.torsion.nieh_yan import NyVariant

# S³×S¹ with squash (equivalent to paper02 S3S1Engine):
cfg = DOFConfig(
    topology=TopologyType.S3,
    enable_squash=True,
    torsion_mode=Mode.AX,
    ny_variant=NyVariant.FULL,
)
engine = UnifiedEngine(cfg)
engine.run()

Veff = engine.data['potential']
fp   = engine.get_free_params()   # dict of active SymPy Symbols
```

### Legacy engine preset names

```python
cfg = DOFConfig.from_engine('S3S1Engine')    # S³×S¹ squash
cfg = DOFConfig.from_engine('T3S1Engine')    # T³×S¹
cfg = DOFConfig.from_engine('Nil3S1Engine')  # Nil³×S¹ squash
cfg.torsion_mode = Mode.MX
engine = UnifiedEngine(cfg)
engine.run()
```

### Pontryagin Protection Proofs

```bash
# Algebraic proof of P = 0 for AX/VT modes
python scripts/proofs/analytic_proof_P_equals_zero.py

# Plücker cancellation (AX mode protection mechanism)
python scripts/proofs/plucker_cancellation_proof.py

# Pfaffian decomposability (VT mode protection mechanism)
python scripts/proofs/pfaffian_decomposability.py

# Symbolic block-structure verification for all topologies
python scripts/proofs/symbolic_block_all_topologies.py
```

### Spin-2 Analysis

```bash
# Schur degeneracy at isotropic point
python scripts/paper03/spin2_schur_degeneracy.py

# SO(3)→SO(2) splitting and phase transitions ε_c±
python scripts/paper03/spin2_spectrum_splitting.py

# H₁₁^q(ε), H₂₂^q(ε) via finite differences
python scripts/paper03/spin2_diagonal_hessians.py
```

### Spin-1 Analysis

```bash
# Field-space metric G (velocity mode, ~300–600 s)
python scripts/paper03/spin1_kinetic_metric.py

# Null mode classification
python scripts/paper03/spin1_null_mode_analysis.py

# Cross-Hessian decoupling check
python scripts/paper03/spin1_cross_hessian.py
```

---

## paper03 Main Results and Corresponding Scripts

| Result | Paper Section | Content (Proof Type) | Script |
|--------|--------------|----------------------|--------|
| **Pontryagin Protection Theorem** | §3 | P ≡ 0 for AX/VT — **analytic proof** (Plücker/Pfaffian cancellation) | `proofs/plucker_cancellation_proof.py`, `proofs/pfaffian_decomposability.py` |
| **P ≠ 0 for MX mode** | §3 | P₀ = 2Vη(V²r²+9η²−36)/(9r³); source is torsion mixing Vη, not twist | `paper03/torsion_mode_analysis.py` |
| **Schur 5-fold degeneracy** | §4 | H_AB = μ²_q δ_AB, μ²_q = 48π² at isotropic point — **symbolic proof** | `paper03/spin2_schur_degeneracy.py` |
| **Phase transitions ε_c±** | §4 | ε_c+ ≈ 0.483 (squash), ε_c− ≈ −0.293 (shear) — **numerical analysis** | `paper03/spin2_spectrum_splitting.py`, `paper03/spin2_diagonal_hessians.py` |
| **Spin-1 physical/auxiliary decomp.** | §5 | X_k physical, Y_k auxiliary (Case B confirmed) — **analytic + numerical** | `paper03/spin1_kinetic_metric.py`, `paper03/spin1_null_mode_analysis.py` |
| **λ_phys ≈ −0.4183 < 0** | §5 | Schur-reduced finite generalized eigenvalue: tachyonic but non-ghost spin-1 mode | `paper03/spin1_kinetic_metric.py` |
| **Quadratic sector decoupling** | §6 | Cross-Hessians vanish at isotropic point — **numerical verification** | `paper03/spin1_cross_hessian.py`, `paper03/spin0_torsion_cross_hessian.py` |
| **EFT invariant basis** | §6 | Landau EFT selection rules and forbidden terms — **group-theoretic** | `paper03/eft_invariant_generation.py` |

---

## Topology Reference

| Topology | Lie Group | Background Curvature | Koszul Formula |
|----------|-----------|----------------------|----------------|
| S³×S¹ | SU(2) | +24/r² | Simplified (bi-invariant) |
| T³×S¹ | U(1)³ | 0 (flat) | Trivial |
| Nil³×S¹ | Heisenberg group | −1/(2R²) | General (bi-invariant does not hold) |

---

## Checkpoint Feature

Long-running computations can be interrupted and resumed:

```bash
# Run with checkpointing (--checkpoint-dir is independent of --output-dir)
python scripts/paper03/spin1_kinetic_metric.py \
    --output-dir output/ --checkpoint-dir ./checkpoints

# Resume after interruption (simply re-specify the checkpoint directory)
python scripts/paper03/spin1_kinetic_metric.py \
    --output-dir output/ --checkpoint-dir ./checkpoints
```

---

## License

See LICENSE file in the repository root.

---

**Author**: Muacca
**Version**: DPPUv2 Engine v5
**Date**: 2026-03
