# Calculation code and reproduction

[日本語](README_ja.md) / [Conventions](docs/CONVENTIONS.md) / [Implementation guide](docs/IMPLEMENTATION.md)

For Paper03 v3 (2026-11-01). Python 3.10 or later.
EC denotes Einstein–Cartan; NY denotes Nieh–Yan. The model is spatially
homogeneous and Euclidean, with constant complete NY coupling and alpha=0.
The manuscript is supplied as [LaTeX source](../LaTeX/main.tex).

## Reproduce the figures and data

Run from the `script` directory:

```powershell
python -m pip install -r requirements.txt
python -B -m paper_results render --output-dir ../rendered
python -B -m paper_results recompute --output-dir ../recomputed-data
python -B -m paper_results render --data-dir ../recomputed-data --output-dir ../recomputed-figures
```

Each output directory must be new. `render` reads [data/revised](../data/revised),
verifies the eleven dataset hashes in [manifest.json](../data/revised/manifest.json),
and produces ten PNG/PDF pairs. It does not call the geometric calculation.
`recompute` evaluates the fixed grids and writes eleven data files and a manifest.
The supplied figure pairs are in [LaTeX/figures/v3](../LaTeX/figures/v3).

The [requirements file](requirements.txt) lists NumPy and Matplotlib for
calculation and rendering, and IPython and ipywidgets for notebook display.
The supplement consists of `paper_results`, the figure notebook and the saved-data viewer.

## Computation interfaces

| Entry point | Result | Manuscript |
|---|---|---|
| [model.cartan](paper_results/model.py) | Coframe, connection, full EC curvature, NY, P and B | Section 2; Appendices A–B |
| [model.potential](paper_results/model.py) | Closed static potential U | (M3) |
| [model.scalar_data / scalar_bulk](paper_results/model.py) | Scalar gradient/Hessian and round scalar circle density | Section 3 |
| [model.twist_bulk](paper_results/model.py) | Static twist cost, quartic derivative term and velocity Hessian | Section 4 |
| [model.shear_curvatures / shape_minimum](paper_results/model.py) | Path curvatures and constrained doublet minimum | Section 5 |
| [circle.joint_bulk](paper_results/circle.py) | Joint bulk density and static, quadratic and linear contributions | (M4)–(M6); Appendix J |
| [circle.diagonal_cross_coefficients](paper_results/circle.py) | Five velocity cross coefficients | (X2) |
| [circle.coupled_periodic_family](paper_results/circle.py) | Fixed quadrature for the periodic shape–twist sequence | (X4) |
| [data / figures](paper_results/data.py) / [renderer](paper_results/figures.py) | Saved arrays, manifest and figure generation | Appendix H |

`eta,V` in `paper_results` are geometric torsion coefficients. `joint_bulk`
accepts instantaneous values and first derivatives of r, L, M and omega;
M has determinant one and `tr(inv(M) @ M_prime)=0`. Its density has periodic
total derivatives removed. Static minimization and variation over periodic
profiles have the separate domains specified in the manuscript.

## Notebook and viewer

Open [the figure notebook](DPPUv2_Paper03_Figures_v3.ipynb) in a Jupyter environment,
with the working directory at the repository root or `script`. Its default
`RECOMPUTE=False` renders the saved data; choose `True` to recompute explicitly.
Each run writes to new time-stamped folders next to `script`
(`notebook-figures-YYYYMMDD-HHMMSS`, and `notebook-data-...` when recomputing).
The notebook also loads the [saved-data viewer](paper03_saved_data_viewer.py) in the same directory.
The viewer's `show(data_dir=None)` returns a dropdown and output widget for eight
saved-data views, including static slices, velocity cross coefficients and the
periodic action sequence. It does not solve new geometries or vary the dataset grid.
Jupyter itself must be available in the environment; the widget uses IPython and ipywidgets.

The conventions are in [CONVENTIONS](docs/CONVENTIONS.md). APIs and return values are described in the [implementation guide](docs/IMPLEMENTATION.md).
