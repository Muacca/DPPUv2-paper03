# Spin-0, Spin-1 and Spin-2 Deformations of S³×S¹ in Einstein–Cartan–Nieh–Yan Theory

Muacca. First release: v1, 2026-03-21. Previous revision: v2, 2026-06-20.
Current revision: **v3, 2026-11-01**.

The paper organizes static potentials and periodic circle actions by spatial
spin-0, spin-1 and spin-2 representations. It derives scalar couplings, a
conditional radius minimum, positive restricted twist and shape derivative
terms, and a continuous bifurcation of constrained static shape minima.
At fixed scalar controls, minimizing the static potential over constant
volume-preserving shapes and twist selects the round, zero-twist geometry.

An exact joint circle action exhibits nonlinear derivative cross terms and
shift couplings. A smooth coupled shape-twist family makes this Euclidean
functional unbounded below even at fixed scalar controls. The static minimum
and the circle functional therefore have distinct variational properties.
The model has alpha=0 and constant complete Nieh-Yan coupling.

- [English PDF](DPPUv2-paper03.pdf)
- [English LaTeX](LaTeX/main.tex)
- [Japanese manuscript](DPPUv2-paper03_sec01.md)
- [Computation and figure instructions](script/README.md)
- [Saved data and manifest](data/revised/manifest.json)
- [Figure notebook](script/DPPUv2_Paper03_Figures_v3.ipynb)
- [Revision history](CHANGELOG.md)
- [Japanese README](README_ja.md)

Install `script/requirements.txt`, then run from `script`:

    python -B -m paper_results render --output-dir ../rendered
    python -B -m paper_results recompute --output-dir ../recomputed-data
    python -B -m paper_results render --data-dir ../recomputed-data --output-dir ../recomputed-figures

Use new output directories. The first command only renders supplied data;
recalculation is explicit. The notebook defaults to saved-data rendering.
The saved-data viewer is described in the computation instructions.
The code supplement consists of paper_results, the figure notebook and the saved-data viewer.
[Conventions](script/docs/CONVENTIONS.md) and an [implementation guide](script/docs/IMPLEMENTATION.md) are included.

The associated archive is [Zenodo](https://doi.org/10.5281/zenodo.19144481).
The supplied datasets are described in [their manifest](data/revised/manifest.json).
Earlier archival versions may use different definitions.
Related work: [Paper01](https://github.com/Muacca/DPPUv2-paper01) and
[Paper02](https://github.com/Muacca/DPPUv2-paper02).
Code and data use the [MIT License](LICENSE.txt).
Contact: muacca@dmwp.jp or [Issues](https://github.com/Muacca/DPPUv2-paper03/issues).

The joint circle-action derivation is Appendix J of the [PDF](DPPUv2-paper03.pdf)
([LaTeX](LaTeX/appendices/appJ.tex); [Japanese](DPPUv2-paper03_appJ.md)).
