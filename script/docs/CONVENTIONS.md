# Geometric and computational conventions

[日本語](CONVENTIONS_ja.md) / [Reproduction](../README.md) / [Implementation](IMPLEMENTATION.md)

The definitions are given in manuscript Section 2 and Appendices A, B, I and J.
The attached sources are [Section 2](../../LaTeX/sections/sec2.tex),
[Appendix A](../../LaTeX/appendices/appA.tex),
[Appendix B](../../LaTeX/appendices/appB.tex),
[Appendix I](../../LaTeX/appendices/appI.tex) and
[Appendix J](../../LaTeX/appendices/appJ.tex).

## 1. Model and global geometry

The calculation uses spatially homogeneous Euclidean EC with constant complete
NY coupling and no curvature-squared term. Spatial indices are 0,1,2; the
circle index is 3; the orientation is epsilon0123=+1.
The circle coordinate y has period one. At constant L its circumference is L;
for varying profiles the length at a fixed spatial point is the integral of L(y).
The round sphere radius is r/2 and the spatially integrated four-volume density is
$\pi^2 Lr^3/4$. M is a real 3 by 3 matrix with determinant one.

## 2. Component arrays and static curvature

`C[a,b,c]`, `Gamma[a,b,c]` and `R[a,b,c,d]` store
$`C^a{}_{bc}`$, $`\Gamma^a{}_{bc}`$ and $`R^a{}_{bcd}`$.
The last connection index is the differentiation direction.
Frame indices are raised and lowered with the Euclidean identity metric.

$$
\begin{aligned}
de^a&=\frac12 C^a{}_{bc}e^b\wedge e^c,\qquad
\Gamma^a{}_{b}=\Gamma^a{}_{bc}e^c,\\
T^a{}_{bc}&=C^a{}_{bc}+\Gamma^a{}_{cb}-\Gamma^a{}_{bc},\\
\Gamma^{\mathrm{LC}}_{abc}&=\frac12(C_{abc}+C_{cba}-C_{bac}),\\
K_{abc}[T]&=-\frac12(T_{abc}+T_{bca}-T_{cab}),\\
R^a{}_{bcd}&=\Gamma^a{}_{ec}\Gamma^e{}_{bd}
-\Gamma^a{}_{ed}\Gamma^e{}_{bc}+\Gamma^a{}_{be}C^e{}_{cd}.
\end{aligned}
$$

The last formula assumes constant frame components. Circle-dependent
connection derivatives are treated by the separate bulk formula (M6).
The lowered EC curvature is antisymmetric within each pair; pair exchange is
not assumed. For a constant coframe change e=A ebar, transform the structure
constants with A on the output and the actual inverse on each input:
`A[a,d] inv[e,b] inv[f,c] Cbar[d,e,f]`.

## 3. Torsion input

The manuscript and `paper_results` use geometric T:
$`T_{ijk}=(2\eta/r)\varepsilon_{ijk}`$ and
$`T^{(V)}_{abc}=(\delta_{ac}V_b-\delta_{ab}V_c)/3`$, with
$`V_b=V\delta_{b3}`$.
AX sets V=0; VT sets eta=0; MX allows both. All include zero amplitudes.
For comparison, the input array in the code accompanying Papers01/02 is
tau=-T; both input amplitudes change sign:
`eta_code=-eta` and `V_code=-V`.
`cartan` returns geometric T directly. Its `torsion_residual` checks the
Cartan reconstruction including C.

## 4. NY and the two curvature contractions

$$
\begin{aligned}
N&=TT-REE
=\frac14\varepsilon^{abcd}T^e{}_{ab}T_{ecd}
-\frac12\varepsilon^{abcd}R_{abcd},\\
P&=\frac14\varepsilon^{cdef}R_{abcd}R_{abef},\\
B&=\sum_{c<d}(R_{01cd}R_{23cd}-R_{02cd}R_{13cd}+R_{03cd}R_{12cd}).
\end{aligned}
$$

N is a scalar density; N times the volume form is the NY four-form.
The standard form-index Pontryagin density P differs from the internal
contraction B. The full component Hodge inner product is 2P.
`cartan` returns `TT, REE, NY, P, B` with these normalizations.
On the constant unimodular family N=P=0; smooth periodic profiles make the
complete NY integral vanish even when its local density is nonzero.

## 5. Action sign and variation domains

The paper defines $\mathcal S=-\int\mathrm{vol}_4(R/(2\kappa^2)+\theta N)$.
At constant profiles $\mathcal S=U$. This is the sign of U=-S in Papers01/02.
The whole positive-integral action S in those papers has the opposite sign.

Static minimization over shape and twist at fixed scalars selects the round,
zero-twist geometry. With nonzero twist fixed, the static minimizing shape can
shift. Variation over smooth periodic profiles is a different problem: the
coupled shape–twist sequence is unbounded below without a frequency cutoff.
`joint_bulk` returns the periodic bulk density, with total derivatives removed;
it includes the anisotropic shift and velocity-linear term.
Its inputs obey `tr(inv(M) @ M_prime)=0`.

## 6. Shape paths and restricted functions

The shear basis has $`\mathrm{tr}(T_AT_B)=\delta_{AB}`$.
`shear_curvatures` uses `M=M0 exp(q T_A)` by default; its
`linear_normalized` option denotes a different path. At a nonstationary point,
second derivatives depend on that path.
The coordinate relation to Paper02 at static zero twist is
$x=a^{-2}=(1+\epsilon)^{-2/3}$.

`scalar_bulk` has round shape and zero twist, with periodic scalar profiles.
`twist_bulk` has round shape and fixed scalar controls. Its quartic term is
fourth order in amplitude and second order in circle derivatives.
The reported Hessians describe these specified geometric variations.
