## Appendix A. Conventions and Computational Pipeline

本付録では、本論の計算と記述で用いる記号規約、および結果を再現する計算パイプラインを簡潔にまとめる。

### A.1 記号規約

- 空間 tensor 指数は $i,j=1,2,3$ とする。
- spin-2 基底ラベルは $A=1,\dots,5$ とする。
- spin-1 triplet ラベルは script と整合して $k=0,1,2$ と書く。
- 体積因子は ${\rm Vol}=2\pi^2Lr^3$ とする。
- 有効ポテンシャルは $V_{\rm eff}=-\left(\frac{R}{2\kappa^2}+\theta_{\rm NY}N\right)\,{\rm Vol}$ で定義する。
- Pontryagin 密度は $P=\langle R,\star R\rangle=\frac12\,\varepsilon^{abcd}R_{abef}R_{cd}{}^{ef}$ で定義する。

### A.2 Where the mode dictionary is fixed

mode dictionary 自体は本論と数値付録で役割分担して提示する。

1. §2.5 では、raw 幾何変数から最終基底への見取り図だけを与える。
2. §5 で spin-1 の physical / auxiliary 判定を確立する。
3. §6.4 で、本論として採用する最終辞書を確定する。
4. Appendix H.3 では、その辞書に対応する数値アンカーを参照用にまとめる。

### A.3 Computational pipeline

本論の数値・記号計算は、概念的には以下の順序で構成される。

1. homogeneous frame と構造定数 $C^i{}_{jk}$ を与える。
2. Einstein-Cartan 接続から curvature、Ricci scalar、Nieh-Yan 密度、Pontryagin 密度を計算する。
3. $V_{\rm eff}$ を組み立て、等方点と変形点での Hessian / generalized eigenvalue を評価する。
4. full shear basis $T_A$、spin-1 physical/auxiliary basis $(X,Y)$、torsion mode を最終基底に整理する。
5. 追加の consistency check として cross-Hessian、homogeneous consistency、不変量生成を独立 script で検証する。

### A.4 Scripts

主要スクリプトを用途別に整理する。主要な Appendix (B, D, E, F, G, H) の末尾にも対応するスクリプト参照を記載している。

**Pontryagin 保護（Appendix B 対応）**

- `proofs/riemann_block_classification.py`
- `proofs/pontryagin_block_decomposition.py`
- `proofs/plucker_cancellation_proof.py`
- `proofs/pfaffian_decomposability.py`
- `proofs/vt_mode_cancellation.py`
- `proofs/analytic_proof_P_equals_zero.py`
- `paper03/pontryagin_scan_s3.py`
- `paper03/pontryagin_scan_nil3.py`
- `paper03/pontryagin_anisotropic_twist.py`
- `paper03/torsion_mode_analysis.py`

**spin-2 sector（Appendix D 対応）**

- `paper03/spin2_schur_degeneracy.py`
- `paper03/spin2_spectrum_splitting.py`
- `paper03/spin2_diagonal_hessians.py`
- `paper03/spin2_offdiag_hessian.py`

**spin-1 Schur reduction and auxiliary analysis（§5, Appendix E, Appendix I 対応）**

- `paper03/spin1_symbolic_firstprinciples.py`
- `paper03/spin1_kinetic_metric.py`
- `paper03/spin1_null_mode_analysis.py`
- `paper03/spin1_auxiliary_exclusion.py`
- `paper03/spin1_appI_projection_artifact.py`

**cross-Hessian と最終基底整合性（Appendix F, G 対応）**

- `paper03/spin1_cross_hessian.py`
- `paper03/homogeneous_consistency.py`
- `paper03/eft_invariant_generation.py`

**Scalar/torsion sector と数値表（Appendix H 対応）**

- `paper03/spin0_torsion_cross_hessian.py`
- `paper03/spin0_reta_block_eigenvalues.py`
- `paper03/spin0_rLetaV_hessian.py`

### A.5 Sanity checks

本論で本文に用いる結論は、少なくとも以下の sanity check を通過している。

1. raw spin-1 二次係数 $(\mu_\omega^2,\mu_\delta^2,m_{\omega\delta})$ の再現。
2. 全 21 cross-Hessian が 0 であること。
3. homogeneous Maurer-Cartan consistency。
4. invariant basis $\{I_2,I_3,I_4,W_2\}$ の exact check。

これらは本論の主要命題を補助するための再現性情報である。
