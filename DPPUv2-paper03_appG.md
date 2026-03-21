## Appendix G. Landau-EFT Selection Rules and Forbidden Terms

本付録では、paper03 で採用する最終不変量基底と、その上での selection rule を簡潔にまとめる。

### G.1 Final invariant basis

採用する低次不変量は

$$
I_2=\mathrm{tr}(h^2),\qquad
I_3=h_{ij}\omega_i\omega_j,\qquad
I_4=\mathrm{tr}(h^3),\qquad
W_2=\omega_i\omega_i
$$

である。quartic helper として

$$
I_2^2,\qquad I_2W_2,\qquad W_2^2
$$

を用いる。

### G.2 Explicit form of $I_3$

最終 $T_A$ 規約では

$$
I_3=
\frac{q_1}{\sqrt6}\left(\omega_1^2+\omega_2^2-2\omega_3^2\right)+
\frac{q_2}{\sqrt2}\left(\omega_1^2-\omega_2^2\right)+
\sqrt2\,q_3\omega_1\omega_2+
\sqrt2\,q_4\omega_1\omega_3+
\sqrt2\,q_5\omega_2\omega_3.
$$

この成分表現が

$$
I_3=q_A F_A(\omega)
$$

として exact に成立することを確認した。

### G.3 Allowed and forbidden structures

最終基底での整理により、spin-2 と spin-1 の cubic coupling は Clebsch-Gordan 分解

$$
\text{spin-2}\otimes\text{spin-1}\otimes\text{spin-1}
$$

の spin-0 射影として判定できる。 $\omega_i\omega_j$ や $\delta_i\delta_j$ は可換変数の対称積なので

$$
\mathrm{Sym}^2(\text{spin-1})=\text{spin-0}\oplus\text{spin-2}
$$

だけを含み、antisymmetric な spin-1 成分は現れない。したがって spin-2 テンソル $h_{ij}$ と縮約してスカラーを作る独立な cubic coupling は、その spin-2 射影

$$
h_{ij}\omega_i\omega_j,\qquad h_{ij}\delta_i\delta_j
$$

に尽きる。

この事情を final basis の言葉で言い直すと、

- $h_{ij}\omega_i\omega_j$
- $h_{ij}\delta_i\delta_j$

は代表的な許容 cubic term とみなせる。一方、既知の禁制項 $H_{\rm sym}$ と $K_{\rm sym}$ は、この symmetric spin-2 射影に還元されないため、final EFT の独立不変量としては採用しない。言い換えれば、許容されるのは対称 2 階テンソルとしての spin-2 成分を持つ bilinear だけであり、禁制項はその条件を満たさない。

### G.4 Boundary of the claim

本付録が主張するのは、本論で採用する final basis に対する low-order invariant organization である。ここから

- arbitrary non-invariant polynomial の完全生成
- あらゆる非線形結合の一意分類

までは主張しない。本論では、この限定を明示したうえで unified Landau-EFT というラベルを用いる。

### G.5 Computational verification

本付録の不変量基底と $I_3=q_AF_A(\omega)$ の exact check は次のスクリプトで確認されている。

- **不変量生成と exact assertion pass**:
  - `paper03/eft_invariant_generation.py`
