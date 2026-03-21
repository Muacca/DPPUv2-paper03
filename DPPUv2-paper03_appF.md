## Appendix F. Cross-Hessian Checks and Quadratic Decoupling

本付録では、physical spin-1 mode $X_k$ が二次レベルで他 sector と decouple することを示す。

### F.1 Setup

physical basis は

$$
X_k=\frac{-2\omega_k+3\delta_k}{\sqrt{13}}
$$

である。各 mixed Hessian は中央差分

$$
\frac{\partial^2 V_{\rm eff}}{\partial p\,\partial X_k}
\approx
\frac{V(p_0+h,+X_k)-V(p_0+h,-X_k)-V(p_0-h,+X_k)+V(p_0-h,-X_k)}{4h^2}
$$

で評価した。

### F.2 spin-1 x spin-2

$$
M_{Ai}^{(qX)}=
\left.\frac{\partial^2V_{\rm eff}}{\partial q_A\partial X_i}\right|_{\rm iso}=0
$$

が全 15 成分で成立した。

| variable | $X_0$ | $X_1$ | $X_2$ |
|---|---:|---:|---:|
| $\varepsilon$ | 0.0000 | 0.0000 | 0.0000 |
| $s$ | 0.0000 | 0.0000 | 0.0000 |
| $q_3$ | 0.0000 | 0.0000 | 0.0000 |
| $q_4$ | 0.0000 | 0.0000 | 0.0000 |
| $q_5$ | 0.0000 | 0.0000 | 0.0000 |

### F.3 spin-1 x radial and torsion

radial と torsion についても

$$
M_i^{(rX)}=0,
\qquad
M_i^{(\eta X)}=0
$$

が成立した。

| sector | $X_0$ | $X_1$ | $X_2$ |
|---|---:|---:|---:|
| $r$ | 0.0000 | 0.0000 | 0.0000 |
| $\eta$ | 0.0000 | 0.0000 | 0.0000 |

### F.4 Step-size stability

step size を変えてもゼロ判定は不変であった。

| $h$ | max abs value over all 21 components |
|---:|---:|
| $10^{-1}$ | 0.0 |
| $10^{-2}$ | 0.0 |
| $10^{-3}$ | 0.0 |
| $5\times10^{-4}$ | 0.0 |
| $10^{-4}$ | 0.0 |

### F.5 Evenness under $X_k\to -X_k$

代表的な非零背景でも

$$
V(p,+X_k)-V(p,-X_k)=0
$$

が数値的に成立した。したがって、少なくとも現実装では $V_{\rm eff}$ は $X_k$ の符号反転に対して偶関数として振る舞い、その結果として cross-Hessian が消えている。

### F.6 Scope boundary

この結果が主張するのは、final basis における **quadratic decoupling** までである。ここから直ちに、非線形 order や任意背景での完全 decoupling までを主張するわけではない。本論ではこの scope を明示したうえで §6 の block-diagonal EFT を述べる。

### F.7 Computational verification

本付録の cross-Hessian 評価は次のスクリプトで確認されている。

- **全 21 成分の差分評価（step size 安定性・偶性確認を含む）**:
  - `paper03/spin1_cross_hessian.py`
- **homogeneous consistency**（antisymmetry・Jacobi pass; AppC.4 と連携）:
  - `paper03/homogeneous_consistency.py`
- **不変量基底の exact check**（AppG と連携）:
  - `paper03/eft_invariant_generation.py`
