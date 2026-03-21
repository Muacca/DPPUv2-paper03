## Appendix D. Mass Splitting Formulae in the Spin-2 Sector

本付録では、spin-2 quintet の residual SO(2) 分裂、等方点近傍での branch splitting、臨界点近傍の相図情報を整理する。

### D.1 Residual SO(2) decomposition

正側 squash 背景では $SO(3)\to SO(2)$ となり、spin-2 quintet は

| mode | basis | residual charge |
|---|---|---|
| singlet | $T_1$ | $m=0$ |
| doublet | $(T_2,T_3)$ | $m=\pm2$ |
| doublet | $(T_4,T_5)$ | $m=\pm1$ |

に分裂する。

### D.2 Isotropic point

等方点では Schur の補題により

$$
H_{AB}=2\lambda\,\delta_{AB},
\qquad
\mu_q^2=48\pi^2
$$

が成立し、5 成分は完全に縮退する。

### D.3 Local splitting around the isotropic point

軸対称 slice $s=0$ 上では、3 本の代表 branch を

$$
\mu^2_{m=0}\equiv \mu^2_\varepsilon,\qquad
\mu^2_{m=\pm2}\equiv \mu^2_s,\qquad
\mu^2_{m=\pm1}\equiv \mu^2_\delta
$$

として追跡できる。

ここで $\mu_s^2$ は old variable $s$ を用いて追跡した $(T_2,T_3)$ doublet の radial branch を表す。Appendix C.6 で見たように、 $s\neq0$ 相の角度 $\phi$ は同じ $(q_2,q_3)$ 平面内の SO(2) 角変数であり、独立な mass branch label ではない。したがって本付録では、 $m=\pm2$ sector の代表値として radial direction の固有値を $\mu^2_{m=\pm2}\equiv\mu_s^2$ と記す。

等方点近傍の数値スキャンから、局所線形分裂は

$$
\mu^2_{m=0}(\varepsilon)\approx 48\pi^2 - 1.37\times10^3\,\varepsilon + O(\varepsilon^2),
$$

$$
\mu^2_{m=\pm2}(\varepsilon)\approx 48\pi^2 + 7.95\times10^2\,\varepsilon + O(\varepsilon^2),
$$

$$
\mu^2_{m=\pm1}(\varepsilon)\approx 48\pi^2 - 3.0\times10^2\,\varepsilon + O(\varepsilon^2)
$$

で近似される。ここで係数は $|\varepsilon|\lesssim0.1$ の数値データから得た局所傾きであり、closed-form の厳密不変量ではなく、等方点近傍の branch ordering を要約するための実用式である。

この分裂の向きは

- 正の $\varepsilon$ では $m=0$ singlet が最も速く softening する。
- $m=\pm2$ doublet は同じ領域で硬化し、正側の第一不安定化には関与しない。
- $m=\pm1$ doublet はその中間に位置し、臨界点近傍でも正の質量を保つ。

という §4 の記述と整合する。

### D.4 Closed-form Hessians for all spin-2 modes

以下の式は spin-2 quintet の 3 代表モードに対する $H_{ii}^q$ の閉形式であり、 $a=(1+\varepsilon)^{1/3}$ を用いる（物理パラメータ $\eta=L=\kappa=1$, $\alpha=\theta_\text{NY}=0$, Mode.MX）。

**T₁ (m=0 squash singlet)**:

$$
H_{11}^q(a,r) = \frac{-160\pi^2 a^6 r + 352\pi^2 r}{3a^8}, \qquad a=(1+\varepsilon)^{1/3}
$$

**T₂ (m=±2 shear doublet)**:

$$
H_{22}^q(a,r) = \frac{128\pi^2 a^6 r - 64\pi^2 r}{a^2}, \qquad a=(1+\varepsilon)^{1/3}
$$

**T₄ (m=±1 off-diagonal doublet)**:

$$
H_{44}^q(f,r) = \frac{32\pi^2 f^{12} r - 48\pi^2 f^{6} r + 32\pi^2 f^{3} r + 48\pi^2 r}{f^{8}}, \qquad f=(1+\varepsilon)^{2/3}=a^2
$$

等方点 $\varepsilon=0$（ $a=f=1$ ）では Schur の補題から

$$
H_{11}^q = H_{22}^q = H_{44}^q = 64\pi^2 r_0
$$

が成立し（D.2 節と整合）、スクリプトにより数値的にも確認されている。

閉形式から得られる臨界点（ともに $r_0$ 非依存の厳密解）：

$$
H_{22}^q = 0 \;\iff\; 2a^6=1 \;\Rightarrow\; \varepsilon_{c-} = \frac{1}{\sqrt{2}}-1 \approx -0.292893 \quad\text{(exact)}
$$

$$
H_{11}^q = 0 \;\iff\; 5a^6=11 \;\Rightarrow\; \varepsilon_{c+} = \sqrt{\frac{11}{5}}-1 \approx 0.48324 \quad\text{(exact)}
$$

スクリプト参照：`paper03/spin2_diagonal_hessians.py`（H₁₁, H₂₂ 閉形式）、`paper03/spin2_offdiag_hessian.py`（H₄₄ 閉形式; Method 2 rationalization）。

### D.5 Representative spectrum values

代表点での branch 値は次の通りである。

| $\varepsilon$ | $\mu^2_{m=0}$ | $\mu^2_{m=\pm2}$ | $\mu^2_{m=\pm1}$ |
|---|---:|---:|---:|
| $0.00$ | 473.74 | 473.74 | 473.74 |
| $0.10$ | 336.83 | 589.55 | 442.41 |
| $0.20$ | 219.78 | 685.05 | 413.32 |
| $0.40$ | 44.91 | 842.19 | 359.80 |
| $0.482$ | 0.49 | 891.77 | 341.58 |

この表が示すのは、soft branch が $m=0$ のみであること、および $m=\pm2$ branch が負側転移に向けて別個に動くことである。

### D.6 Positive side

正側では $T_1$ singlet が最初に soft となり

$$
\varepsilon_{c+}=\sqrt{\tfrac{11}{5}}-1\approx 0.48324 \quad\text{(exact)},
\qquad
r_0(\varepsilon_{c+})\approx 2.8114
$$

で $\mu_\varepsilon^2\to0$ となる（§4.4 の記述と対応）。 $\varepsilon_{c+}$ は $r_0$ 非依存の厳密解であり（D.4 節参照）、数値スキャン推定値 $0.48324$ と整合する。本節はその数値的裏付けを提供する。 $\varepsilon_{c+}$ は positive- $\varepsilon$ 側での first softening threshold であって、full condensed phase の確立点ではない点に注意（§4.4 参照）。

### D.7 Negative side

負側では $m=\pm2$ sector が $s\neq0$ phase への遷移を支配する（§4.5 参照）。閉形式から得られる厳密臨界点は

$$
\varepsilon_{c-}=\tfrac{1}{\sqrt{2}}-1\approx -0.292893 \quad\text{(exact)}
$$

であり（D.4 節参照）、three-indicator 数値解析の推定値 $-0.29309\pm0.00019$ と整合する。Landau 展開

$$
F(s;\varepsilon)=a(\varepsilon)s^2+b(\varepsilon)s^4+O(s^6)
$$

に対して

$$
b(\varepsilon)\in[2288,7805]>0
$$

が確認されるため、相転移は 2 次である。3 指標の収束値は $-0.29303,\,-0.29335,\,-0.29289$ であり（§4.5 参照）、これを数値的に詳細化したものが本節である。うち第三指標 $-0.29289$ は D.4 節の閉形式厳密解 $\tfrac{1}{\sqrt{2}}-1$ と一致しており、この値が H₂₂^q=0 条件から直接導かれた結果であることを示している。

### D.8 Order parameter fit

order parameter は

$$
s_0\propto|\varepsilon-\varepsilon_c|^\beta,
\qquad
\beta=0.4678,\quad R^2=0.987
$$

でフィットされ、連続相転移の記述と整合する（§4.5 参照）。 $\beta$ は mean-field 値 $1/2$ をわずかに下回る（ $1/2-\beta\approx0.032$ ）が、有限レンジでのフィットと高次項の寄与で容易に動く量であり、本論はこの minisuperspace EFT から universality class を抽出することまでは主張しない。

### D.9 Interpretation

spin-2 sector の最終像は、単一の squash 変数ではなく residual representation ごとの mass splitting として理解されるべきである。本文ではこの representation-theoretic language を優先して用いる。

### D.10 Computational verification

本付録の数値は以下のスクリプトで確認されている。

- **等方点 Schur 縮退（ $\mu_q^2=48\pi^2$ ）**:
  - `paper03/spin2_schur_degeneracy.py`
- **SO(2) 分裂・軟化・相転移（三指標解析、 $\varepsilon_{c\pm}$ 数値推定）**:
  - `paper03/spin2_spectrum_splitting.py`
- **H₁₁^q, H₂₂^q 閉形式導出（D.4a, D.4b; $\varepsilon_{c\pm}$ 厳密解）**:
  - `paper03/spin2_diagonal_hessians.py`
- **H₄₄^q 閉形式導出（D.4c; Method 2 rationalization）**:
  - `paper03/spin2_offdiag_hessian.py`
