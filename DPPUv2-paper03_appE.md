# 付録E. 非線形twist作用の導出

第4節「有限振幅の円方向依存」のスカラー固定・round形状の族を回転前coframeで表すと、左不変基底での空間metric、lapse、shiftは

式(E1)：

$$
\begin{aligned}
\gamma&=r^2I+L^2\omega\omega^{\mathsf T},&
\mathcal N&=\frac{Lr}{\sqrt{r^2+L^2\omega^2}},\\
\beta^i&=\frac{L^2\omega_i}{r^2+L^2\omega^2},&
\gamma^{-1}&=\frac1{r^2}\left(I-\frac{L^2\omega\omega^{\mathsf T}}{r^2+L^2\omega^2}\right).
\end{aligned}
$$

各円座標で空間shiftは $\gamma$ の対称軸に沿い、 $\gamma$ の空間Lie微分は零となる。従って法線の全体符号を除き $K=\gamma'/(2\mathcal N)$ である。直接の行列traceから、

式(E2)：

$$
\begin{aligned}
\mathrm{tr}\left[(\gamma^{-1}K)^2\right]
-\left[\mathrm{tr}(\gamma^{-1}K)\right]^2
&=\frac{L^2}{2r^4}Q,\\
{}^{(3)}R&=\frac{24}{r^2}-\frac{8L^2}{r^4}\omega^2.
\end{aligned}
$$

Euclidean曲率の分解では第一行の負値と発散項 $-2\partial_y(\omega\cdot\omega')/r^2$ を $^{(3)}R$ に加え、式(V2)を得る。式(A1)–(A2)による計算もmetric分解を使わず同じ結果を与える。この比較で四次項の符号と正規化が定まる。

任意の滑らかなmixing回転に対し、 $e'^3$ に双対な最終frameベクトルの円成分を $E_3^y$ とする。 $r,L,\eta,V$ を固定するとtorsion分解は

式(E3)：

$$
\begin{aligned}
R_{\mathrm{EC}}&=R_{\mathrm{LC}}-\frac{6\eta^2}{r^2}-\frac{2V^2}{3}
+2V\partial_yE_3^y,\\
N&=-\frac{6\eta}{r}\partial_yE_3^y.
\end{aligned}
$$

体積密度は一定で、空間の左不変発散は零である。従って周期profileでは最後の項の積分が零となる。単一の03回転、 $\omega=(w,0,0)$ では $E_3^y=\cos\chi/L+w\sin\chi/r$ となり、境界currentを明示できる。これにより円方向族でのmixingに対するscalar作用の独立性が示される。

さらに $p=\omega'$ 、 $s=\omega\cdot p$ と置くと、微分項の運動量とEuler式は

式(E4)：

$$
\begin{aligned}
\frac{\partial\mathcal L}{\partial p}
&=\frac{C_\omega}{8}(\omega^2p-s\omega),\\
E&=2C_\omega\omega+\frac{C_\omega}{8}
\left[(2p^2+\omega\cdot\omega'')\omega
-2sp-\omega^2\omega''\right].
\end{aligned}
$$

$\omega$ で縮約すると二階微分が消え、式(V5)を得る。これは点ごとの証明であり、極座標 $\rho,n$ が特異になる点でも有効である。振幅で割る必要がなく、スカラー固定・round形状のこの族に属する実の滑らかなtwist profile全体で、twistゼロの停留結果を示す。

[前の節・付録](DPPUv2-paper03_appD.md) / [次の節・付録](DPPUv2-paper03_appF.md) / [第1節](DPPUv2-paper03_sec01.md)
