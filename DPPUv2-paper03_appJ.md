# 付録J 一般円方向作用と速度構造

## metricによる導出と境界current

誘導計量・lapse・shiftは式(M4)を用いる。空間的に一様なshiftは式(M5)を満たし、外部曲率行列は $K=(\gamma'-\ell)/(2\mathcal N)$ である。Euclidean曲率分解は、内在三次元曲率に $[\mathrm{tr}(\gamma^{-1}K)]^2-\mathrm{tr}[(\gamma^{-1}K)^2]$ と発散を加える。このunimodular代数では不変ベクトルの空間発散は零である。残る法線currentは $-2\sqrt{\det\gamma}\thinspace \mathrm{tr}(\gamma^{-1}K)$ であり、 $\mathrm{tr}Y=0$ と $\mathcal N\sqrt{\det\gamma}=J$ により、

式(J1)：

$$
\begin{aligned}
R_{\mathrm{LC}}&=R_{\mathrm{stat,LC}}
-\frac{\mathrm{tr}(X^2)-(\mathrm{tr}X)^2-2\mathrm{tr}(XY)}
{4\mathcal N^2}+\frac{B_R'}{J},\\
B_R&=-\frac{J}{\mathcal N^2}\mathrm{tr}X.
\end{aligned}
$$

$R_{\mathrm{stat,LC}}$ は式(A3)で $\eta=V=0$ としたものである。 $X-Y$ の二次形式を展開すると、 $Y^2$ は速度ゼロの寄与であり、内在三次元曲率と合わせて $R_{\mathrm{stat,LC}}$ になる。 $-\pi^2J/(8\kappa^2)$ を掛け、周期境界積分を除けば式(M6)を得る。非等方shapeでKillingでないshiftを保持した導出である。lapse・shiftの幾何は [6] を参照し、本Euclidean規約の符号はCartan曲率式(A1)からも定まる。

## torsion微分と周期回転

滑らかなcoframe回転O(y)に対し、回転後の双対frameベクトルの円方向成分を $E=E_3^y$ とする。指定したaxial・trace torsionを計量整合接続へ代入すると、

式(J2)：

$$
\begin{aligned}
R_{\mathrm{EC}}&=R_{\mathrm{LC}}-\frac{6\eta^2}{r^2}
-\frac{2V^2}{3}+\frac2J\partial_y(JVE),\\
N_{\mathrm{NY}}&=-\frac6J\partial_y\left(J\frac{\eta}{r}E\right).
\end{aligned}
$$

torsion normの二次項から二つの代数項、trace-torsionの発散から第一行の最後の項が生じ、axial NY currentから第二行を得る。非回転coframeでは $E=1/L$ なので、twistゼロで式(A4)–(A5)、スカラー固定・round形状で式(E3)へ還元する。体積密度を掛けると両currentは全微分であり、周期積分は零になる。従って $\eta',V'$ とOの微分を含めたうえで、独立なbulk速度項が生じないと分かる。一般のNY密度が各点で零である必要はない。

## 交叉項と速度rank

一定round・twistゼロの点の周りで $M=e^h$ と式(M4)–(M6)を展開すると、式(T6)のround近傍での二次形式と、代表的な三次交叉項(X1)を得る。対角族では $\gamma=\mathrm{diag}(r^2e^{2s}+L^2w^2,r^2)$ であり、直接微分すると式(X2)・(X3)となる。式(V6)の回転shapeでは $\gamma'$ と $\ell$ が同じ01回転生成子とのcommutatorとなり、 $-2\mathrm{tr}(XY)$ の評価から表示した係数が得られる。一般作用はdoubletとvectorの全成分を保持し、単一の静的Hessianから結合を推定してはいない。

速度二次形式では $\gamma^{-1}\gamma'$ が対称行列と相似であり、tracefree部分に正の5方向、traceに負の1方向がある。半径と5つのshape座標で任意の対称な $\delta\gamma$ を生成できるため、 $(r,L,\eta,V,q_1,\ldots,q_5,\omega_0,\omega_1,\omega_2)$ のHessianはrank 6、正・負・零の方向数は $(5,1,6)$ となる。非零twistでの零方向はこれらの座標の混合である。4つのスカラー制御変数を固定した8次元shape–twist Hessianでは、twistゼロで $(5,0,3)$ 、非零twistで $(5,1,2)$ となる。後者ではtwist変分が誘導計量の不足するtrace方向を与える。これは縮約前のEuclidean速度形式であり、lapse・運動量拘束を課した後の物理的偏極数ではない。速度一次項はこのHessianを変えずにEuler式へ寄与する。

[前の節・付録](DPPUv2-paper03_appI.md) / [第1節](DPPUv2-paper03_sec01.md)
