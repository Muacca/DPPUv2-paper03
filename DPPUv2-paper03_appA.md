# 付録A. Cartan接続とスカラー円方向作用

$de^a=\tfrac12C^a{}_{bc}e^b\wedge e^c$ 、 $`\Gamma^a{}_b=\Gamma^a{}_{bc}e^c`$ とする。計量整合性と指定した幾何学的torsionから、

式(A1)：

$$
\begin{aligned}
\Gamma^{\mathrm{LC}}_{abc}&=\tfrac12(C_{abc}+C_{cba}-C_{bac}),\\
K_{abc}&=-\tfrac12(T_{abc}+T_{bca}-T_{cab}),\qquad
\Gamma^{\mathrm{EC}}=\Gamma^{\mathrm{LC}}+K,\\
R^a{}_{bcd}&=e_c(\Gamma^a{}_{bd})-e_d(\Gamma^a{}_{bc})
+\Gamma^a{}_{ec}\Gamma^e{}_{bd}-\Gamma^a{}_{ed}\Gamma^e{}_{bc}
+\Gamma^a{}_{be}C^e{}_{cd}.
\end{aligned}
$$

代入により $`C^a{}_{bc}+\Gamma^a{}_{cb}-\Gamma^a{}_{bc}=T^a{}_{bc}`$ を確認できる。 $\bar e=(\sigma^0,\sigma^1,\sigma^2,dy)$ 、 $e=A(y)\bar e$ に対する完全な外微分は

式(A2)：

$$
\begin{aligned}
C^a{}_{bc}
&=A^a{}_d(A^{-1})^e{}_b(A^{-1})^f{}_c\bar C^d{}_{ef}\\
&\quad +(\partial_y A^a{}_f)
\left[(A^{-1})^3{}_b(A^{-1})^f{}_c
-(A^{-1})^3{}_c(A^{-1})^f{}_b\right],\\
e_c(\Gamma^a{}_{bd})&=(A^{-1})^3{}_c\partial_y\Gamma^a{}_{bd}.
\end{aligned}
$$

$`A=O\left(\begin{smallmatrix}rM&0\\L\omega^{\mathsf T}&L\end{smallmatrix}\right)`$ が一定なら、これらは一定coframeの代数に帰着する。 $S=MM^{\mathsf T}$ によるscalar縮約は

式(A3)：

$$
R=\frac{8}{r^2}\lbrace[\mathrm{tr}S]^2-2\mathrm{tr}(S^2)\rbrace
-\frac{8L^2}{r^4}\lVert M\omega\rVert^2-\frac{6\eta^2}{r^2}-\frac{2V^2}{3}.
$$

対角 $M=\mathrm{diag}(f_0,f_1,f_2)$ 、 $\prod_i f_i=1$ では空間構造係数は $4f_i^2\epsilon_{ijk}/r$ である。これを代入すると式(A3)のshape項が得られ、式(A2)の行列変換が一般のcoframeとtwistを与える。 $-\pi^2Lr^3/(8\kappa^2)$ を掛けると式(M3)となる。スカラーの円方向族では式(A1)–(A2)から直接、

式(A4)：

$$
\begin{aligned}
R_{\mathrm{EC}}
&=\frac{24-6\eta^2}{r^2}-\frac{2V^2}{3}
-\frac{6r''}{L^2r}-\frac{6(r')^2}{L^2r^2}
+\frac{6r'L'}{L^3r}
+\frac{2}{L}\left(V'+\frac{3r'V}{r}\right),\\
N&=-\frac{6}{Lr^3}\partial_y(\eta r^2).
\end{aligned}
$$

bulkを分離するため体積密度を掛ける。 $R_0=(24-6\eta^2)/r^2-2V^2/3$ と置くと、

式(A5)：

$$
\begin{aligned}
Lr^3R_{\mathrm{EC}}
&=Lr^3R_0+\frac{6r(r')^2}{L}
+\partial_y\left(-\frac{6r^2r'}{L}+2r^3V\right),\\
Lr^3N&=-6\partial_y(\eta r^2).
\end{aligned}
$$

滑らかな周期profileでは境界項の積分は零となる。これにより符号・正規化を含めて式(S4)が証明される。torsion微分項を省略したのではなく、その寄与が表示した全微分に含まれている。

[前の節・付録](DPPUv2-paper03_sec08.md) / [次の節・付録](DPPUv2-paper03_appB.md) / [第1節](DPPUv2-paper03_sec01.md)
