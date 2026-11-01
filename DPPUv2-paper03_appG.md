# 付録G. 不変量展開

対称tracelessな3×3行列では、Cayley–Hamiltonは $h^3-I_2h/2-I_3I/3=0$ となる。hを掛けてtraceを取ると $\mathrm{tr}(h^4)=I_2^2/2$ を得る。指数関数を展開すると、

式(G1)：

$$
\begin{aligned}
\mathrm{tr}(e^{2h})&=3+2I_2+\frac43I_3+\frac13I_2^2+O(h^5),\\
\mathrm{tr}(e^{4h})&=3+8I_2+\frac{32}{3}I_3+\frac{16}{3}I_2^2+O(h^5).
\end{aligned}
$$

第二式の二倍から第一式の二乗を引けば式(I1)となる。hは対称なので $\lVert e^h\omega\rVert^2=\omega^{\mathsf T}e^{2h}\omega$ であり、twist–tensor項も得られる。展開のスカラー係数は式(M3)を通じて $r,L,\eta,V$ に依存する。これらの係数自体を展開すれば、高次のスカラー結合も生じる。

独立な二ベクトルu,vがあれば、回転不変性は $h_{ij}u_iv_j$ を許す。実際にその不変量が現れるかは、作用と変数の意味によって決まる。本論のmixing座標はcoframe回転を表し、scalar作用から消える。式(I1)にmixingがないことは、この作用に固有の結果である。第4節の微分不変量Qは、別に指定した円方向ベクトル族に属する。

[前の節・付録](DPPUv2-paper03_appF.md) / [次の節・付録](DPPUv2-paper03_appH.md) / [第1節](DPPUv2-paper03_sec01.md)
