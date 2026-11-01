# 付録C. 正規化tensor基底と回転軌道

行列単位を $E_{ij}$ とすると、基底は

式(C1)：

$$
\begin{aligned}
T_1&=\frac{\mathrm{diag}(1,1,-2)}{\sqrt6},&
T_2&=\frac{\mathrm{diag}(1,-1,0)}{\sqrt2},\\
T_3&=\frac{E_{01}+E_{10}}{\sqrt2},&
T_4&=\frac{E_{02}+E_{20}}{\sqrt2},&
T_5&=\frac{E_{12}+E_{21}}{\sqrt2}.
\end{aligned}
$$

$\mathrm{tr}(T_AT_B)=\delta_{AB}$ を満たす。01平面の回転による共役でT1は不変、 $(T_2,T_3)$ は二倍角、 $(T_4,T_5)$ は同じ角度で回転する。従って実spin-2表現は軸対称背景で $1+2+2$ に分解する。doubletの場合、hの01 blockの固有値は $\pm\rho/\sqrt2$ である。そのblockで $M_0$ と可換なので、式(T2)直前の固有値表示が得られる。

shape_matrix関数の対角座標sは $M=\mathrm{diag}(a(1+s),a/(1+s),a^{-2})$ を用いる。正規化したT2座標とのexactな関係は $1+s=e^{q_2/\sqrt2}$ である。断面 $q_3=0$ では90度回転が二軸を交換して $q_2\mapsto-q_2$ とする。従って符号の異なる二枝は同じ残余回転軌道に属する。

[前の節・付録](DPPUv2-paper03_appB.md) / [次の節・付録](DPPUv2-paper03_appD.md) / [第1節](DPPUv2-paper03_sec01.md)
