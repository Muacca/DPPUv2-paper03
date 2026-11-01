# 付録D. shape微分と臨界近傍

$D=M_0^2$ 、 $M=M_0e^{qT}$ とする。traceの循環性から $A(q)=\mathrm{tr}(De^{2qT})$ 、 $B(q)=\mathrm{tr}(De^{2qT}De^{2qT})$ であり、 $q=0$ で

式(D1)：

$$
\begin{aligned}
A_0&=\mathrm{tr}D,& A'_0&=2\mathrm{tr}(DT),&
A''_0&=4\mathrm{tr}(DT^2),\\
B'_0&=4\mathrm{tr}(D^2T),&
B''_0&=8\mathrm{tr}(D^2T^2)+8\mathrm{tr}(DTDT),\\
F''_0&=2B''_0-2(A'_0)^2-2A_0A''_0.
\end{aligned}
$$

式(C1)を代入すると式(T1)を得る。正規化した線形経路 $M=M_0(I+qT)/\det(I+qT)^{1/3}$ では $\log(M_0^{-1}M)=qT-q^2(T^2-I/3)/2+O(q^3)$ となり、曲率は

式(D2)：

$$
\begin{aligned}
\widetilde H_0&=\frac{4\mathcal K(7-a^6)}{3a^8},\\
\widetilde H_2&=\frac{4\mathcal K(4a^6-1)(3a^6-1)}{3a^8},\\
\widetilde H_1&=\frac{2\mathcal K(7a^6+5)}{3a^8}.
\end{aligned}
$$

両経路は一次の接線が同じでも加速度が異なる。連鎖律 $d^2U/dq^2=U_{,\mu\nu}x'^\mu x'^\nu+U_{,\mu}x''^\mu$ がこの差を説明する。shapeが停留するround点では加速度項が消える。

doubletについて $u=a^6=(1+\epsilon)^2$ とする。式(T2)の凸性と $z\ge2$ により、 $u\ge1/2$ では $z=2$ 、それ以外では $z=u^{-1}$ が極小となる。 $z=2\cosh(\sqrt2\rho)$ の展開から、

式(D3)：

$$
F(a,\rho)=F(a,0)
+\frac{4(2a^6-1)}{a^2}\rho^2
+\frac{2(8a^6-1)}{3a^2}\rho^4+O(\rho^6).
$$

臨界点では二次係数が零、四次係数が $2\thinspace 2^{1/3}$ となる。 $1/[2(1+\epsilon)^2]=1+2\sqrt2(\epsilon_c-\epsilon)+O((\epsilon_c-\epsilon)^2)$ なので、exact停留条件から平方根則が得られる。式(T3)の二枝を微分すれば一次微分の一致と式(T4)の跳びを得る。臨界点より下でのポテンシャル差は $F_{\min}-F(a,0)=-(1-2a^6)^2/a^8$ である。

[前の節・付録](DPPUv2-paper03_appC.md) / [次の節・付録](DPPUv2-paper03_appE.md) / [第1節](DPPUv2-paper03_sec01.md)
