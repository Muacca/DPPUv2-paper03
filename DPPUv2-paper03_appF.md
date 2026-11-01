# 付録F. スカラー微分と二次の表現block

$x=(r,L,\eta,V)$ として式(S1)を微分すると、

式(F1)：

$$
\nabla U_0=\frac{\pi^2}{\kappa^2}
\begin{pmatrix}
L(3V^2r^2+9\eta^2-36)/12\\
r(V^2r^2+9\eta^2-36)/12\\
3Lr\eta/2\\ Lr^3V/6
\end{pmatrix}.
$$

式(F2)：

$$
H=\frac{\pi^2}{\kappa^2}
\begin{pmatrix}
LV^2r/2 &(3V^2r^2+9\eta^2-36)/12&3L\eta/2&LVr^2/2\\
(3V^2r^2+9\eta^2-36)/12&0&3r\eta/2&Vr^3/6\\
3L\eta/2&3r\eta/2&3Lr/2&0\\
LVr^2/2&Vr^3/6&0&Lr^3/6
\end{pmatrix}.
$$

MX点 $(r,L,\eta,V,\kappa)=(3,1,1,1,1)$ の勾配は $\pi^2(0,-9/2,9/2,9/2)$ 、AX点 $(3,1,1,0,1)$ では $\pi^2(-9/4,-27/4,9/2,0)$ である。これらの行列は非零の勾配と併せて解釈する。特にAXの $(r,\eta)$ 行列式は負であるが、その点は全スカラーポテンシャルの停留鞍点ではない。このような点でのHessianの通常固有値は選んだ座標に依存し、粒子質量ではない。

$M=I,\omega=0$ では、式(M3)のshape一次微分は全てのスカラー座標について零である。これを各スカラーで再微分すると、scalar–shape交差項が全て零と分かる。twist項は $\omega$ について偶かつ二次なので、その一次微分もこの点で全てのscalarとshapeに対して零である。mixingは静的scalar作用に入らない。従って既約な二次構成は直接、

式(F3)：

$$
(r,L,\eta,V)\ \oplus\ (\omega_i,\delta_i)\ \oplus\ (q_A),
$$

となる。第一block内のスカラー混合は残り、第二blockは式(V1)、第三blockは $8\mathcal K I_5$ である。これは静的二次形式の整理であり、各変数の独立な運動項を仮定しない。式(S4)のbulkスカラー速度Hessianは $\mathrm{diag}(-3\pi^2r/(2\kappa^2L),0,0,0)$ であり、座標ごとの役割の違いと整合する。

[前の節・付録](DPPUv2-paper03_appE.md) / [次の節・付録](DPPUv2-paper03_appG.md) / [第1節](DPPUv2-paper03_sec01.md)
