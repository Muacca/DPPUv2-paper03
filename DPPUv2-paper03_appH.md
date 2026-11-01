# 付録H. 数値と図の再現

図示の基準点は $r=3,L=\kappa=\eta=V=1$ とし、AX・VTでは無効なtorsion係数を零にする。解析式から、

式(H1)：

$$
\begin{aligned}
8\mathcal K&=24\pi^2,\qquad
H_{\omega\omega}=\frac{2\pi^2}{3},\qquad C_\omega=\frac{\pi^2}{3},\\
H_{\eta\eta}&=H_{VV}=\frac{9\pi^2}{2},\qquad
H_{r\eta}=\frac{3\pi^2}{2},\qquad
H_{rV}=\frac{9\pi^2}{2}\quad(\mathrm{MX}),\\
\epsilon_c&=\frac1{\sqrt2}-1\simeq-0.292893218.
\end{aligned}
$$

paper_results パッケージは静的Cartan計算と、制限関数shape_minimum・twist_bulk・scalar_bulkを含む。circle.joint_bulk関数は式(M4)–(M6)を実装し、静的・速度二次・速度一次の寄与を別々に返す。coupled_periodic_family関数は固定した周期求積で式(X4)を評価する。10組の図には、静的・制限セクターの結果に加えて交叉係数と連成系列を含む。data manifestは入力パラメータ、grid、source hash、data hashを記録する。

リポジトリの script ディレクトリから次を実行する。各コマンドは一行で入力し、出力先には新しいディレクトリを指定する。

    python -B -m paper_results render --output-dir ../rendered
    python -B -m paper_results recompute --output-dir ../recomputed-data
    python -B -m paper_results render --data-dir ../recomputed-data --output-dir ../recomputed-figures

描画は保存配列を使い、幾何計算を呼び出さない。notebookも保存データ経路を既定とし、再計算は明示的に選ぶ。添付はdata・LaTeX・scriptディレクトリに収録する。
[再現手順](script/README_ja.md)、[幾何・計算規約](script/docs/CONVENTIONS_ja.md)、
[実装ガイド](script/docs/IMPLEMENTATION_ja.md) に依存関係・API・viewerの使い方を示す。
[データ一覧](data/revised/manifest.json)、[図notebook](script/DPPUv2_Paper03_Figures_v3.ipynb)、
[LaTeX原稿](LaTeX/main.tex) を添付する。

リポジトリは [Muacca/DPPUv2-paper03](https://github.com/Muacca/DPPUv2-paper03)、保存記録は [Zenodo](https://doi.org/10.5281/zenodo.19144481) にある。

[前の節・付録](DPPUv2-paper03_appG.md) / [次の節・付録](DPPUv2-paper03_appI.md) / [第1節](DPPUv2-paper03_sec01.md)
