# Einstein–Cartan–Nieh–Yan理論におけるS³×S¹のspin-0・spin-1・spin-2変形

Muacca。初版：v1、2026-03-21。前回改訂：v2、2026-06-20。
現行版：**v3、2026-11-01**。

空間spin-0・spin-1・spin-2表現に沿って、静的ポテンシャルと周期的な円方向作用を整理する。
スカラー結合、条件付き半径極小、制限族での正のtwist・shape微分項、
制限された静的shape極小の連続分岐を導く。
スカラー制御変数を固定し、一定の体積保存shapeとtwistを静的に最小化すると、round・twistゼロが選ばれる。

exactな一般円方向作用には非線形微分交叉項とshift結合が現れる。
滑らかなshape–twist連成族により、スカラー固定でもEuclidean作用汎関数が下に有界でないことを示す。
静的極小と円方向作用汎関数の変分上の性質を区別する。
モデルはalpha=0、完全Nieh–Yanの係数は一定である。

- [英語PDF](DPPUv2-paper03.pdf)
- [英語LaTeX](LaTeX/main.tex)
- [日本語本文・第1節](DPPUv2-paper03_sec01.md)
- [計算と図の再現手順](script/README_ja.md)
- [保存データとmanifest](data/revised/manifest.json)
- [図生成notebook](script/DPPUv2_Paper03_Figures_v3.ipynb)
- [改訂履歴](CHANGELOG.md)
- [English README](README.md)

`script/requirements.txt` の依存関係を導入し、`script` から実行する。

    python -B -m paper_results render --output-dir ../rendered
    python -B -m paper_results recompute --output-dir ../recomputed-data
    python -B -m paper_results render --data-dir ../recomputed-data --output-dir ../recomputed-figures

出力先には新しいディレクトリを指定する。最初のコマンドは保存データの描画のみを行い、
再計算は明示的に指定する。notebookも保存データ描画が既定である。
保存データviewerの使い方は計算手順に記載する。
添付コードはpaper_results、図notebook、保存データviewerで構成する。
[規約](script/docs/CONVENTIONS_ja.md) と [実装ガイド](script/docs/IMPLEMENTATION_ja.md) を添付する。

関連する保存記録は [Zenodo](https://doi.org/10.5281/zenodo.19144481) にある。
添付データのパラメータと格子は [manifest](data/revised/manifest.json) に記載する。
以前の保存版には異なる定義による内容が含まれる場合がある。
関連論文：[Paper01](https://github.com/Muacca/DPPUv2-paper01)、
[Paper02](https://github.com/Muacca/DPPUv2-paper02)。
コードとデータは [MIT License](LICENSE.txt) に従う。
連絡先：muacca@dmwp.jp または [Issues](https://github.com/Muacca/DPPUv2-paper03/issues)。

## 日本語本文

- [第1節](DPPUv2-paper03_sec01.md)
- [第2節](DPPUv2-paper03_sec02.md)
- [第3節](DPPUv2-paper03_sec03.md)
- [第4節](DPPUv2-paper03_sec04.md)
- [第5節](DPPUv2-paper03_sec05.md)
- [第6節](DPPUv2-paper03_sec06.md)
- [第7節](DPPUv2-paper03_sec07.md)
- [第8節](DPPUv2-paper03_sec08.md)

- [付録A](DPPUv2-paper03_appA.md)
- [付録B](DPPUv2-paper03_appB.md)
- [付録C](DPPUv2-paper03_appC.md)
- [付録D](DPPUv2-paper03_appD.md)
- [付録E](DPPUv2-paper03_appE.md)
- [付録F](DPPUv2-paper03_appF.md)
- [付録G](DPPUv2-paper03_appG.md)
- [付録H](DPPUv2-paper03_appH.md)
- [付録I](DPPUv2-paper03_appI.md)
- [付録J：一般円方向作用](DPPUv2-paper03_appJ.md)
