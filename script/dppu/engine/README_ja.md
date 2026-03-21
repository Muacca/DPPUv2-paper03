# Engine Layer

⇒ [English](README.md)

計算パイプラインの基盤を提供するモジュール群。

## 概要

15ステップの計算パイプライン、ロギング、チェックポイント機能を提供。

## モジュール

### pipeline.py

計算パイプラインの基底クラス。

**BaseFrameEngine:**

全トポロジーエンジンの抽象基底クラス。15ステップの計算を順次実行:

1. **E4.1**: パラメータ設定
2. **E4.2**: 計量とフレーム定義
3. **E4.3**: Christoffel記号（一般Koszul公式）
4. **E4.4**: トーションテンソル構成
5. **E4.5**: Contortionテンソル
6. **E4.6**: Einstein-Cartan接続
7. **E4.7**: Riemannテンソル（厳密検証付き）
8. **E4.8**: Ricciスカラー
9. **E4.9**: トーションスカラー
10. **E4.10**: Nieh-Yan密度（全バリアント）
11. **E4.11**: ラグランジアン構成
12. **E4.12**: 角度積分
13. **E4.13**: 有効ポテンシャル抽出
14. **E4.14**: 安定性解析
15. **E4.15**: サマリー

**抽象メソッド（サブクラスで実装）:**

```python
@abstractmethod
def _define_structure_constants(self):
    """構造定数 C^a_{bc} を定義"""

@abstractmethod
def _define_metric(self):
    """フレーム計量 η_{ab} を定義"""

@abstractmethod
def _get_volume_factor(self):
    """体積因子を返す"""
```

### logger.py

計算進捗のロギング。

**ComputationLogger:**

ファイルおよびコンソールへのログ出力:

```python
from dppu.engine import ComputationLogger

logger = ComputationLogger('dppu_run.log')
logger.step(1, "E4.1", "Parameter setup")
logger.info("Setting r = 1.0")
logger.success("Step completed")
```

**NullLogger:**

ログ出力を無効化（テスト用）:

```python
from dppu.engine import NullLogger

logger = NullLogger()  # 何も出力しない
```

### checkpoint.py

計算状態の保存・復元。

**CheckpointManager:**

各ステップの状態をpickleで保存:

```python
from dppu.engine import CheckpointManager

ckpt = CheckpointManager('checkpoints/', enabled=True)

# 保存
ckpt.save('E4.7', {'R': riemann_tensor})

# 復元
data = ckpt.load('E4.7')
```

**チェックポイントの有効化:**

```bash
python scripts/run_s3s1.py --mode MX --ny-variant FULL \
    --checkpoint-dir ./checkpoints
```

## 使用例

```python
from dppu.topology import S3S1Engine
from dppu.torsion import Mode, NyVariant
from dppu.engine import ComputationLogger, CheckpointManager

# ロガーとチェックポイントを設定
logger = ComputationLogger('run.log')
ckpt = CheckpointManager('ckpt/', enabled=True)

# エンジン実行
engine = S3S1Engine(Mode.MX, NyVariant.FULL, logger, ckpt)
engine.run()
```

## パイプラインのカスタマイズ

特定ステップのみ実行:

```python
engine = S3S1Engine(Mode.MX, NyVariant.FULL)

# ステップ1-7のみ実行（曲率計算まで）
engine.run(stop_after='E4.7')

# ステップ7から再開
engine.run(resume_from='E4.7')
```

## 依存関係

- SymPy (シンボリック計算)
- pickle (チェックポイント)

## 関連モジュール

- [topology](../topology/README_ja.md): トポロジー固有エンジン
