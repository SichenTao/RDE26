# RDE26-CMOP：IEEE WCCI/CEC 2026 数値最適化コンペティション Competition 4-CMOPs で第 1 位（[公式結果](https://github.com/P-N-Suganthan/2026-CEC)）

<p align="center">
  <a href="README.md">English</a> |
  <a href="README.zh-CN.md">中文</a> |
  <strong>日本語</strong>
</p>

RDE26-CMOP は、IEEE WCCI 2026 / IEEE CEC 2026「精度と速度を考慮した数値最適化コンペティション」の 1 トラックである制約付き多目的最適化に対する RDE26 の出場アルゴリズムです。このトラック（Competition 4-CMOPs）で第 1 位でした。

得点、必要環境、ほかのトラックは[トップレベルの README](../README.ja.md)を参照してください。

## 問題と設定

- テスト問題：決定空間がスケーラブルな制約付き多目的問題 15 問（SDC1～SDC15）。C++ 版は `benchmark/cpp/cmop_cec2026/` にあり、`benchmark/SDC/` の MATLAB ファイルは参考用です。
- D = 30、個体数 N = 100、問題ごとに独立 15 試行、1 試行あたり 200,000 回の関数評価。
- ソルバーは C++17、実行スクリプトは Python 3（サードパーティのパッケージも MATLAB も不要）。

## 内容

- `benchmark/SDC/`：SDC1～SDC15 と PlatEMO の基底クラス（MATLAB）。`benchmark/cpp/cmop_cec2026/`：テストスイートの C++ 版。
- `code/rde26-cmop/main.py`：ビルド・実行スクリプト。`code/rde26-cmop/rde_cmop_native/`：ソルバーのソース。

## 実行

リポジトリのルートから実行します（またはフルパスで指定）：

```bash
python3 RDE26_CMOP/code/rde26-cmop/main.py    # 既定設定での全体実行
python3 RDE26_CMOP/code/rde26-cmop/main.py --problems 1 --runs 1 --maxFE 5000    # 短い動作確認
```

## `main.py` のオプション

| オプション | 既定値 | 意味 |
| --- | --- | --- |
| `--problems` | `1:15` | 実行する問題。リストまたは範囲 |
| `--runs` | `15` | 問題ごとの試行回数 |
| `--seed-base` | `20260608` | 乱数シードの基準値 |
| `--N` | `100` | 個体数 |
| `--D` | `30` | 決定変数の数 |
| `--maxFE` | `200000` | 1 試行あたりの評価回数 |
| `--save` | `1000` | 1 試行あたりの記録点数 |
| `--pf-root` | (from the benchmark folder) | 参照パレートフロントのフォルダ |
| `--out-dir` | `outputs/rde26-cmop_<time>` | 出力フォルダ |
| `--paper-id` | `RDE26` | 出力ファイル名の接頭辞 |
| `--cxx` | `g++` | C++ コンパイラ |
| `--skip-build` | `off` | コンパイル済みプログラムを再利用 |

## 出力

問題ごとに 1 ファイル `RDE26_CEC26_CMOP_F<k>.txt`：1000 行、試行ごとに 4 つの値（タブ区切り：IGD、制約違反スコア、平均制約違反量、実行可能率。定義できない箇所は `nan`）。生のログは `_raw/` に残ります。
