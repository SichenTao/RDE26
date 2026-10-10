# RDE26-SOP：IEEE WCCI/CEC 2026 数値最適化コンペティション Competition 1-BC-SOPs で第 1 位（[公式結果](https://github.com/P-N-Suganthan/2026-CEC)）

<p align="center">
  <a href="README.md">English</a> |
  <a href="README.zh-CN.md">中文</a> |
  <strong>日本語</strong>
</p>

RDE26-SOP は、IEEE WCCI 2026 / IEEE CEC 2026「精度と速度を考慮した数値最適化コンペティション」の 1 トラックである有界制約付き単目的最適化に対する RDE26 の出場アルゴリズムです。このトラック（Competition 1-BC-SOPs）で第 1 位でした。

得点、必要環境、ほかのトラックは[トップレベルの README](../README.ja.md)を参照してください。

## 問題と設定

- テスト問題：CEC 2017 実数値単目的関数 29 問（実行スクリプトでは 1～29 番）。
- D = 30、関数ごとに独立 25 試行、1 試行あたり 10000 x D = 300,000 回の関数評価。
- ソルバーは C++14、実行スクリプトは Python 3（サードパーティのパッケージは不要）。

## 内容

- `benchmark/CEC2017_SOP/`：`cec17_test_func.cpp` と `input_data/`（CEC 2017 の評価プログラムとそのデータ）。
- `code/rde26-sop/main.py`：ビルド・実行スクリプト。`code/rde26-sop/src/`：ソルバーのソースと読み込む入力データ。

## 実行

リポジトリのルートから実行します（またはフルパスで指定）：

```bash
python3 RDE26_SOP/code/rde26-sop/main.py    # 既定設定での全体実行
python3 RDE26_SOP/code/rde26-sop/main.py --functions 1,4 --runs 1    # 短い動作確認
```

## `main.py` のオプション

| オプション | 既定値 | 意味 |
| --- | --- | --- |
| `--functions` | `1:29` | 実行する問題。リスト（`1,4`）または範囲（`1:29`） |
| `--runs` | `25` | 問題ごとの試行回数 |
| `--seed-base` | `20260601` | 乱数シードの基準値 |
| `--max-evaluations` | `300000` | 1 試行あたりの評価回数 |
| `--out-dir` | `outputs/rde26-sop_<time>` | 出力フォルダ |
| `--paper-id` | `RDE26` | 出力ファイル名の接頭辞 |
| `--cxx` | `g++` | C++ コンパイラ |
| `--skip-build` | `off` | コンパイル済みプログラムを再利用 |

## 出力

関数ごとに 1 ファイル `RDE26_CEC26_SOP_F<k>.txt`：1000 行（300 評価ごとに 1 行）、試行ごとに 1 列（タブ区切り）。それまでに得られた最良解の誤差を記録します。
