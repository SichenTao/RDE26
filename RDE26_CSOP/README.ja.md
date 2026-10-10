# RDE26-CSOP：IEEE WCCI/CEC 2026 数値最適化コンペティション Competition 2-CSOPs で第 1 位（[公式結果](https://github.com/P-N-Suganthan/2026-CEC)）

<p align="center">
  <a href="README.md">English</a> |
  <a href="README.zh-CN.md">中文</a> |
  <strong>日本語</strong>
</p>

RDE26-CSOP は、IEEE WCCI 2026 / IEEE CEC 2026「精度と速度を考慮した数値最適化コンペティション」の 1 トラックである制約付き単目的最適化に対する RDE26 の出場アルゴリズムです。このトラック（Competition 2-CSOPs）で第 1 位でした。

得点、必要環境、ほかのトラックは[トップレベルの README](../README.ja.md)を参照してください。

## 問題と設定

- テスト問題：CEC 2017 制約付き実数値問題 28 問（実行スクリプトでは 1～28 番）。
- D = 30、問題ごとに独立 25 試行、1 試行あたり 20000 x D = 600,000 回の関数評価。
- ソルバーは C++11、実行スクリプトは Python 3（サードパーティのパッケージは不要）。

## 内容

- `benchmark/CEC2017_CSOP/`：`cec17_test_COP.cpp` と `inputData/`（CEC 2017 制約付き問題の評価プログラムとデータ）。
- `code/rde26-csop/`：`main.py`、ソルバー `RDEx.cpp`、評価プログラム `cec17_test_COP.cpp`、`inputData/`。

## 実行

リポジトリのルートから実行します（またはフルパスで指定）：

```bash
python3 RDE26_CSOP/code/rde26-csop/main.py    # 既定設定での全体実行
python3 RDE26_CSOP/code/rde26-csop/main.py --functions 1,7 --runs 1    # 短い動作確認
```

## `main.py` のオプション

| オプション | 既定値 | 意味 |
| --- | --- | --- |
| `--functions` | `1:28` | 実行する問題。リストまたは範囲 |
| `--runs` | `25` | 問題ごとの試行回数 |
| `--seed-base` | `20260611` | 乱数シードの基準値 |
| `--out-dir` | `outputs/rde26-csop_<time>` | 出力フォルダ |
| `--paper-id` | `RDE26` | 出力ファイル名の接頭辞 |
| `--cxx` | `g++` | C++ コンパイラ |
| `--skip-build` | `off` | コンパイル済みプログラムを再利用 |

## 出力

問題ごとに 1 ファイル `RDE26_CEC26_CSOP_F<k>.txt`：2000 行（チェックポイント）、試行ごとに 2 つの値（タブ区切り：そのチェックポイントでの最良解の目的関数値と制約の指標）。プログラムの生出力は `_raw/` に残ります。
