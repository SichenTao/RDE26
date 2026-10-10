# RDE26-MOP：IEEE WCCI/CEC 2026 数値最適化コンペティション Competition 3-BC-MOPs で第 1 位（[公式結果](https://github.com/P-N-Suganthan/2026-CEC)）

<p align="center">
  <a href="README.md">English</a> |
  <a href="README.zh-CN.md">中文</a> |
  <strong>日本語</strong>
</p>

RDE26-MOP は、IEEE WCCI 2026 / IEEE CEC 2026「精度と速度を考慮した数値最適化コンペティション」の 1 トラックである有界制約付き多目的最適化に対する RDE26 の出場アルゴリズムです。このトラック（Competition 3-BC-MOPs）で第 1 位でした。

得点、必要環境、ほかのトラックは[トップレベルの README](../README.ja.md)を参照してください。

## 問題と設定

- テスト問題：多目的問題 10 問（MaOP1～MaOP10）。参照パレートフロントは `benchmark/SEC2018_MaOP_M3_D10/POF` にあります。
- `main.py` の既定値：目的数 M = 3、決定変数 D = 10、個体数 N = 100、問題ごとに独立 30 試行、1 試行あたり 100,000 回の関数評価。C++ プログラム自体の既定は D = 7 ですが、`main.py` は常に自分の値を渡します。
- コンペティションで使う次元・目的数・予算は主催者が定めます。各回の要項に合わせるには、`--D 7 --M 3 --maxFE 100000` のようにオプションで指定してください。
- ソルバーには C++17 と `make` が必要で、実行スクリプトは Python 3 です（サードパーティのパッケージは不要）。

## 内容

- `benchmark/SEC2018_MaOP_M3_D10/`：`MaOP/`（10 問）、`POF/`（参照フロント）、`PROBLEM.m`、`SOLUTION.m`。
- `code/rde26-mop/main.py`：ビルド・実行スクリプト。`code/rde26-mop/src/`：ソルバーのソースと `Makefile`。

## 実行

リポジトリのルートから実行します（またはフルパスで指定）：

```bash
python3 RDE26_MOP/code/rde26-mop/main.py    # 既定設定での全体実行
python3 RDE26_MOP/code/rde26-mop/main.py --problems 1 --runs 1 --maxFE 5000    # 短い動作確認
```

## `main.py` のオプション

| オプション | 既定値 | 意味 |
| --- | --- | --- |
| `--problems` | `1:10` | 実行する問題。リストまたは範囲 |
| `--runs` | `30` | 問題ごとの試行回数 |
| `--seed-base` | `2025170000` | 乱数シードの基準値 |
| `--N` | `100` | 個体数 |
| `--M` | `3` | 目的数 |
| `--D` | `10` | 決定変数の数 |
| `--maxFE` | `100000` | 1 試行あたりの評価回数 |
| `--save` | `501` | 保存する個体群の数 + 1（1 試行あたり IGD 500 個） |
| `--pf-root` | `benchmark/SEC2018_MaOP_M3_D10/POF` | 参照パレートフロントのフォルダ |
| `--out-dir` | `outputs/rde26-mop_<time>` | 出力フォルダ |
| `--paper-id` | `RDE26` | 出力ファイル名の接頭辞 |
| `--skip-build` | `off` | コンパイル済みプログラムを再利用 |

## 出力

問題ごとに 1 ファイル `RDE26_CEC26_MOP_F<k>.txt`：500 行、試行ごとに 1 列（タブ区切り）。保存した個体群の IGD を記録します。生のトレースは `_raw/` に残ります。
