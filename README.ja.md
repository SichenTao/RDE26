# RDE26：IEEE WCCI/CEC 2026 数値最適化コンペティション 4 トラックすべてで第 1 位（[公式結果](https://github.com/P-N-Suganthan/2026-CEC)）

<p align="center">
  <a href="README.md">English</a> |
  <a href="README.zh-CN.md">中文</a> |
  <strong>日本語</strong>
</p>

RDE26 は、東北大学を中心とするチームが IEEE WCCI 2026 / IEEE CEC 2026「精度と速度を考慮した数値最適化コンペティション」の 4 トラックに出場させた差分進化ソルバーのシリーズです。主催者の比較では 4 トラックすべてで第 1 位でした。このリポジトリには、4 つの出場アルゴリズムのソースコード、実行に必要なベンチマーク資料、実行スクリプトが、トラックごとのディレクトリに入っています。

## 結果

| トラック | 問題の種類 | 順位 | sum / RS | ディレクトリ |
| --- | --- | --- | --- | --- |
| Competition 1-BC-SOPs | 有界制約付き単目的 | 1 | 281518 / 78 | [RDE26_SOP](./RDE26_SOP) |
| Competition 2-CSOPs | 制約付き単目的 | 1 | 192754 / 55 | [RDE26_CSOP](./RDE26_CSOP) |
| Competition 3-BC-MOPs | 有界制約付き多目的 | 1 | 87261.5 / 16 | [RDE26_MOP](./RDE26_MOP) |
| Competition 4-CMOPs | 制約付き多目的 | 1 | 188840 / 38 | [RDE26_CMOP](./RDE26_CMOP) |

出典：主催者の比較スライド（[P-N-Suganthan/2026-CEC](https://github.com/P-N-Suganthan/2026-CEC)）。sum は総得点、RS は順位和で、スライドの記載どおりです。BC-SOP の行はスライドの V6 版に従っています。得点を記載するのはこの場所だけです。

## テスト問題とデフォルト設定

| ディレクトリ | テスト問題 | 次元 | 独立試行回数 | 1 試行あたりの予算 | 言語 |
| --- | --- | --- | --- | --- | --- |
| `RDE26_SOP` | CEC 2017 実数値単目的関数 29 問 | D = 30 | 25 | 10000 x D = 300,000 評価 | C++14、Python |
| `RDE26_CSOP` | CEC 2017 制約付き実数値問題 28 問 | D = 30 | 25 | 20000 x D = 600,000 評価 | C++11、Python |
| `RDE26_MOP` | 多目的問題 10 問（MaOP1–MaOP10） | 目的数 M = 3、`main.py` では D = 10、個体数 100 | 30 | 100,000 評価 | C++17、make、Python |
| `RDE26_CMOP` | 決定空間がスケーラブルな制約付き多目的問題 15 問（SDC1–SDC15） | D = 30、個体数 100 | 15 | 200,000 評価 | C++17、Python |

`RDE26_MOP` の実行スクリプトは既定で D = 10 です（フォルダ名 `benchmark/SEC2018_MaOP_M3_D10` と同じ）。コンペティションで使う次元・目的数・予算は主催者が定めます。各回の要項に合わせて `main.py` に `--D`、`--M`、`--maxFE` で指定してください。

## ディレクトリ構成

```text
.
├── README.md, README.zh-CN.md, README.ja.md
├── CHANGELOG.md, CITATION.cff, LICENSE, THIRD_PARTY_NOTICES.md
├── RDE26_SOP/       有界制約付き単目的最適化
│   ├── benchmark/CEC2017_SOP/        cec17_test_func.cpp, input_data/
│   └── code/rde26-sop/               main.py, src/
├── RDE26_CSOP/      制約付き単目的最適化
│   ├── benchmark/CEC2017_CSOP/       cec17_test_COP.cpp, inputData/
│   └── code/rde26-csop/              main.py, RDEx.cpp, cec17_test_COP.cpp, inputData/
├── RDE26_MOP/       有界制約付き多目的最適化
│   ├── benchmark/SEC2018_MaOP_M3_D10/   MaOP/, POF/, PROBLEM.m, SOLUTION.m
│   └── code/rde26-mop/               main.py, src/
└── RDE26_CMOP/      制約付き多目的最適化
    ├── benchmark/SDC/                SDC 問題（MATLAB）
    ├── benchmark/cpp/cmop_cec2026/   SDC テストスイートの C++ 版
    └── code/rde26-cmop/              main.py, rde_cmop_native/
```

各ディレクトリに、オプションと出力形式を説明した README があります。

## 必要環境とクイックスタート

Python 3 と C++ コンパイラが必要です（既定は `g++`。SOP・CSOP・CMOP は `--cxx` で変更できます。MOP は `make` も必要です）。サードパーティの Python パッケージは不要です。実行スクリプトは、`--out-dir` を指定しない限り `main.py` の隣の `outputs/` に結果を書き出します。

```bash
# 既定設定での全体実行
python3 RDE26_SOP/code/rde26-sop/main.py     # 29 関数 x 25 試行
python3 RDE26_CSOP/code/rde26-csop/main.py   # 28 問題 x 25 試行
python3 RDE26_MOP/code/rde26-mop/main.py     # 10 問題 x 30 試行
python3 RDE26_CMOP/code/rde26-cmop/main.py   # 15 問題 x 15 試行

# 短い動作確認：1～2 問題、1 試行
python3 RDE26_SOP/code/rde26-sop/main.py --functions 1,4 --runs 1
python3 RDE26_CSOP/code/rde26-csop/main.py --functions 1,7 --runs 1
python3 RDE26_MOP/code/rde26-mop/main.py --problems 1 --runs 1 --maxFE 5000
python3 RDE26_CMOP/code/rde26-cmop/main.py --problems 1 --runs 1 --maxFE 5000
```

短い動作確認は macOS（arm64）、Apple clang、Python 3 で実行しました。Linux と Windows では試していません。

## このリポジトリの経緯

バージョン 1.0.0 は、RDE26-SOP、RDE26-CSOP、RDE26-MOP、RDE26-CMOP の 4 つの独立したリポジトリの 2026-07-09 時点の内容を、トラックごとのディレクトリにして 1 つにまとめたものです。アルゴリズムのコードは変更していません。[CHANGELOG.md](CHANGELOG.md) を参照してください。RDE26-SOP は RDE26 に改名しました（旧アドレスはここへ転送されます）。RDE26-CSOP、RDE26-MOP、RDE26-CMOP はアーカイブされ、上記のディレクトリを案内しています。

## 論文

RDE26 の論文はまだありません。<!-- TODO: 論文が公開されたら、ここと CITATION.cff にリンクを追加する -->ソルバーは RDE と RDEx の論文に基づいています（[シリーズ概要](https://github.com/SichenTao/reconstructed-differential-evolution)を参照）。主催者によるコンペティションの報告：K. Qiao, X. Ban, P. Chen, K. V. Price, P. N. Suganthan, J. Liang, C. Yue, and G. Wu, "Performance comparison of CEC 2026 competition entries on numerical optimization considering accuracy and speed," IEEE WCCI/CEC 2026 competition slides（[P-N-Suganthan/2026-CEC](https://github.com/P-N-Suganthan/2026-CEC)）。

## チーム

Sichen Tao, Hanyu Hu, Ruihan Zhao, Qingke Zhang, Yifei Yang, Jian Wang, Masatoshi Kawai, and Hiroyuki Takizawa.

## 関連リポジトリ

- [reconstructed-differential-evolution](https://github.com/SichenTao/reconstructed-differential-evolution)：RDE シリーズの概要（2024 RDE、2025 RDEx、2026 RDE26）。
- [RDE-SOR](https://github.com/SichenTao/RDE-SOR)：別のチームによる別のアルゴリズムで、2026 年に同じく BC-MOP トラックに出場しました。独立したリポジトリで管理しています。

## 引用とライセンス

「Cite this repository」ボタンを使ってください（メタデータは [CITATION.cff](CITATION.cff)）。各リリースをアーカイブした後に DOI を追記します。

コードは [MIT ライセンス](LICENSE)で公開しています。`benchmark/` 以下のベンチマークプログラムとデータ、および `code/` 以下の一部のファイルは他の著者のものです。[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) を参照してください。
