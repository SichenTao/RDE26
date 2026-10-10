# RDE26：IEEE WCCI/CEC 2026 数值优化竞赛四个赛道全部第 1 名（[官方结果](https://github.com/P-N-Suganthan/2026-CEC)）

<p align="center">
  <a href="README.md">English</a> |
  <strong>中文</strong> |
  <a href="README.ja.md">日本語</a>
</p>

RDE26 是以东北大学为中心的团队参加 IEEE WCCI 2026 / IEEE CEC 2026 “兼顾精度与速度的数值优化竞赛”四个赛道的差分进化求解器系列。在组织方的比较中，四个赛道均排名第 1。本仓库包含四个参赛算法的源代码、运行所需的基准材料和运行脚本，每个赛道一个目录。

## 结果

| 赛道 | 问题类别 | 名次 | sum / RS | 目录 |
| --- | --- | --- | --- | --- |
| Competition 1-BC-SOPs | 有界约束单目标 | 1 | 281518 / 78 | [RDE26_SOP](./RDE26_SOP) |
| Competition 2-CSOPs | 约束单目标 | 1 | 192754 / 55 | [RDE26_CSOP](./RDE26_CSOP) |
| Competition 3-BC-MOPs | 有界约束多目标 | 1 | 87261.5 / 16 | [RDE26_MOP](./RDE26_MOP) |
| Competition 4-CMOPs | 约束多目标 | 1 | 188840 / 38 | [RDE26_CMOP](./RDE26_CMOP) |

出处：组织方的比较幻灯片（[P-N-Suganthan/2026-CEC](https://github.com/P-N-Suganthan/2026-CEC)）；sum 为总分，RS 为名次和，均按幻灯片所给。BC-SOP 一行采用幻灯片 V6 版。全仓库只在这里给出分数。

## 测试集与默认设置

| 目录 | 测试集 | 维度 | 独立运行次数 | 每次运行预算 | 语言 |
| --- | --- | --- | --- | --- | --- |
| `RDE26_SOP` | 29 个 CEC 2017 实参数单目标函数 | D = 30 | 25 | 10000 x D = 300,000 次评价 | C++14、Python |
| `RDE26_CSOP` | 28 个 CEC 2017 约束实参数问题 | D = 30 | 25 | 20000 x D = 600,000 次评价 | C++11、Python |
| `RDE26_MOP` | 10 个多目标问题（MaOP1–MaOP10） | 目标数 M = 3，`main.py` 中 D = 10，种群规模 100 | 30 | 100,000 次评价 | C++17、make、Python |
| `RDE26_CMOP` | 15 个可扩展决策空间约束多目标问题（SDC1–SDC15） | D = 30，种群规模 100 | 15 | 200,000 次评价 | C++17、Python |

`RDE26_MOP` 的运行脚本默认 D = 10（与文件夹 `benchmark/SEC2018_MaOP_M3_D10` 同名）。比赛使用的维度、目标数和预算由组织方规定；复现某届比赛时，请按该届通知用 `--D`、`--M`、`--maxFE` 传给 `main.py`。

## 目录结构

```text
.
├── README.md, README.zh-CN.md, README.ja.md
├── CHANGELOG.md, CITATION.cff, LICENSE, THIRD_PARTY_NOTICES.md
├── RDE26_SOP/       有界约束单目标优化
│   ├── benchmark/CEC2017_SOP/        cec17_test_func.cpp, input_data/
│   └── code/rde26-sop/               main.py, src/
├── RDE26_CSOP/      约束单目标优化
│   ├── benchmark/CEC2017_CSOP/       cec17_test_COP.cpp, inputData/
│   └── code/rde26-csop/              main.py, RDEx.cpp, cec17_test_COP.cpp, inputData/
├── RDE26_MOP/       有界约束多目标优化
│   ├── benchmark/SEC2018_MaOP_M3_D10/   MaOP/, POF/, PROBLEM.m, SOLUTION.m
│   └── code/rde26-mop/               main.py, src/
└── RDE26_CMOP/      约束多目标优化
    ├── benchmark/SDC/                SDC 问题（MATLAB）
    ├── benchmark/cpp/cmop_cec2026/   SDC 测试集的 C++ 版本
    └── code/rde26-cmop/              main.py, rde_cmop_native/
```

每个目录有自己的 README，写明选项和输出格式。

## 环境要求与快速开始

需要 Python 3 和 C++ 编译器（默认调用 `g++`；SOP、CSOP、CMOP 可用 `--cxx` 换编译器；MOP 另需 `make`）。不需要第三方 Python 包。运行脚本默认把结果写到 `main.py` 旁的 `outputs/`，可用 `--out-dir` 改。

```bash
# 默认设置的完整运行
python3 RDE26_SOP/code/rde26-sop/main.py     # 29 个函数 x 25 次
python3 RDE26_CSOP/code/rde26-csop/main.py   # 28 个问题 x 25 次
python3 RDE26_MOP/code/rde26-mop/main.py     # 10 个问题 x 30 次
python3 RDE26_CMOP/code/rde26-cmop/main.py   # 15 个问题 x 15 次

# 简短检查：一两个问题，运行 1 次
python3 RDE26_SOP/code/rde26-sop/main.py --functions 1,4 --runs 1
python3 RDE26_CSOP/code/rde26-csop/main.py --functions 1,7 --runs 1
python3 RDE26_MOP/code/rde26-mop/main.py --problems 1 --runs 1 --maxFE 5000
python3 RDE26_CMOP/code/rde26-cmop/main.py --problems 1 --runs 1 --maxFE 5000
```

简短检查在 macOS（arm64）、Apple clang 和 Python 3 下运行过；Linux 和 Windows 没有试过。

## 本仓库的来历

1.0.0 版把 RDE26-SOP、RDE26-CSOP、RDE26-MOP、RDE26-CMOP 四个独立仓库在 2026-07-09 的状态合并到一个仓库，每个赛道一个目录。算法代码没有改动，见 [CHANGELOG.md](CHANGELOG.md)。RDE26-SOP 已改名为 RDE26（旧地址会自动跳转到这里）。RDE26-CSOP、RDE26-MOP、RDE26-CMOP 已归档，并指向上面的各目录。

## 论文

目前还没有 RDE26 的论文。<!-- TODO: 论文可得后，在这里和 CITATION.cff 中补上链接 -->求解器沿用 RDE 与 RDEx 的论文，见[系列总览](https://github.com/SichenTao/reconstructed-differential-evolution)。组织方关于本届竞赛的报告：K. Qiao, X. Ban, P. Chen, K. V. Price, P. N. Suganthan, J. Liang, C. Yue, and G. Wu, “Performance comparison of CEC 2026 competition entries on numerical optimization considering accuracy and speed,” IEEE WCCI/CEC 2026 competition slides（[P-N-Suganthan/2026-CEC](https://github.com/P-N-Suganthan/2026-CEC)）。

## 团队

Sichen Tao, Hanyu Hu, Ruihan Zhao, Qingke Zhang, Yifei Yang, Jian Wang, Masatoshi Kawai, and Hiroyuki Takizawa.

## 相关仓库

- [reconstructed-differential-evolution](https://github.com/SichenTao/reconstructed-differential-evolution)：RDE 系列总览（2024 RDE、2025 RDEx、2026 RDE26）。
- [RDE-SOR](https://github.com/SichenTao/RDE-SOR)：另一支团队的另一种算法，2026 年同样参加 BC-MOP 赛道；保留在独立仓库中。

## 引用与许可证

请使用“Cite this repository”按钮（元数据见 [CITATION.cff](CITATION.cff)）。每个版本归档后会补充对应的 DOI。

代码以 [MIT 许可证](LICENSE)发布。`benchmark/` 下的基准程序与数据以及 `code/` 下的少数文件来自其他作者，见 [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)。
