# RDE26-CSOP：IEEE WCCI/CEC 2026 数值优化竞赛 Competition 2-CSOPs 第 1 名（[官方结果](https://github.com/P-N-Suganthan/2026-CEC)）

<p align="center">
  <a href="README.md">English</a> |
  <strong>中文</strong> |
  <a href="README.ja.md">日本語</a>
</p>

RDE26-CSOP 是 RDE26 在约束单目标优化赛道的参赛算法，该赛道属于 IEEE WCCI 2026 / IEEE CEC 2026 “兼顾精度与速度的数值优化竞赛”。它在本赛道（Competition 2-CSOPs）排名第 1。

分数、环境要求和其他赛道见[顶层 README](../README.zh-CN.md)。

## 问题与设置

- 测试集：28 个 CEC 2017 约束实参数问题（运行脚本中编号 1–28）。
- D = 30，每个问题独立运行 25 次，每次 20000 x D = 600,000 次函数评价。
- 求解器为 C++11，运行脚本为 Python 3（不需要第三方包）。

## 内容

- `benchmark/CEC2017_CSOP/`：`cec17_test_COP.cpp` 和 `inputData/`（CEC 2017 约束问题评价程序及其数据）。
- `code/rde26-csop/`：`main.py`、求解器 `RDEx.cpp`、评价程序 `cec17_test_COP.cpp` 和 `inputData/`。

## 运行

在仓库根目录运行（或写完整路径）：

```bash
python3 RDE26_CSOP/code/rde26-csop/main.py    # 默认设置的完整运行
python3 RDE26_CSOP/code/rde26-csop/main.py --functions 1,7 --runs 1    # 简短检查
```

## `main.py` 的选项

| 选项 | 默认值 | 含义 |
| --- | --- | --- |
| `--functions` | `1:28` | 要运行的问题，列表或区间 |
| `--runs` | `25` | 每个问题的运行次数 |
| `--seed-base` | `20260611` | 随机种子的基数 |
| `--out-dir` | `outputs/rde26-csop_<time>` | 输出目录 |
| `--paper-id` | `RDE26` | 输出文件名的前缀 |
| `--cxx` | `g++` | C++ 编译器 |
| `--skip-build` | `off` | 沿用已编译的程序 |

## 输出

每个问题一个文件 `RDE26_CEC26_CSOP_F<k>.txt`：2000 行（检查点），每次运行两个值（制表符分隔：该检查点处最优解的目标值和约束度量）。程序原始输出保存在 `_raw/`。
