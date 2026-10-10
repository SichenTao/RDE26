# RDE26-SOP：IEEE WCCI/CEC 2026 数值优化竞赛 Competition 1-BC-SOPs 第 1 名（[官方结果](https://github.com/P-N-Suganthan/2026-CEC)）

<p align="center">
  <a href="README.md">English</a> |
  <strong>中文</strong> |
  <a href="README.ja.md">日本語</a>
</p>

RDE26-SOP 是 RDE26 在有界约束单目标优化赛道的参赛算法，该赛道属于 IEEE WCCI 2026 / IEEE CEC 2026 “兼顾精度与速度的数值优化竞赛”。它在本赛道（Competition 1-BC-SOPs）排名第 1。

分数、环境要求和其他赛道见[顶层 README](../README.zh-CN.md)。

## 问题与设置

- 测试集：29 个 CEC 2017 实参数单目标函数（运行脚本中编号 1–29）。
- D = 30，每个函数独立运行 25 次，每次 10000 x D = 300,000 次函数评价。
- 求解器为 C++14，运行脚本为 Python 3（不需要第三方包）。

## 内容

- `benchmark/CEC2017_SOP/`：`cec17_test_func.cpp` 和 `input_data/`（CEC 2017 评价程序及其数据）。
- `code/rde26-sop/main.py`：编译并运行的脚本。`code/rde26-sop/src/`：求解器源码和它读取的输入数据。

## 运行

在仓库根目录运行（或写完整路径）：

```bash
python3 RDE26_SOP/code/rde26-sop/main.py    # 默认设置的完整运行
python3 RDE26_SOP/code/rde26-sop/main.py --functions 1,4 --runs 1    # 简短检查
```

## `main.py` 的选项

| 选项 | 默认值 | 含义 |
| --- | --- | --- |
| `--functions` | `1:29` | 要运行的问题，列表（`1,4`）或区间（`1:29`） |
| `--runs` | `25` | 每个问题的运行次数 |
| `--seed-base` | `20260601` | 随机种子的基数 |
| `--max-evaluations` | `300000` | 每次运行的评价次数 |
| `--out-dir` | `outputs/rde26-sop_<time>` | 输出目录 |
| `--paper-id` | `RDE26` | 输出文件名的前缀 |
| `--cxx` | `g++` | C++ 编译器 |
| `--skip-build` | `off` | 沿用已编译的程序 |

## 输出

每个函数一个文件 `RDE26_CEC26_SOP_F<k>.txt`：1000 行（每 300 次评价一行），每次运行一列（制表符分隔），记录到当前为止最优解的误差。
