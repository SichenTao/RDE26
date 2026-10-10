# RDE26-MOP：IEEE WCCI/CEC 2026 数值优化竞赛 Competition 3-BC-MOPs 第 1 名（[官方结果](https://github.com/P-N-Suganthan/2026-CEC)）

<p align="center">
  <a href="README.md">English</a> |
  <strong>中文</strong> |
  <a href="README.ja.md">日本語</a>
</p>

RDE26-MOP 是 RDE26 在有界约束多目标优化赛道的参赛算法，该赛道属于 IEEE WCCI 2026 / IEEE CEC 2026 “兼顾精度与速度的数值优化竞赛”。它在本赛道（Competition 3-BC-MOPs）排名第 1。

分数、环境要求和其他赛道见[顶层 README](../README.zh-CN.md)。

## 问题与设置

- 测试集：10 个多目标问题 MaOP1–MaOP10，参考帕累托前沿在 `benchmark/SEC2018_MaOP_M3_D10/POF`。
- `main.py` 的默认值：目标数 M = 3，决策变量 D = 10，种群规模 N = 100，每个问题独立运行 30 次，每次 100,000 次函数评价。C++ 程序自身的默认 D 是 7，但 `main.py` 总是传入自己的值。
- 比赛使用的维度、目标数和预算由组织方规定。要与某届通知一致，请用选项传入，例如 `--D 7 --M 3 --maxFE 100000`。
- 求解器需要 C++17 和 `make`，运行脚本为 Python 3（不需要第三方包）。

## 内容

- `benchmark/SEC2018_MaOP_M3_D10/`：`MaOP/`（十个问题）、`POF/`（参考前沿）、`PROBLEM.m`、`SOLUTION.m`。
- `code/rde26-mop/main.py`：编译并运行的脚本。`code/rde26-mop/src/`：求解器源码和 `Makefile`。

## 运行

在仓库根目录运行（或写完整路径）：

```bash
python3 RDE26_MOP/code/rde26-mop/main.py    # 默认设置的完整运行
python3 RDE26_MOP/code/rde26-mop/main.py --problems 1 --runs 1 --maxFE 5000    # 简短检查
```

## `main.py` 的选项

| 选项 | 默认值 | 含义 |
| --- | --- | --- |
| `--problems` | `1:10` | 要运行的问题，列表或区间 |
| `--runs` | `30` | 每个问题的运行次数 |
| `--seed-base` | `2025170000` | 随机种子的基数 |
| `--N` | `100` | 种群规模 |
| `--M` | `3` | 目标数 |
| `--D` | `10` | 决策变量个数 |
| `--maxFE` | `100000` | 每次运行的评价次数 |
| `--save` | `501` | 保存的种群数加 1（每次运行 500 个 IGD 值） |
| `--pf-root` | `benchmark/SEC2018_MaOP_M3_D10/POF` | 参考帕累托前沿所在目录 |
| `--out-dir` | `outputs/rde26-mop_<time>` | 输出目录 |
| `--paper-id` | `RDE26` | 输出文件名的前缀 |
| `--skip-build` | `off` | 沿用已编译的程序 |

## 输出

每个问题一个文件 `RDE26_CEC26_MOP_F<k>.txt`：500 行，每次运行一列（制表符分隔），记录所保存种群的 IGD。原始轨迹保存在 `_raw/`。
