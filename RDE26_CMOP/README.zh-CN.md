# RDE26-CMOP：IEEE WCCI/CEC 2026 数值优化竞赛 Competition 4-CMOPs 第 1 名（[官方结果](https://github.com/P-N-Suganthan/2026-CEC)）

<p align="center">
  <a href="README.md">English</a> |
  <strong>中文</strong> |
  <a href="README.ja.md">日本語</a>
</p>

RDE26-CMOP 是 RDE26 在约束多目标优化赛道的参赛算法，该赛道属于 IEEE WCCI 2026 / IEEE CEC 2026 “兼顾精度与速度的数值优化竞赛”。它在本赛道（Competition 4-CMOPs）排名第 1。

分数、环境要求和其他赛道见[顶层 README](../README.zh-CN.md)。

## 问题与设置

- 测试集：15 个可扩展决策空间约束多目标问题 SDC1–SDC15（C++ 版在 `benchmark/cpp/cmop_cec2026/`，`benchmark/SDC/` 中的 MATLAB 文件供参考）。
- D = 30，种群规模 N = 100，每个问题独立运行 15 次，每次 200,000 次函数评价。
- 求解器为 C++17，运行脚本为 Python 3（不需要第三方包，也不需要 MATLAB）。

## 内容

- `benchmark/SDC/`：SDC1–SDC15 和 PlatEMO 基类（MATLAB）。`benchmark/cpp/cmop_cec2026/`：测试集的 C++ 版本。
- `code/rde26-cmop/main.py`：编译并运行的脚本。`code/rde26-cmop/rde_cmop_native/`：求解器源码。

## 运行

在仓库根目录运行（或写完整路径）：

```bash
python3 RDE26_CMOP/code/rde26-cmop/main.py    # 默认设置的完整运行
python3 RDE26_CMOP/code/rde26-cmop/main.py --problems 1 --runs 1 --maxFE 5000    # 简短检查
```

## `main.py` 的选项

| 选项 | 默认值 | 含义 |
| --- | --- | --- |
| `--problems` | `1:15` | 要运行的问题，列表或区间 |
| `--runs` | `15` | 每个问题的运行次数 |
| `--seed-base` | `20260608` | 随机种子的基数 |
| `--N` | `100` | 种群规模 |
| `--D` | `30` | 决策变量个数 |
| `--maxFE` | `200000` | 每次运行的评价次数 |
| `--save` | `1000` | 每次运行记录的点数 |
| `--pf-root` | (from the benchmark folder) | 参考帕累托前沿所在目录 |
| `--out-dir` | `outputs/rde26-cmop_<time>` | 输出目录 |
| `--paper-id` | `RDE26` | 输出文件名的前缀 |
| `--cxx` | `g++` | C++ 编译器 |
| `--skip-build` | `off` | 沿用已编译的程序 |

## 输出

每个问题一个文件 `RDE26_CEC26_CMOP_F<k>.txt`：1000 行，每次运行四个值（制表符分隔：IGD、约束违反得分、平均约束违反量、可行率；无定义处为 `nan`）。原始日志保存在 `_raw/`。
