# RDE26: Rank 1 in all four tracks of the IEEE WCCI/CEC 2026 numerical optimization competition ([official results](https://github.com/P-N-Suganthan/2026-CEC))

<p align="center">
  <strong>English</strong> |
  <a href="README.zh-CN.md">中文</a> |
  <a href="README.ja.md">日本語</a>
</p>

RDE26 is the series of differential evolution solvers that the Tohoku University-centered team entered in the four tracks of the IEEE WCCI 2026 / IEEE CEC 2026 Competition on Numerical Optimization Considering Accuracy and Speed. It ranked first in every track of the organizers' comparison. This repository contains the source code of the four entries, the benchmark material needed to run them, and run scripts, one directory per track.

## Results

| Track | Problem class | Rank | sum / RS | Directory |
| --- | --- | --- | --- | --- |
| Competition 1-BC-SOPs | Bound-constrained single-objective | 1 | 281518 / 78 | [RDE26_SOP](./RDE26_SOP) |
| Competition 2-CSOPs | Constrained single-objective | 1 | 192754 / 55 | [RDE26_CSOP](./RDE26_CSOP) |
| Competition 3-BC-MOPs | Bound-constrained multi-objective | 1 | 87261.5 / 16 | [RDE26_MOP](./RDE26_MOP) |
| Competition 4-CMOPs | Constrained multi-objective | 1 | 188840 / 38 | [RDE26_CMOP](./RDE26_CMOP) |

Source: the organizers' comparison slides ([P-N-Suganthan/2026-CEC](https://github.com/P-N-Suganthan/2026-CEC)); sum is the total score and RS the rank sum, as given there. The BC-SOP row follows version V6 of the slides. This is the only place in the repository where scores are given.

## Test sets and default settings

| Directory | Test set | Dimension | Independent runs | Budget per run | Language |
| --- | --- | --- | --- | --- | --- |
| `RDE26_SOP` | 29 CEC 2017 real-parameter single-objective functions | D = 30 | 25 | 10000 x D = 300,000 evaluations | C++14, Python |
| `RDE26_CSOP` | 28 CEC 2017 constrained real-parameter problems | D = 30 | 25 | 20000 x D = 600,000 evaluations | C++11, Python |
| `RDE26_MOP` | 10 multi-objective problems (MaOP1-MaOP10) | M = 3 objectives, D = 10 in `main.py`, population size 100 | 30 | 100,000 evaluations | C++17, make, Python |
| `RDE26_CMOP` | 15 scalable decision-space constrained multi-objective problems (SDC1-SDC15) | D = 30, population size 100 | 15 | 200,000 evaluations | C++17, Python |

For `RDE26_MOP`, the run script uses D = 10 by default (the folder `benchmark/SEC2018_MaOP_M3_D10` has the same name). The dimension, number of objectives and budget used in a competition are set by the organizers; pass them to `main.py` (`--D`, `--M`, `--maxFE`) to match the call of the competition you run.

## Repository layout

```text
.
├── README.md, README.zh-CN.md, README.ja.md
├── CHANGELOG.md, CITATION.cff, LICENSE, THIRD_PARTY_NOTICES.md
├── RDE26_SOP/       bound-constrained single-objective optimization
│   ├── benchmark/CEC2017_SOP/        cec17_test_func.cpp, input_data/
│   └── code/rde26-sop/               main.py, src/
├── RDE26_CSOP/      constrained single-objective optimization
│   ├── benchmark/CEC2017_CSOP/       cec17_test_COP.cpp, inputData/
│   └── code/rde26-csop/              main.py, RDEx.cpp, cec17_test_COP.cpp, inputData/
├── RDE26_MOP/       bound-constrained multi-objective optimization
│   ├── benchmark/SEC2018_MaOP_M3_D10/   MaOP/, POF/, PROBLEM.m, SOLUTION.m
│   └── code/rde26-mop/               main.py, src/
└── RDE26_CMOP/      constrained multi-objective optimization
    ├── benchmark/SDC/                SDC problems (MATLAB)
    ├── benchmark/cpp/cmop_cec2026/   C++ version of the SDC suite
    └── code/rde26-cmop/              main.py, rde_cmop_native/
```

Each directory has its own README with the options and the output format.

## Requirements and quick start

Python 3 and a C++ compiler (called as `g++` by default; `--cxx` selects another one for SOP, CSOP and CMOP; MOP also needs `make`). No third-party Python package is needed. Every run script writes to `outputs/` next to `main.py` unless `--out-dir` is given.

```bash
# full runs with the default settings
python3 RDE26_SOP/code/rde26-sop/main.py     # 29 functions x 25 runs
python3 RDE26_CSOP/code/rde26-csop/main.py   # 28 problems x 25 runs
python3 RDE26_MOP/code/rde26-mop/main.py     # 10 problems x 30 runs
python3 RDE26_CMOP/code/rde26-cmop/main.py   # 15 problems x 15 runs

# short checks: one or two problems, one run
python3 RDE26_SOP/code/rde26-sop/main.py --functions 1,4 --runs 1
python3 RDE26_CSOP/code/rde26-csop/main.py --functions 1,7 --runs 1
python3 RDE26_MOP/code/rde26-mop/main.py --problems 1 --runs 1 --maxFE 5000
python3 RDE26_CMOP/code/rde26-cmop/main.py --problems 1 --runs 1 --maxFE 5000
```

The short checks were run on macOS (arm64) with Apple clang and Python 3; Linux and Windows have not been tried.

## History of this repository

Version 1.0.0 combines the four stand-alone repositories RDE26-SOP, RDE26-CSOP, RDE26-MOP and RDE26-CMOP, as they were on 2026-07-09, in one repository with one directory per track. The algorithm code is unchanged; see [CHANGELOG.md](CHANGELOG.md). RDE26-SOP was renamed to RDE26 (the old address redirects here). RDE26-CSOP, RDE26-MOP and RDE26-CMOP are archived and point to the directories above.

## Papers

There is no RDE26 paper yet. <!-- TODO: add the RDE26 paper link here and in CITATION.cff when it is available -->The solvers follow the RDE and RDEx papers listed in the [series overview](https://github.com/SichenTao/reconstructed-differential-evolution). The organizers' report on the competition: K. Qiao, X. Ban, P. Chen, K. V. Price, P. N. Suganthan, J. Liang, C. Yue, and G. Wu, "Performance comparison of CEC 2026 competition entries on numerical optimization considering accuracy and speed," IEEE WCCI/CEC 2026 competition slides ([P-N-Suganthan/2026-CEC](https://github.com/P-N-Suganthan/2026-CEC)).

## Team

Sichen Tao, Hanyu Hu, Ruihan Zhao, Qingke Zhang, Yifei Yang, Jian Wang, Masatoshi Kawai, and Hiroyuki Takizawa.

## Related repositories

- [reconstructed-differential-evolution](https://github.com/SichenTao/reconstructed-differential-evolution): overview of the RDE line (2024 RDE, 2025 RDEx, 2026 RDE26).
- [RDE-SOR](https://github.com/SichenTao/RDE-SOR): a separate algorithm by a different team, also entered in the BC-MOP track in 2026; kept in its own repository.

## Citation and license

Use the "Cite this repository" button (metadata in [CITATION.cff](CITATION.cff)). A DOI for each release will be added after the release is archived.

The code is released under the [MIT License](LICENSE). The benchmark programs and data under `benchmark/` and a few files under `code/` come from other authors; see [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
