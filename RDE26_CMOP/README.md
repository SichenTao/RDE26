# RDE26-CMOP: Rank 1 in Competition 4-CMOPs of the IEEE WCCI/CEC 2026 numerical optimization competition ([official results](https://github.com/P-N-Suganthan/2026-CEC))

<p align="center">
  <strong>English</strong> |
  <a href="README.zh-CN.md">中文</a> |
  <a href="README.ja.md">日本語</a>
</p>

RDE26-CMOP is the RDE26 entry for constrained multi-objective optimization, one of the four tracks of the IEEE WCCI 2026 / IEEE CEC 2026 Competition on Numerical Optimization Considering Accuracy and Speed. It ranked first in this track (Competition 4-CMOPs).

Scores, requirements and the other tracks: see the [top-level README](../README.md).

## Problem and settings

- Test set: 15 scalable decision-space constrained multi-objective problems, SDC1-SDC15 (a C++ version is in `benchmark/cpp/cmop_cec2026/`; the MATLAB files in `benchmark/SDC/` are for reference).
- D = 30, population size N = 100, 15 independent runs per problem, 200,000 function evaluations per run.
- C++17 for the solver, Python 3 for the run script (no third-party package, no MATLAB needed).

## Contents

- `benchmark/SDC/`: SDC1-SDC15 and the PlatEMO base classes (MATLAB). `benchmark/cpp/cmop_cec2026/`: the C++ version of the suite.
- `code/rde26-cmop/main.py`: build and run script. `code/rde26-cmop/rde_cmop_native/`: solver source.

## Run

Run from the repository root (or give the full path):

```bash
python3 RDE26_CMOP/code/rde26-cmop/main.py    # full run with the default settings
python3 RDE26_CMOP/code/rde26-cmop/main.py --problems 1 --runs 1 --maxFE 5000    # short check
```

## Options of `main.py`

| Option | Default | Meaning |
| --- | --- | --- |
| `--problems` | `1:15` | problems to run, as a list or a range |
| `--runs` | `15` | runs per problem |
| `--seed-base` | `20260608` | base of the random seeds |
| `--N` | `100` | population size |
| `--D` | `30` | number of decision variables |
| `--maxFE` | `200000` | budget per run (evaluations) |
| `--save` | `1000` | number of recorded points per run |
| `--pf-root` | (from the benchmark folder) | folder with the reference Pareto fronts |
| `--out-dir` | `outputs/rde26-cmop_<time>` | output folder |
| `--paper-id` | `RDE26` | prefix of the output file names |
| `--cxx` | `g++` | C++ compiler |
| `--skip-build` | `off` | reuse the compiled program |

## Output

One file per problem, `RDE26_CEC26_CMOP_F<k>.txt`: 1000 rows with four tab-separated values per run (IGD, score of constraint violation, mean constraint violation, feasible rate; `nan` where a value is not defined). The raw logs are kept in `_raw/`.
