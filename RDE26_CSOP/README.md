# RDE26-CSOP: Rank 1 in Competition 2-CSOPs of the IEEE WCCI/CEC 2026 numerical optimization competition ([official results](https://github.com/P-N-Suganthan/2026-CEC))

<p align="center">
  <strong>English</strong> |
  <a href="README.zh-CN.md">中文</a> |
  <a href="README.ja.md">日本語</a>
</p>

RDE26-CSOP is the RDE26 entry for constrained single-objective optimization, one of the four tracks of the IEEE WCCI 2026 / IEEE CEC 2026 Competition on Numerical Optimization Considering Accuracy and Speed. It ranked first in this track (Competition 2-CSOPs).

Scores, requirements and the other tracks: see the [top-level README](../README.md).

## Problem and settings

- Test set: 28 CEC 2017 constrained real-parameter problems (numbered 1-28 in the run script).
- D = 30, 25 independent runs per problem, 20000 x D = 600,000 function evaluations per run.
- C++11 for the solver, Python 3 for the run script (no third-party package).

## Contents

- `benchmark/CEC2017_CSOP/`: `cec17_test_COP.cpp` and `inputData/` (the CEC 2017 constrained evaluator and its data).
- `code/rde26-csop/`: `main.py`, the solver `RDEx.cpp`, the evaluator `cec17_test_COP.cpp` and `inputData/`.

## Run

Run from the repository root (or give the full path):

```bash
python3 RDE26_CSOP/code/rde26-csop/main.py    # full run with the default settings
python3 RDE26_CSOP/code/rde26-csop/main.py --functions 1,7 --runs 1    # short check
```

## Options of `main.py`

| Option | Default | Meaning |
| --- | --- | --- |
| `--functions` | `1:28` | problems to run, as a list or a range |
| `--runs` | `25` | runs per problem |
| `--seed-base` | `20260611` | base of the random seeds |
| `--out-dir` | `outputs/rde26-csop_<time>` | output folder |
| `--paper-id` | `RDE26` | prefix of the output file names |
| `--cxx` | `g++` | C++ compiler |
| `--skip-build` | `off` | reuse the compiled program |

## Output

One file per problem, `RDE26_CEC26_CSOP_F<k>.txt`: 2000 rows (checkpoints) with two tab-separated values per run (objective value and constraint measure of the best solution at that checkpoint). The raw program output is kept in `_raw/`.
