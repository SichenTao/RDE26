# RDE26-SOP: Rank 1 in Competition 1-BC-SOPs of the IEEE WCCI/CEC 2026 numerical optimization competition ([official results](https://github.com/P-N-Suganthan/2026-CEC))

<p align="center">
  <strong>English</strong> |
  <a href="README.zh-CN.md">中文</a> |
  <a href="README.ja.md">日本語</a>
</p>

RDE26-SOP is the RDE26 entry for bound-constrained single-objective optimization, one of the four tracks of the IEEE WCCI 2026 / IEEE CEC 2026 Competition on Numerical Optimization Considering Accuracy and Speed. It ranked first in this track (Competition 1-BC-SOPs).

Scores, requirements and the other tracks: see the [top-level README](../README.md).

## Problem and settings

- Test set: 29 CEC 2017 real-parameter single-objective functions (numbered 1-29 in the run script).
- D = 30, 25 independent runs per function, 10000 x D = 300,000 function evaluations per run.
- C++14 for the solver, Python 3 for the run script (no third-party package).

## Contents

- `benchmark/CEC2017_SOP/`: `cec17_test_func.cpp` and `input_data/` (the CEC 2017 evaluator and its data).
- `code/rde26-sop/main.py`: build and run script. `code/rde26-sop/src/`: solver source and the input data it reads.

## Run

Run from the repository root (or give the full path):

```bash
python3 RDE26_SOP/code/rde26-sop/main.py    # full run with the default settings
python3 RDE26_SOP/code/rde26-sop/main.py --functions 1,4 --runs 1    # short check
```

## Options of `main.py`

| Option | Default | Meaning |
| --- | --- | --- |
| `--functions` | `1:29` | problems to run, as a list (`1,4`) or a range (`1:29`) |
| `--runs` | `25` | runs per problem |
| `--seed-base` | `20260601` | base of the random seeds |
| `--max-evaluations` | `300000` | budget per run (evaluations) |
| `--out-dir` | `outputs/rde26-sop_<time>` | output folder |
| `--paper-id` | `RDE26` | prefix of the output file names |
| `--cxx` | `g++` | C++ compiler |
| `--skip-build` | `off` | reuse the compiled program |

## Output

One file per function, `RDE26_CEC26_SOP_F<k>.txt`: 1000 rows (one every 300 evaluations) with one tab-separated column per run, holding the error of the best solution found so far.
