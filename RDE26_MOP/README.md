# RDE26-MOP: Rank 1 in Competition 3-BC-MOPs of the IEEE WCCI/CEC 2026 numerical optimization competition ([official results](https://github.com/P-N-Suganthan/2026-CEC))

<p align="center">
  <strong>English</strong> |
  <a href="README.zh-CN.md">中文</a> |
  <a href="README.ja.md">日本語</a>
</p>

RDE26-MOP is the RDE26 entry for bound-constrained multi-objective optimization, one of the four tracks of the IEEE WCCI 2026 / IEEE CEC 2026 Competition on Numerical Optimization Considering Accuracy and Speed. It ranked first in this track (Competition 3-BC-MOPs).

Scores, requirements and the other tracks: see the [top-level README](../README.md).

## Problem and settings

- Test set: 10 multi-objective problems, MaOP1-MaOP10, with their reference Pareto fronts (`benchmark/SEC2018_MaOP_M3_D10/POF`).
- Defaults of `main.py`: M = 3 objectives, D = 10 decision variables, population size N = 100, 30 independent runs per problem, 100,000 function evaluations per run. The built-in default of the C++ program is D = 7; `main.py` always passes its own value.
- The dimension, number of objectives and budget of a competition are set by the organizers. To match a call, pass them as options, for example `--D 7 --M 3 --maxFE 100000`.
- C++17 and `make` for the solver, Python 3 for the run script (no third-party package).

## Contents

- `benchmark/SEC2018_MaOP_M3_D10/`: `MaOP/` (the ten problems), `POF/` (reference fronts), `PROBLEM.m`, `SOLUTION.m`.
- `code/rde26-mop/main.py`: build and run script. `code/rde26-mop/src/`: solver source and `Makefile`.

## Run

Run from the repository root (or give the full path):

```bash
python3 RDE26_MOP/code/rde26-mop/main.py    # full run with the default settings
python3 RDE26_MOP/code/rde26-mop/main.py --problems 1 --runs 1 --maxFE 5000    # short check
```

## Options of `main.py`

| Option | Default | Meaning |
| --- | --- | --- |
| `--problems` | `1:10` | problems to run, as a list or a range |
| `--runs` | `30` | runs per problem |
| `--seed-base` | `2025170000` | base of the random seeds |
| `--N` | `100` | population size |
| `--M` | `3` | number of objectives |
| `--D` | `10` | number of decision variables |
| `--maxFE` | `100000` | budget per run (evaluations) |
| `--save` | `501` | number of saved populations plus one (500 IGD values per run) |
| `--pf-root` | `benchmark/SEC2018_MaOP_M3_D10/POF` | folder with the reference Pareto fronts |
| `--out-dir` | `outputs/rde26-mop_<time>` | output folder |
| `--paper-id` | `RDE26` | prefix of the output file names |
| `--skip-build` | `off` | reuse the compiled program |

## Output

One file per problem, `RDE26_CEC26_MOP_F<k>.txt`: 500 rows with one tab-separated column per run, holding the IGD of the saved populations. The raw traces are kept in `_raw/`.
