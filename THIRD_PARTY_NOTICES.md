# Third-party material

The MIT License of this repository covers the solvers and run scripts written by the authors listed in [LICENSE](LICENSE). The files below were written by other authors. They are included only so that the solvers can be built and run as in the competition; they stay under the terms of their original authors, and the MIT License does not apply to them. Please cite the original sources when you use them.

| Where | What | Origin |
| --- | --- | --- |
| `RDE26_SOP/benchmark/CEC2017_SOP/`, `RDE26_SOP/code/rde26-sop/src/` (evaluator and `input_data`) | CEC 2017 real-parameter single-objective test suite | N. H. Awad, M. Z. Ali, J. J. Liang, B. Y. Qu and P. N. Suganthan, CEC 2017 competition on single-objective real-parameter optimization (file header: Noor Awad, 2013) |
| `RDE26_CSOP/benchmark/CEC2017_CSOP/`, `RDE26_CSOP/code/rde26-csop/` (`cec17_test_COP.cpp`, `inputData/`) | CEC 2017 constrained real-parameter test suite | G. Wu, R. Mallipeddi and P. N. Suganthan, CEC 2017 competition on constrained real-parameter optimization (file header: Guohua Wu, 2016) |
| `RDE26_MOP/benchmark/SEC2018_MaOP_M3_D10/` | MaOP1-MaOP10 test problems and reference Pareto fronts | as provided for the competition (MATLAB, PlatEMO format) |
| `RDE26_CMOP/benchmark/SDC/`, `RDE26_CMOP/benchmark/cpp/` | SDC1-SDC15 constrained multi-objective test suite (MATLAB) and a C++ port of it | as provided for the competition; the C++ port follows the MATLAB suite |
| `PROBLEM.m`, `SOLUTION.m`, `UserProblem.m` and utility files under `benchmark/` | base classes of the PlatEMO platform | PlatEMO, https://github.com/BIMK/PlatEMO |

If you are the author of one of these files and want the notice changed, please open an issue.
