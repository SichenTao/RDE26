# Third-party material

The MIT License of this repository covers the solvers and run scripts written by the authors listed in [LICENSE](LICENSE). The files below were written by other authors. They are included only so that the solvers can be built and run as in the competition; they stay under the terms of their original authors, and the MIT License does not apply to them. Please cite the original sources when you use them.

| Where | What | Origin |
| --- | --- | --- |
| `RDE26_SOP/benchmark/CEC2017_SOP/`, `RDE26_SOP/code/rde26-sop/src/` (evaluator and `input_data`) | CEC 2017 real-parameter single-objective test suite | N. H. Awad, M. Z. Ali, J. J. Liang, B. Y. Qu and P. N. Suganthan, CEC 2017 competition on single-objective real-parameter optimization (file header: Noor Awad, 2013) |
| `RDE26_CSOP/benchmark/CEC2017_CSOP/`, `RDE26_CSOP/code/rde26-csop/` (`cec17_test_COP.cpp`, `inputData/`) | CEC 2017 constrained real-parameter test suite | G. Wu, R. Mallipeddi and P. N. Suganthan, CEC 2017 competition on constrained real-parameter optimization (file header: Guohua Wu, 2016) |
| `RDE26_MOP/benchmark/SEC2018_MaOP_M3_D10/` | MaOP1-MaOP10 test problems and reference Pareto fronts | as provided for the competition (MATLAB, PlatEMO format) |
| `RDE26_CMOP/benchmark/SDC/Multi-objective optimization/SDC/`, `RDE26_CMOP/benchmark/cpp/` | SDC1-SDC15 constrained multi-objective test suite (MATLAB) and a C++ port of it | K. Qiao, J. Liang, K. Yu, C. Yue, H. Lin, D. Zhang and B. Qu, "Evolutionary constrained multiobjective optimization: scalable high-dimensional constraint benchmarks and algorithm," IEEE Transactions on Evolutionary Computation, 2023 (the reference given in the file headers); the C++ port follows the MATLAB suite |
| `PROBLEM.m`, `SOLUTION.m`, `UserProblem.m` and other PlatEMO files under `benchmark/` | base classes of the PlatEMO platform | PlatEMO, https://github.com/BIMK/PlatEMO. The file headers carry the PlatEMO notice: PlatEMO may be used for research purposes, and publications that use it should acknowledge PlatEMO and cite Y. Tian, R. Cheng, X. Zhang and Y. Jin, "PlatEMO: A MATLAB platform for evolutionary multi-objective optimization [educational forum]," IEEE Computational Intelligence Magazine, 12(4): 73-87, 2017 |

If you are the author of one of these files and want the notice changed, please open an issue.
