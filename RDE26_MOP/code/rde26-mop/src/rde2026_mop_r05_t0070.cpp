#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

static const double PI = 3.141592653589793238462643383279502884;
static const double EPS = 1.0e-12;

struct Config {
    int problem = 1;
    int run = 1;
    int N = 100;
    int M = 3;
    int D = 7;
    int maxFE = 100000;
    int save = 501;
    unsigned seed = 2025170000u;
    string algName = "RDE2026MOPTPE_R05_T0070";
    string outDir = "results_cpp";
    string pfRoot = "../../../benchmark/SEC2018_MaOP_M3_D10/POF";
    string evalFile;
    string evalOut;
};

struct Rng {
    mt19937 gen;
    uniform_real_distribution<double> uni;
    normal_distribution<double> norm;
    cauchy_distribution<double> cauchy;

    explicit Rng(unsigned seed) : gen(seed), uni(0.0, 1.0), norm(0.0, 1.0), cauchy(0.0, 1.0) {}

    int randint(int n) {
        if (n <= 0) return 0;
        uniform_int_distribution<int> dist(0, n - 1);
        return dist(gen);
    }
    double rand01() { return uni(gen); }
    double randn() { return norm(gen); }
    double cauchy01() { return cauchy(gen); }
    double matlabCauchy01() { return tan(PI * (rand01() - 0.5)); }
};

struct Solution {
    vector<double> dec;
    vector<double> obj;
};

struct FrontState {
    double frontRatio = 0.0;
    double nnCV = 0.0;
    double nnMean = 0.0;
    double keepRatio = 1.0;
    double spanRatio = 1.0;
    bool needsPolish = false;
};

struct ESDState {
    int overflowCount = 0;
    int convergenceCount = 0;
    int activations = 0;
    int cfesActivations = 0;
    int dfesActivations = 0;
    vector<double> ideal;
};

struct RootParams {
    double traceArchiveCapFactor = 3.0;
    double guideArchiveCapFactor = 3.0;
    double archiveOutputProgress = 0.82;
    double archiveRefreshFEFrac = 1.0 / 40.0;
    double espStartProgress = 0.25;
    int espOverflowThreshold = 18;
    double espCoverageFloorEarly = 0.965;
    double espCoverageFloorLate = 0.925;
    double injectStartProgress = 0.45;
    double injectBiasThreshold = 0.10;
    double injectQBase = 0.05;
    double injectQGain = 0.10;
    double elitePbestBase = 0.14;
    double elitePbestDecay = 0.07;
    double eliteArchiveBase = 0.18;
    double eliteArchiveGain = 0.10;
    double eliteFarchBase = 0.18;
    double eliteFarchGain = 0.12;
};

struct RootOptions {
    bool useESD = false;
    bool useESP = false;
    bool useDirectional = true;
    bool useExplore = true;
    bool useEliteBranch = true;
    bool useEliteOperator = true;
    bool useRefine = true;
    bool useMemory = true;
    bool useRepair = true;
    bool useArchiveUpdate = true;
    bool usePolish = true;
    bool useStableOutput = true;
    bool useSidecar = true;
    int archiveMode = 0;
    bool useArchiveOutput = false;
    RootParams params;
};

struct Problem {
    int id;
    int N;
    int M;
    int D;
    int maxFE;
    int FE = 0;
    vector<double> lower;
    vector<double> upper;

    Problem(int problemId, int n, int m, int d, int maxFe)
        : id(problemId), N(n), M(m), D(d), maxFE(maxFe), lower(d, 0.0), upper(d, 1.0) {}

    Solution evaluateOne(vector<double> x) {
        for (int i = 0; i < D; ++i) {
            if (x[i] < lower[i]) x[i] = lower[i];
            if (x[i] > upper[i]) x[i] = upper[i];
        }
        Solution s;
        s.dec = std::move(x);
        s.obj = calObj(s.dec);
        ++FE;
        return s;
    }

    vector<Solution> evaluate(const vector<vector<double>>& xs) {
        vector<Solution> out;
        out.reserve(xs.size());
        for (const auto& x : xs) out.push_back(evaluateOne(x));
        return out;
    }

    vector<double> calObj(const vector<double>& x) const {
        vector<double> f(M, 0.0);
        const int nobj = M;
        const int nvar = D;
        double tmp = 1.0;
        for (int j = 0; j < nobj - 1; ++j) tmp *= sin(0.5 * PI * x[j]);

        if (id == 1) {
            double g = 0.0;
            for (int n = nobj; n <= nvar; ++n) {
                double z = x[n - 1] - 0.5;
                g += z * z + (1.0 - cos(20.0 * PI * z));
            }
            g /= nvar;
            double prod = 1.0;
            for (int m = nobj; m >= 1; --m) {
                int out = m - 1;
                int idx = nobj - m;
                if (m > 1) {
                    f[out] = (1.0 + g) * (1.0 - prod * (1.0 - x[idx]));
                    prod *= x[idx];
                } else {
                    f[out] = (1.0 + g) * (1.0 - prod);
                }
                f[out] *= 0.1 + 10.0 * m;
            }
            return f;
        }

        if (id == 2) {
            double g = 0.0;
            for (int n = nobj; n <= nvar; ++n) {
                double target = (n % 5 == 0) ? tmp : 0.5;
                double z = x[n - 1] - target;
                g += z * z;
            }
            g *= 200.0;
            double prod = 1.0;
            for (int m = nobj; m >= 1; --m) {
                int out = m - 1;
                double p = pow(2.0, (m % 2) + 1);
                if (m == nobj) {
                    f[out] = (1.0 + g) * pow(sin(0.5 * x[0] * PI), p);
                } else if (m >= 2) {
                    prod *= cos(0.5 * PI * x[nobj - m - 1]);
                    f[out] = (1.0 + g) * pow(prod * sin(0.5 * PI * x[nobj - m]), p);
                } else {
                    f[out] = (1.0 + g) * pow(prod * cos(0.5 * PI * x[nobj - 2]), p);
                }
            }
            return f;
        }

        if (id == 3 || id == 4) {
            double g = 0.0;
            for (int n = nobj; n <= nvar; ++n) {
                double target = (n % 5 == 0) ? tmp : 0.5;
                double z = fabs(x[n - 1] - target);
                if (id == 3) {
                    g += n * pow(z, 0.1);
                } else {
                    g += 2.0 * sin(x[0] * PI) * (fabs(-0.9 * z * z) + pow(z, 0.6));
                }
            }
            if (id == 4) g *= 10.0;
            double prod = 1.0;
            for (int m = nobj; m >= 1; --m) {
                int out = m - 1;
                if (m == nobj) {
                    f[out] = (1.0 + g) * sin(0.5 * x[0] * PI);
                } else if (m >= 2) {
                    prod *= cos(0.5 * PI * x[nobj - m - 1]);
                    f[out] = (1.0 + g) * prod * sin(0.5 * PI * x[nobj - m]);
                } else {
                    f[out] = (1.0 + g) * prod * cos(0.5 * PI * x[nobj - 2]);
                }
            }
            return f;
        }

        if (id == 5 || id == 6) {
            vector<double> g(nobj, 0.0);
            for (int m = 1; m <= nobj; ++m) {
                if (m <= 3) {
                    double term = 0.0;
                    if (id == 5) term = max(0.0, -1.4 * cos(2.0 * x[0] * PI));
                    else term = max(0.0, 1.4 * sin(4.0 * x[0] * PI));
                    for (int j = 2; j < nvar; ++j) {
                        double z = x[j] - x[0] * x[1];
                        term += z * z;
                    }
                    g[m - 1] = 10.0 * term;
                } else {
                    double z = x[m - 1] - x[0] * x[1];
                    g[m - 1] = 10.0 * (exp(z * z) - 1.0);
                }
            }
            if (id == 5) {
                double a1 = cos(0.5 * PI * x[0]) * cos(0.5 * PI * x[1]);
                double a2 = cos(0.5 * PI * x[0]) * sin(0.5 * PI * x[1]);
                double a3 = sin(0.5 * PI * x[0]);
                f[0] = (1.0 + g[0]) * a1;
                f[1] = 4.0 * (1.0 + g[1]) * a2;
                f[2] = (1.0 + g[2]) * a3;
            } else {
                double a1 = x[0] * x[1];
                double a2 = (1.0 - x[1]) * x[0];
                double a3 = 1.0 - x[0];
                f[0] = (1.0 + g[0]) * a1;
                f[1] = 2.0 * (1.0 + g[1]) * a2;
                f[2] = 6.0 * (1.0 + g[2]) * a3;
            }
            return f;
        }

        if (id >= 7 && id <= 10) {
            double g = 0.0;
            for (int n = nobj; n <= nvar; ++n) {
                double target = (n % 5 == 0) ? tmp : 0.5;
                double z = x[n - 1] - target;
                g += z * z;
            }
            g *= 100.0;
            vector<double> alpha(nobj, 0.0);
            double tau = sqrt(2.0) / 2.0;
            alpha[0] = -pow(2.0 * x[0] - 1.0, 3.0) + 1.0;
            int T = (nobj - 1) / 2;
            for (int j = 1; j <= T; ++j) {
                double y = x[j];
                double z = 2.0 * (2.0 * y - floor(2.0 * y)) - 1.0;
                double p;
                if (id == 7 || id == 9) p = 0.5 + x[0];
                else if (id == 8) p = 1.0 - 0.5 * sin(4.0 * PI * x[0]);
                else p = (y < 0.5) ? 0.5 + x[0] : 1.5 - x[0];
                if (id == 7 || id == 8) z = 2.0 * y - 1.0;
                alpha[2 * j - 1] = x[0] + 2.0 * y * tau + tau * pow(fabs(z), p);
                alpha[2 * j] = x[0] - (2.0 * y - 2.0) * tau + tau * pow(fabs(z), p);
            }
            if (nobj % 2 == 0) alpha[nobj - 1] = 1.0 - alpha[0];
            for (int m = 0; m < nobj; ++m) f[m] = (1.0 + g) * alpha[m];
            return f;
        }

        return f;
    }
};

static Config parseArgs(int argc, char** argv) {
    Config cfg;
    unordered_map<string, string*> sargs = {
        {"--out-dir", &cfg.outDir}, {"--pf-root", &cfg.pfRoot}, {"--alg-name", &cfg.algName},
        {"--eval-file", &cfg.evalFile}, {"--eval-out", &cfg.evalOut}
    };
    for (int i = 1; i < argc; ++i) {
        string k = argv[i];
        if (i + 1 >= argc) break;
        string v = argv[++i];
        if (k == "--problem") cfg.problem = stoi(v);
        else if (k == "--run") cfg.run = stoi(v);
        else if (k == "--N") cfg.N = stoi(v);
        else if (k == "--M") cfg.M = stoi(v);
        else if (k == "--D") cfg.D = stoi(v);
        else if (k == "--maxFE") cfg.maxFE = stoi(v);
        else if (k == "--save") cfg.save = stoi(v);
        else if (k == "--seed") cfg.seed = static_cast<unsigned>(stoul(v));
        else if (sargs.count(k)) *sargs[k] = v;
    }
    return cfg;
}

static RootParams fixedR05T0070RootParams() {
    RootParams p;
    p.traceArchiveCapFactor = 2.0;
    p.guideArchiveCapFactor = 3.25;
    p.archiveOutputProgress = 0.78;
    p.archiveRefreshFEFrac = 0.025;
    p.espStartProgress = 0.17;
    p.espOverflowThreshold = 37;
    p.espCoverageFloorEarly = 0.985;
    p.espCoverageFloorLate = 0.895;
    p.injectStartProgress = 0.44;
    p.injectBiasThreshold = 0.05;
    p.injectQBase = 0.09;
    p.injectQGain = 0.07;
    p.elitePbestBase = 0.08;
    p.elitePbestDecay = 0.04;
    p.eliteArchiveBase = 0.16;
    p.eliteArchiveGain = 0.12;
    p.eliteFarchBase = 0.20;
    p.eliteFarchGain = 0.14;
    return p;
}

static void append(vector<Solution>& dst, const vector<Solution>& src) {
    dst.insert(dst.end(), src.begin(), src.end());
}

static double sqrDist(const vector<double>& a, const vector<double>& b) {
    double s = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        double d = a[i] - b[i];
        s += d * d;
    }
    return s;
}

static bool dominates(const vector<double>& a, const vector<double>& b) {
    bool strict = false;
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i] > b[i]) return false;
        if (a[i] < b[i]) strict = true;
    }
    return strict;
}

static vector<int> firstFrontIndices(const vector<Solution>& pop) {
    int n = static_cast<int>(pop.size());
    vector<int> front;
    front.reserve(n);
    for (int i = 0; i < n; ++i) {
        bool dominated = false;
        for (int j = 0; j < n && !dominated; ++j) {
            if (i != j && dominates(pop[j].obj, pop[i].obj)) dominated = true;
        }
        if (!dominated) front.push_back(i);
    }
    return front;
}

static vector<int> firstFrontIndicesObj(const vector<vector<double>>& obj) {
    int n = static_cast<int>(obj.size());
    vector<int> front;
    front.reserve(n);
    for (int i = 0; i < n; ++i) {
        bool dominated = false;
        for (int j = 0; j < n && !dominated; ++j) {
            if (i != j && dominates(obj[j], obj[i])) dominated = true;
        }
        if (!dominated) front.push_back(i);
    }
    return front;
}

static vector<Solution> subset(const vector<Solution>& pop, const vector<int>& idx) {
    vector<Solution> out;
    out.reserve(idx.size());
    for (int i : idx) out.push_back(pop[i]);
    return out;
}

static vector<vector<double>> normalizeObj(const vector<Solution>& pop, vector<double>* spanOut = nullptr) {
    if (pop.empty()) return {};
    int n = static_cast<int>(pop.size());
    int m = static_cast<int>(pop[0].obj.size());
    vector<double> mn(m, numeric_limits<double>::infinity()), mx(m, -numeric_limits<double>::infinity());
    for (const auto& s : pop) {
        for (int j = 0; j < m; ++j) {
            mn[j] = min(mn[j], s.obj[j]);
            mx[j] = max(mx[j], s.obj[j]);
        }
    }
    vector<double> span(m);
    for (int j = 0; j < m; ++j) {
        span[j] = mx[j] - mn[j];
        if (span[j] == 0.0) span[j] = 1.0;
    }
    if (spanOut) *spanOut = span;
    vector<vector<double>> out(n, vector<double>(m, 0.0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) out[i][j] = (pop[i].obj[j] - mn[j]) / span[j];
    }
    return out;
}

static vector<double> nearestDistances(const vector<vector<double>>& x, int k) {
    int n = static_cast<int>(x.size());
    vector<double> out(n, 0.0);
    if (n <= 1) return out;
    k = max(1, min(k, n - 1));
    vector<double> d(n - 1);
    for (int i = 0; i < n; ++i) {
        int p = 0;
        for (int j = 0; j < n; ++j) if (j != i) d[p++] = sqrt(sqrDist(x[i], x[j]));
        nth_element(d.begin(), d.begin() + k - 1, d.end());
        out[i] = d[k - 1];
    }
    return out;
}

static FrontState estimateState(const vector<Solution>& pop, int N, double progress) {
    FrontState st;
    st.frontRatio = static_cast<double>(pop.size()) / max(1, N);
    if (pop.size() < 3) return st;
    vector<double> span;
    auto norm = normalizeObj(pop, &span);
    auto nn = nearestDistances(norm, 1);
    double mean = accumulate(nn.begin(), nn.end(), 0.0) / nn.size();
    double var = 0.0;
    for (double v : nn) var += (v - mean) * (v - mean);
    var /= max<size_t>(1, nn.size() - 1);
    st.nnMean = mean;
    st.nnCV = sqrt(var) / max(mean, EPS);
    auto minmax = minmax_element(span.begin(), span.end());
    st.spanRatio = *minmax.second / max(*minmax.first, EPS);
    st.needsPolish = progress >= 0.90 && (st.nnCV >= 0.55 || st.frontRatio >= 1.9);
    return st;
}

static vector<vector<double>> uniformPoints(int N, int M) {
    auto comb = [](int n, int k) -> long long {
        if (k < 0 || k > n) return 0;
        k = min(k, n - k);
        long long r = 1;
        for (int i = 1; i <= k; ++i) r = r * (n - k + i) / i;
        return r;
    };
    auto lattice3 = [](int H) {
        vector<vector<double>> out;
        if (H <= 0) return out;
        for (int a = 0; a <= H; ++a) {
            for (int b = 0; b <= H - a; ++b) {
                int c = H - a - b;
                out.push_back({a / static_cast<double>(H), b / static_cast<double>(H), c / static_cast<double>(H)});
            }
        }
        return out;
    };
    if (M != 3) {
        vector<vector<double>> W;
        for (int i = 0; i < N; ++i) {
            vector<double> w(M, 1.0 / M);
            w[i % M] = 0.75;
            double rest = (1.0 - w[i % M]) / (M - 1);
            for (int j = 0; j < M; ++j) if (j != i % M) w[j] = rest;
            W.push_back(w);
        }
        return W;
    }

    int H1 = 1;
    while (comb(H1 + M, M - 1) <= N) ++H1;
    vector<vector<double>> W = lattice3(H1);
    if (H1 < M) {
        int H2 = 0;
        while (comb(H1 + M - 1, M - 1) + comb(H2 + M, M - 1) <= N) ++H2;
        if (H2 > 0) {
            vector<vector<double>> W2 = lattice3(H2);
            for (auto& w : W2) {
                for (double& v : w) v = v / 2.0 + 1.0 / (2.0 * M);
                W.push_back(w);
            }
        }
    }
    for (auto& w : W) for (double& v : w) v = max(v, 1.0e-6);
    return W;
}

static vector<double> normalize01(const vector<double>& x) {
    vector<double> y(x.size(), 1.0);
    if (x.empty()) return y;
    auto [mnIt, mxIt] = minmax_element(x.begin(), x.end());
    double span = *mxIt - *mnIt;
    if (span < EPS) return y;
    for (size_t i = 0; i < x.size(); ++i) y[i] = (x[i] - *mnIt) / span;
    return y;
}

static vector<int> pickReferenceRepresentativesNorm(const vector<vector<double>>& normObj, int quota, int M) {
    vector<int> picked;
    if (normObj.empty() || quota <= 0) return picked;
    auto W = uniformPoints(max(quota, M), M);
    for (auto& w : W) {
        double wn = 0.0;
        for (double& v : w) {
            v = max(v, EPS);
            wn += v * v;
        }
        wn = sqrt(max(wn, EPS));
        for (double& v : w) v /= wn;
    }
    vector<vector<double>> y = normObj;
    for (auto& row : y) {
        double s = 0.0;
        for (double& v : row) {
            v = max(v, EPS);
            s += v;
        }
        for (double& v : row) v /= max(s, EPS);
    }
    vector<int> nearest;
    nearest.reserve(W.size());
    for (const auto& w : W) {
        int best = 0;
        double bestCos = -numeric_limits<double>::infinity();
        for (int i = 0; i < static_cast<int>(y.size()); ++i) {
            double yn = sqrt(inner_product(y[i].begin(), y[i].end(), y[i].begin(), 0.0));
            double c = inner_product(y[i].begin(), y[i].end(), w.begin(), 0.0) / max(yn, EPS);
            if (c > bestCos) {
                bestCos = c;
                best = i;
            }
        }
        if (find(picked.begin(), picked.end(), best) == picked.end()) picked.push_back(best);
    }
    if (static_cast<int>(picked.size()) > quota) picked.resize(quota);
    return picked;
}

static vector<double> calFitness(const vector<Solution>& pop, double kappa,
                                vector<double>* Iout = nullptr, vector<double>* Cout = nullptr) {
    int n = static_cast<int>(pop.size());
    if (n == 0) return {};
    int m = static_cast<int>(pop[0].obj.size());
    auto norm = normalizeObj(pop);
    vector<double> I(n * n, 0.0), C(n, 1.0), fit(n, 1.0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double v = -numeric_limits<double>::infinity();
            for (int k = 0; k < m; ++k) v = max(v, norm[i][k] - norm[j][k]);
            I[i * n + j] = v;
        }
    }
    for (int j = 0; j < n; ++j) {
        double c = 0.0;
        for (int i = 0; i < n; ++i) c = max(c, fabs(I[i * n + j]));
        C[j] = (c == 0.0) ? 1.0 : c;
    }
    for (int j = 0; j < n; ++j) {
        double s = 0.0;
        for (int i = 0; i < n; ++i) s += -exp(-I[i * n + j] / C[j] / kappa);
        fit[j] = s + 1.0;
    }
    if (Iout) *Iout = std::move(I);
    if (Cout) *Cout = std::move(C);
    return fit;
}

static vector<Solution> environmentalSelection(const vector<Solution>& pop, int N, double kappa) {
    int n = static_cast<int>(pop.size());
    if (n <= N) return pop;
    vector<int> next(n);
    iota(next.begin(), next.end(), 0);
    vector<double> I, C;
    vector<double> fit = calFitness(pop, kappa, &I, &C);
    while (static_cast<int>(next.size()) > N) {
        int pos = 0;
        for (int i = 1; i < static_cast<int>(next.size()); ++i) {
            if (fit[next[i]] < fit[next[pos]]) pos = i;
        }
        int x = next[pos];
        for (int j = 0; j < n; ++j) fit[j] += exp(-I[x * n + j] / C[j] / kappa);
        next.erase(next.begin() + pos);
    }
    return subset(pop, next);
}

static vector<Solution> selectionExact(vector<Solution> pc, int N, Rng& rng, bool presorted = false) {
    if (!presorted) pc = subset(pc, firstFrontIndices(pc));
    shuffle(pc.begin(), pc.end(), rng.gen);
    int n = static_cast<int>(pc.size());
    if (n <= N) return pc;
    auto norm = normalizeObj(pc);
    vector<vector<double>> d(n, vector<double>(n, numeric_limits<double>::infinity()));
    vector<double> allk;
    allk.reserve(n);
    for (int i = 0; i < n; ++i) {
        vector<double> row;
        row.reserve(n - 1);
        for (int j = 0; j < n; ++j) if (i != j) {
            d[i][j] = sqrt(sqrDist(norm[i], norm[j]));
            row.push_back(d[i][j]);
        }
        int k = min(3, static_cast<int>(row.size()));
        nth_element(row.begin(), row.begin() + k - 1, row.end());
        allk.push_back(row[k - 1]);
    }
    double r = accumulate(allk.begin(), allk.end(), 0.0) / allk.size();
    r = max(r, EPS);
    vector<int> alive(n, 1);
    int aliveCount = n;
    while (aliveCount > N) {
        double bestScore = -1.0;
        int worst = -1;
        for (int i = 0; i < n; ++i) if (alive[i]) {
            double prod = 1.0;
            for (int j = 0; j < n; ++j) if (alive[j] && i != j) {
                prod *= min(d[i][j] / r, 1.0);
                if (prod == 0.0) break;
            }
            double score = 1.0 - prod;
            if (score > bestScore) {
                bestScore = score;
                worst = i;
            }
        }
        alive[worst] = 0;
        --aliveCount;
    }
    vector<Solution> out;
    out.reserve(N);
    for (int i = 0; i < n; ++i) if (alive[i]) out.push_back(pc[i]);
    return out;
}

static vector<Solution> preReduceFirstFront(const vector<Solution>& pc, int keepN, double progress, const FrontState& st) {
    int n = static_cast<int>(pc.size());
    int m = static_cast<int>(pc[0].obj.size());
    auto norm = normalizeObj(pc);
    vector<int> keep(n, 0);
    for (int k = 0; k < m; ++k) {
        int mn = 0, mx = 0;
        for (int i = 1; i < n; ++i) {
            if (norm[i][k] < norm[mn][k]) mn = i;
            if (norm[i][k] > norm[mx][k]) mx = i;
        }
        keep[mn] = keep[mx] = 1;
    }
    vector<double> score(n, 0.0), conv(n, 0.0), crowd(n, 0.0), nearScore(n, 0.0);
    vector<vector<double>> sortedDist(n);
    for (int i = 0; i < n; ++i) {
        vector<double> dist;
        dist.reserve(n - 1);
        for (int j = 0; j < n; ++j) if (i != j) dist.push_back(sqrt(sqrDist(norm[i], norm[j])));
        sort(dist.begin(), dist.end());
        sortedDist[i] = dist;
        int nnK = min(3, static_cast<int>(dist.size()));
        for (int t = 0; t < nnK; ++t) nearScore[i] += dist[t];
        nearScore[i] /= max(1, nnK);
        int crowdK = min(m, static_cast<int>(dist.size()));
        for (int t = 0; t < crowdK; ++t) crowd[i] += dist[t];
        conv[i] = 1.0 - accumulate(norm[i].begin(), norm[i].end(), 0.0) / m;
        score[i] = 0.62 * nearScore[i] + 0.23 * crowd[i] + 0.15 * conv[i];
    }
    int remaining = keepN - accumulate(keep.begin(), keep.end(), 0);
    if (remaining > 0) {
        vector<int> candidates;
        for (int i = 0; i < n; ++i) if (!keep[i]) candidates.push_back(i);
        double refFrac = 0.10 + 0.16 * min(max(st.nnCV - 0.25, 0.0), 1.0) +
                         0.06 * min(max(st.frontRatio - 1.2, 0.0), 1.0);
        int refQuota = min(remaining, max(0, static_cast<int>(round(refFrac * keepN))));
        if (refQuota > 0 && !candidates.empty()) {
            vector<vector<double>> candNorm;
            candNorm.reserve(candidates.size());
            for (int id : candidates) candNorm.push_back(norm[id]);
            vector<int> picked = pickReferenceRepresentativesNorm(candNorm, refQuota, m);
            for (int local : picked) keep[candidates[local]] = 1;
            remaining = keepN - accumulate(keep.begin(), keep.end(), 0);
        }
        double convFrac = min(0.10, max(0.0, 0.10 * (progress - 0.60) + 0.08 * (0.45 - st.nnCV)));
        int convQuota = min(remaining, max(0, static_cast<int>(round(convFrac * keepN))));
        if (convQuota > 0) {
            candidates.clear();
            for (int i = 0; i < n; ++i) if (!keep[i]) candidates.push_back(i);
            sort(candidates.begin(), candidates.end(), [&](int a, int b) { return conv[a] > conv[b]; });
            for (int i = 0; i < min(convQuota, static_cast<int>(candidates.size())); ++i) keep[candidates[i]] = 1;
            remaining = keepN - accumulate(keep.begin(), keep.end(), 0);
        }
        if (remaining > 0) {
            candidates.clear();
            for (int i = 0; i < n; ++i) if (!keep[i]) candidates.push_back(i);
            sort(candidates.begin(), candidates.end(), [&](int a, int b) { return score[a] > score[b]; });
            for (int i = 0; i < min(remaining, static_cast<int>(candidates.size())); ++i) keep[candidates[i]] = 1;
        }
    }
    vector<Solution> out;
    out.reserve(keepN);
    for (int i = 0; i < n; ++i) if (keep[i]) out.push_back(pc[i]);
    return out;
}

static pair<vector<Solution>, FrontState> selectionTwoStage(vector<Solution> pc, int N, double progress, double eps, Rng& rng) {
    pc = subset(pc, firstFrontIndices(pc));
    FrontState st = estimateState(pc, N, progress);
    if (static_cast<int>(pc.size()) <= N) return {pc, st};
    double preRatio = 1.45 + 0.16 * min(st.frontRatio, 3.0) + 0.32 * min(st.nnCV, 1.0) + 0.10 * progress;
    preRatio = min(2.10, max(1.55, preRatio));
    int preN = min(static_cast<int>(pc.size()), max(static_cast<int>(round(preRatio * N)), N + 24));
    st.keepRatio = static_cast<double>(preN) / pc.size();
    st.needsPolish = progress >= 0.90 && (st.nnCV >= 0.55 || st.frontRatio >= 1.9);
    if (static_cast<int>(pc.size()) > preN) pc = preReduceFirstFront(pc, preN, progress, st);
    pc = selectionExact(pc, N, rng, true);
    return {pc, st};
}

static pair<vector<Solution>, FrontState> selectionDirectional(vector<Solution> pc, int N, double progress, double eps, Rng& rng) {
    pc = subset(pc, firstFrontIndices(pc));
    auto norm = normalizeObj(pc);
    FrontState st = estimateState(pc, N, progress);
    int n = static_cast<int>(pc.size());
    if (n <= N) return {pc, st};
    int m = static_cast<int>(pc[0].obj.size());
    auto W = uniformPoints(max(N, m), m);
    vector<vector<double>> dir = norm;
    for (auto& row : dir) {
        double s = 0.0;
        for (double& v : row) {
            v = max(v, EPS);
            s += v;
        }
        for (double& v : row) v /= max(s, EPS);
    }
    vector<int> assoc(n, 0);
    vector<double> bestCos(n, 0.0);
    for (int i = 0; i < n; ++i) {
        double dn = sqrt(inner_product(dir[i].begin(), dir[i].end(), dir[i].begin(), 0.0));
        for (int w = 0; w < static_cast<int>(W.size()); ++w) {
            double wn = sqrt(inner_product(W[w].begin(), W[w].end(), W[w].begin(), 0.0));
            double c = inner_product(dir[i].begin(), dir[i].end(), W[w].begin(), 0.0) / max(dn * wn, EPS);
            if (w == 0 || c > bestCos[i]) {
                bestCos[i] = c;
                assoc[i] = w;
            }
        }
    }
    auto fit = calFitness(pc, eps);
    auto fitScore = normalize01(fit);
    auto angleScore = normalize01(bestCos);
    vector<double> sumNorm(n, 0.0), convRaw(n, 0.0);
    for (int i = 0; i < n; ++i) sumNorm[i] = accumulate(norm[i].begin(), norm[i].end(), 0.0);
    double maxSum = *max_element(sumNorm.begin(), sumNorm.end());
    for (int i = 0; i < n; ++i) convRaw[i] = maxSum - sumNorm[i];
    auto convScore = normalize01(convRaw);
    vector<int> keep(n, 0);
    for (int k = 0; k < m; ++k) {
        int mn = 0, mx = 0;
        for (int i = 1; i < n; ++i) {
            if (norm[i][k] < norm[mn][k]) mn = i;
            if (norm[i][k] > norm[mx][k]) mx = i;
        }
        keep[mn] = keep[mx] = 1;
    }
    double dirWeight = 0.44 + 0.18 * min(max(st.nnCV - 0.25, 0.0), 1.0);
    double fitWeight = 0.40 + 0.14 * max(progress - 0.55, 0.0);
    double convWeight = max(0.08, 1.0 - dirWeight - fitWeight);
    vector<double> repScore(n);
    for (int i = 0; i < n; ++i) repScore[i] = fitWeight * fitScore[i] + dirWeight * angleScore[i] + convWeight * convScore[i];
    unordered_map<int, int> bestByDir;
    for (int i = 0; i < n; ++i) {
        int a = assoc[i];
        if (!bestByDir.count(a) || repScore[i] > repScore[bestByDir[a]]) bestByDir[a] = i;
    }
    for (auto& kv : bestByDir) keep[kv.second] = 1;
    if (accumulate(keep.begin(), keep.end(), 0) > N) {
        vector<int> idx;
        for (int i = 0; i < n; ++i) if (keep[i]) idx.push_back(i);
        auto out = selectionExact(subset(pc, idx), N, rng, true);
        return {out, st};
    }
    while (accumulate(keep.begin(), keep.end(), 0) < min(N, n)) {
        vector<int> candidates;
        vector<double> d2keep;
        for (int i = 0; i < n; ++i) if (!keep[i]) {
            double dkeep = 1.0;
            bool any = false;
            for (int j = 0; j < n; ++j) if (keep[j]) {
                double d = sqrt(sqrDist(norm[i], norm[j]));
                if (!any || d < dkeep) dkeep = d;
                any = true;
            }
            candidates.push_back(i);
            d2keep.push_back(dkeep);
        }
        if (candidates.empty()) break;
        vector<double> dScore = normalize01(d2keep);
        int best = candidates[0];
        double bestScore = -numeric_limits<double>::infinity();
        for (int c = 0; c < static_cast<int>(candidates.size()); ++c) {
            int i = candidates[c];
            double sc = 0.46 * fitScore[i] + 0.34 * dScore[c] + 0.12 * convScore[i] + 0.08 * angleScore[i];
            if (sc > bestScore) {
                bestScore = sc;
                best = i;
            }
        }
        keep[best] = 1;
    }
    vector<int> idx;
    for (int i = 0; i < n; ++i) if (keep[i]) idx.push_back(i);
    auto out = subset(pc, idx);
    if (static_cast<int>(out.size()) > N) out = selectionExact(out, N, rng, true);
    return {out, st};
}

static pair<vector<Solution>, FrontState> selectionRootAdaptive(vector<Solution> pc, int N, double progress, double eps, Rng& rng) {
    if (progress <= 0.42) return selectionDirectional(std::move(pc), N, progress, eps, rng);
    return selectionTwoStage(std::move(pc), N, progress, eps, rng);
}

static pair<vector<Solution>, FrontState> selectionRootConfigured(vector<Solution> pc, int N, double progress,
                                                                  double eps, Rng& rng, const RootOptions& opt) {
    if (opt.useDirectional) return selectionRootAdaptive(std::move(pc), N, progress, eps, rng);
    return selectionTwoStage(std::move(pc), N, progress, eps, rng);
}

static pair<vector<int>, vector<int>> gnR1R2FastStyle(int N, Rng& rng) {
    vector<int> r1(N), r2(N);
    for (int i = 0; i < N; ++i) r1[i] = rng.randint(N);
    while (true) {
        vector<int> conflict;
        for (int i = 0; i < N; ++i) if (r1[i] == i) conflict.push_back(i);
        if (conflict.empty()) break;
        for (int idx : conflict) r1[idx] = rng.randint(N);
    }
    for (int i = 0; i < N; ++i) r2[i] = rng.randint(N);
    while (true) {
        vector<int> conflict;
        for (int i = 0; i < N; ++i) if (r2[i] == i || r2[i] == r1[i]) conflict.push_back(i);
        if (conflict.empty()) break;
        for (int idx : conflict) r2[idx] = rng.randint(N);
    }
    return {r1, r2};
}

static vector<Solution> operatorDE2026(Problem& prob, const vector<Solution>& pop, const vector<double>& fitness, Rng& rng) {
    int N = static_cast<int>(pop.size());
    int D = prob.D;
    double progress = static_cast<double>(prob.FE) / prob.maxFE;
    vector<int> order(N);
    iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int a, int b) { return fitness[a] < fitness[b]; });
    double pRate = max(0.10, 0.18 - 0.08 * progress);
    int P = max(2, static_cast<int>(round(pRate * N)));
    const double FmBase[3] = {0.6, 0.8, 1.0};
    const double CRmBase[3] = {0.1, 0.2, 1.0};
    vector<int> pbest(N);
    for (int i = 0; i < N; ++i) pbest[i] = order[rng.randint(P)];
    vector<double> F(N), CR(N);
    for (int i = 0; i < N; ++i) F[i] = FmBase[rng.randint(3)];
    for (int i = 0; i < N; ++i) CR[i] = CRmBase[rng.randint(3)];
    for (int i = 0; i < N; ++i) F[i] += 0.08 * (0.5 - progress) + 0.03 * rng.randn();
    for (int i = 0; i < N; ++i) CR[i] += 0.05 * (0.35 - progress) + 0.03 * rng.randn();
    for (int i = 0; i < N; ++i) {
        F[i] = min(max(F[i], 0.45), 1.0);
        CR[i] = min(max(CR[i], 0.05), 1.0);
    }
    auto [r1, r2] = gnR1R2FastStyle(N, rng);
    vector<int> jrand(N);
    for (int i = 0; i < N; ++i) jrand[i] = rng.randint(D);
    vector<vector<char>> crossMask(N, vector<char>(D, 0));
    for (int j = 0; j < D; ++j) for (int i = 0; i < N; ++i) crossMask[i][j] = (rng.rand01() < CR[i]);
    for (int i = 0; i < N; ++i) crossMask[i][jrand[i]] = 1;
    double perturbProb = max(0.08, 0.20 - 0.08 * progress);
    vector<vector<char>> perturbMask(N, vector<char>(D, 0));
    bool anyPerturb = false;
    for (int j = 0; j < D; ++j) {
        for (int i = 0; i < N; ++i) {
            perturbMask[i][j] = (rng.rand01() < perturbProb);
            anyPerturb = anyPerturb || perturbMask[i][j];
        }
    }
    double perturbShift = 0.0;
    if (anyPerturb) {
        double scale = max(0.08, 0.20 - 0.06 * progress);
        perturbShift = scale * rng.matlabCauchy01();
    }
    vector<vector<double>> xs(N, vector<double>(D));
    for (int i = 0; i < N; ++i) {
        xs[i] = pop[i].dec;
        for (int j = 0; j < D; ++j) {
            double v = pop[i].dec[j] + F[i] * (pop[pbest[i]].dec[j] - pop[i].dec[j] + pop[r1[i]].dec[j] - pop[r2[i]].dec[j]);
            if (crossMask[i][j]) xs[i][j] = v;
            if (perturbMask[i][j]) xs[i][j] += perturbShift;
        }
    }
    return prob.evaluate(xs);
}

static vector<Solution> operatorDEHalf2026(Problem& prob, const vector<Solution>& anchors, const vector<Solution>& guides, Rng& rng) {
    int N = static_cast<int>(anchors.size());
    int D = prob.D;
    double progress = static_cast<double>(prob.FE) / prob.maxFE;
    vector<double> F(N), CR(N);
    for (int i = 0; i < N; ++i) F[i] = min(max(0.60 + 0.10 * rng.randn(), 0.25), 0.95);
    for (int i = 0; i < N; ++i) CR[i] = min(max(0.55 - 0.20 * progress + 0.08 * rng.randn(), 0.05), 0.95);
    vector<int> perm(N);
    iota(perm.begin(), perm.end(), 0);
    shuffle(perm.begin(), perm.end(), rng.gen);
    vector<int> jrand(N);
    for (int i = 0; i < N; ++i) jrand[i] = rng.randint(D);
    vector<vector<char>> crossMask(N, vector<char>(D, 0));
    for (int j = 0; j < D; ++j) for (int i = 0; i < N; ++i) crossMask[i][j] = (rng.rand01() < CR[i]);
    for (int i = 0; i < N; ++i) crossMask[i][jrand[i]] = 1;
    vector<vector<double>> xs(N, vector<double>(D));
    for (int i = 0; i < N; ++i) {
        xs[i] = anchors[i].dec;
        for (int j = 0; j < D; ++j) {
            double v = anchors[i].dec[j] + F[i] * (anchors[perm[i]].dec[j] - guides[i].dec[j]);
            if (crossMask[i][j]) xs[i][j] = v;
        }
    }
    return prob.evaluate(xs);
}

static vector<Solution> exploration2026(Problem& prob, const vector<Solution>& population,
                                        const vector<Solution>& popul, int nND, Rng& rng) {
    if (population.size() < 3 || popul.empty()) return {};
    auto pcNorm = normalizeObj(population);
    vector<Solution> refPop = population;
    vector<double> mn(population[0].obj.size(), numeric_limits<double>::infinity());
    vector<double> mx(population[0].obj.size(), -numeric_limits<double>::infinity());
    for (const auto& s : population) for (size_t j = 0; j < s.obj.size(); ++j) {
        mn[j] = min(mn[j], s.obj[j]);
        mx[j] = max(mx[j], s.obj[j]);
    }
    vector<vector<double>> npcNorm(popul.size(), vector<double>(mn.size()));
    for (size_t i = 0; i < popul.size(); ++i) {
        for (size_t j = 0; j < mn.size(); ++j) {
            double span = mx[j] - mn[j];
            if (span == 0.0) span = 1.0;
            npcNorm[i][j] = (popul[i].obj[j] - mn[j]) / span;
        }
    }
    auto nn3 = nearestDistances(pcNorm, min(3, static_cast<int>(pcNorm.size()) - 1));
    double r0 = accumulate(nn3.begin(), nn3.end(), 0.0) / nn3.size();
    double r = max(EPS, static_cast<double>(nND) / prob.N * r0);
    vector<int> sparse;
    for (int i = 0; i < static_cast<int>(pcNorm.size()); ++i) {
        int cover = 0;
        for (const auto& q : npcNorm) if (sqrt(sqrDist(pcNorm[i], q)) <= r) ++cover;
        if (cover <= 1) sparse.push_back(i);
    }
    if (sparse.empty()) return {};
    double progress = static_cast<double>(prob.FE) / prob.maxFE;
    int maxExplore = max(1, static_cast<int>(round(max(0.50, 0.95 - 0.30 * progress) * prob.N)));
    if (static_cast<int>(sparse.size()) > maxExplore) {
        shuffle(sparse.begin(), sparse.end(), rng.gen);
        sparse.resize(maxExplore);
    }
    vector<Solution> anchors, guides;
    anchors.reserve(sparse.size());
    guides.reserve(sparse.size());
    for (int idx : sparse) {
        anchors.push_back(population[idx]);
        guides.push_back(population[rng.randint(static_cast<int>(population.size()))]);
    }
    return operatorDEHalf2026(prob, anchors, guides, rng);
}

static pair<vector<Solution>, int> selectionRDExMOP(vector<Solution> pc, int N, Rng& rng) {
    pc = subset(pc, firstFrontIndices(pc));
    shuffle(pc.begin(), pc.end(), rng.gen);
    int nND = static_cast<int>(pc.size());
    if (static_cast<int>(pc.size()) <= N) return {pc, nND};

    auto norm = normalizeObj(pc);
    int n = static_cast<int>(pc.size());
    vector<vector<double>> d(n, vector<double>(n, numeric_limits<double>::infinity()));
    vector<double> kth(n, 0.0);
    for (int i = 0; i < n; ++i) {
        vector<double> row;
        row.reserve(n - 1);
        for (int j = 0; j < n; ++j) if (i != j) {
            d[i][j] = sqrt(sqrDist(norm[i], norm[j]));
            row.push_back(d[i][j]);
        }
        int k = min(3, static_cast<int>(row.size()));
        nth_element(row.begin(), row.begin() + k - 1, row.end());
        kth[i] = row[k - 1];
    }
    double r = accumulate(kth.begin(), kth.end(), 0.0) / kth.size();
    r = max(r, EPS);
    vector<int> alive(n, 1);
    int aliveCount = n;
    while (aliveCount > N) {
        int worst = -1;
        double worstScore = -numeric_limits<double>::infinity();
        for (int i = 0; i < n; ++i) if (alive[i]) {
            double prod = 1.0;
            for (int j = 0; j < n; ++j) if (alive[j] && i != j) {
                prod *= min(d[i][j] / r, 1.0);
                if (prod == 0.0) break;
            }
            double score = 1.0 - prod;
            if (score > worstScore) {
                worstScore = score;
                worst = i;
            }
        }
        alive[worst] = 0;
        --aliveCount;
    }
    vector<Solution> out;
    out.reserve(N);
    for (int i = 0; i < n; ++i) if (alive[i]) out.push_back(pc[i]);
    return {out, nND};
}

static double angleBetween(const vector<double>& a, const vector<double>& b) {
    double dot = 0.0, na = 0.0, nb = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        double va = max(a[i], EPS);
        double vb = max(b[i], EPS);
        dot += va * vb;
        na += va * va;
        nb += vb * vb;
    }
    double c = dot / max(sqrt(na) * sqrt(nb), EPS);
    c = min(1.0, max(-1.0, c));
    return acos(c);
}

static vector<int> deterministicKMeansLabels(const vector<vector<double>>& norm, int K) {
    int n = static_cast<int>(norm.size());
    int m = static_cast<int>(norm[0].size());
    K = max(1, min(K, n));
    vector<vector<double>> centers;
    centers.reserve(K);
    int first = 0;
    double firstFit = numeric_limits<double>::infinity();
    for (int i = 0; i < n; ++i) {
        double fit = accumulate(norm[i].begin(), norm[i].end(), 0.0);
        if (fit < firstFit) {
            firstFit = fit;
            first = i;
        }
    }
    centers.push_back(norm[first]);
    while (static_cast<int>(centers.size()) < K) {
        int best = 0;
        double bestDist = -1.0;
        for (int i = 0; i < n; ++i) {
            double dmin = numeric_limits<double>::infinity();
            for (const auto& c : centers) dmin = min(dmin, sqrDist(norm[i], c));
            if (dmin > bestDist) {
                bestDist = dmin;
                best = i;
            }
        }
        centers.push_back(norm[best]);
    }

    vector<int> label(n, 0), count(K, 0);
    for (int iter = 0; iter < 10; ++iter) {
        fill(count.begin(), count.end(), 0);
        bool changed = false;
        for (int i = 0; i < n; ++i) {
            int best = 0;
            double bestDist = sqrDist(norm[i], centers[0]);
            for (int k = 1; k < K; ++k) {
                double d = sqrDist(norm[i], centers[k]);
                if (d < bestDist) {
                    bestDist = d;
                    best = k;
                }
            }
            if (label[i] != best) changed = true;
            label[i] = best;
            ++count[best];
        }
        vector<vector<double>> next(K, vector<double>(m, 0.0));
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) next[label[i]][j] += norm[i][j];
        }
        for (int k = 0; k < K; ++k) {
            if (count[k] == 0) {
                int farthest = 0;
                double bestDist = -1.0;
                for (int i = 0; i < n; ++i) {
                    double dmin = numeric_limits<double>::infinity();
                    for (int kk = 0; kk < K; ++kk) if (count[kk] > 0) dmin = min(dmin, sqrDist(norm[i], centers[kk]));
                    if (dmin > bestDist) {
                        bestDist = dmin;
                        farthest = i;
                    }
                }
                centers[k] = norm[farthest];
            } else {
                for (int j = 0; j < m; ++j) next[k][j] /= count[k];
                centers[k] = next[k];
            }
        }
        if (!changed && iter >= 2) break;
    }
    return label;
}

static vector<Solution> esdEnvironmentalSelection(const vector<Solution>& front, int N, ESDState& esd) {
    int n = static_cast<int>(front.size());
    int m = static_cast<int>(front[0].obj.size());
    auto norm = normalizeObj(front);
    if (esd.ideal.size() != static_cast<size_t>(m)) esd.ideal.assign(m, numeric_limits<double>::infinity());
    for (const auto& s : front) for (int j = 0; j < m; ++j) esd.ideal[j] = min(esd.ideal[j], s.obj[j]);

    vector<double> fit(n, 0.0);
    for (int i = 0; i < n; ++i) fit[i] = accumulate(norm[i].begin(), norm[i].end(), 0.0);
    int K = min(6, max(1, min(n, max(2, N / 8))));
    vector<int> label = deterministicKMeansLabels(norm, K);
    vector<int> count(K, 0);
    for (int i = 0; i < n; ++i) {
        ++count[label[i]];
    }
    ++esd.activations;
    ++esd.cfesActivations;
    ++esd.convergenceCount;

    vector<int> protect(n, 0);
    for (int j = 0; j < m; ++j) {
        int best = 0;
        for (int i = 1; i < n; ++i) if (norm[i][j] < norm[best][j]) best = i;
        protect[best] = 1;
    }
    for (int k = 0; k < K; ++k) {
        int best = -1;
        for (int i = 0; i < n; ++i) {
            if (label[i] == k && (best < 0 || fit[i] < fit[best])) best = i;
        }
        if (best >= 0) protect[best] = 1;
    }

    vector<vector<double>> d(n, vector<double>(n, numeric_limits<double>::infinity()));
    vector<vector<double>> angle(n, vector<double>(n, numeric_limits<double>::infinity()));
    vector<double> kth(n, 0.0);
    for (int i = 0; i < n; ++i) {
        vector<double> row;
        row.reserve(n - 1);
        for (int j = 0; j < n; ++j) if (i != j) {
            d[i][j] = sqrt(sqrDist(norm[i], norm[j]));
            angle[i][j] = angleBetween(norm[i], norm[j]);
            row.push_back(d[i][j]);
        }
        int kk = min(3, static_cast<int>(row.size()));
        nth_element(row.begin(), row.begin() + kk - 1, row.end());
        kth[i] = row[kk - 1];
    }
    double r = accumulate(kth.begin(), kth.end(), 0.0) / kth.size();
    r = max(r, EPS);

    vector<int> alive(n, 1), clusterAlive = count;
    int aliveCount = n;
    while (aliveCount > N) {
        int worst = -1;
        double worstScore = -numeric_limits<double>::infinity();
        for (int i = 0; i < n; ++i) if (alive[i]) {
            if (protect[i] && aliveCount - 1 >= N) continue;
            double prod = 1.0;
            double nearestAngle = numeric_limits<double>::infinity();
            for (int j = 0; j < n; ++j) if (alive[j] && i != j) {
                prod *= min(d[i][j] / r, 1.0);
                nearestAngle = min(nearestAngle, angle[i][j]);
                if (prod == 0.0 && nearestAngle < EPS) break;
            }
            double density = 1.0 - prod;
            double target = max(1.0, static_cast<double>(count[label[i]]) * N / n);
            double clusterPressure = max(0.0, clusterAlive[label[i]] / target - 1.0);
            double convBad = fit[i] / max(1, m);
            double angleCrowd = 1.0 - min(nearestAngle / (0.5 * PI), 1.0);
            double score = 0.76 * density + 0.12 * clusterPressure + 0.08 * convBad + 0.04 * angleCrowd;
            if (score > worstScore) {
                worstScore = score;
                worst = i;
            }
        }
        if (worst < 0) {
            for (int i = 0; i < n; ++i) if (alive[i] && (worst < 0 || fit[i] > fit[worst])) worst = i;
        }
        alive[worst] = 0;
        --clusterAlive[label[worst]];
        --aliveCount;
    }
    vector<int> idx;
    idx.reserve(N);
    for (int i = 0; i < n; ++i) if (alive[i]) idx.push_back(i);
    return subset(front, idx);
}

static pair<vector<Solution>, int> selectionRDExMOPESD(vector<Solution> pc, int N, Rng& rng, ESDState& esd) {
    pc = subset(pc, firstFrontIndices(pc));
    int nND = static_cast<int>(pc.size());
    if (nND > N) ++esd.overflowCount;
    else esd.overflowCount = max(0, esd.overflowCount - 1);
    if (nND <= N) return {pc, nND};

    const int alpha = 50;
    if (esd.overflowCount <= alpha) return selectionRDExMOP(std::move(pc), N, rng);
    vector<Solution> out = esdEnvironmentalSelection(pc, N, esd);
    if (static_cast<int>(out.size()) > N) out = selectionExact(out, N, rng, true);
    return {out, nND};
}

struct SelectionQuality {
    double coverage = 0.0;
    double meanSum = 0.0;
    double nnCV = numeric_limits<double>::infinity();
};

static SelectionQuality qualityRelativeToFront(const vector<Solution>& fullFront, const vector<Solution>& pop) {
    SelectionQuality q;
    if (fullFront.empty() || pop.empty()) return q;
    int m = static_cast<int>(fullFront[0].obj.size());
    vector<double> mn(m, numeric_limits<double>::infinity()), mx(m, -numeric_limits<double>::infinity());
    for (const auto& s : fullFront) {
        for (int j = 0; j < m; ++j) {
            mn[j] = min(mn[j], s.obj[j]);
            mx[j] = max(mx[j], s.obj[j]);
        }
    }
    vector<double> span(m, 1.0);
    for (int j = 0; j < m; ++j) span[j] = max(mx[j] - mn[j], EPS);

    vector<vector<double>> norm(pop.size(), vector<double>(m, 0.0));
    vector<double> selMn(m, numeric_limits<double>::infinity()), selMx(m, -numeric_limits<double>::infinity());
    for (int i = 0; i < static_cast<int>(pop.size()); ++i) {
        double sum = 0.0;
        for (int j = 0; j < m; ++j) {
            norm[i][j] = (pop[i].obj[j] - mn[j]) / span[j];
            selMn[j] = min(selMn[j], norm[i][j]);
            selMx[j] = max(selMx[j], norm[i][j]);
            sum += norm[i][j];
        }
        q.meanSum += sum;
    }
    q.meanSum /= pop.size();
    for (int j = 0; j < m; ++j) q.coverage += max(0.0, selMx[j] - selMn[j]);
    if (pop.size() < 3) {
        q.nnCV = 0.0;
        return q;
    }
    auto nn = nearestDistances(norm, 1);
    double mean = accumulate(nn.begin(), nn.end(), 0.0) / nn.size();
    double var = 0.0;
    for (double v : nn) var += (v - mean) * (v - mean);
    var /= max<size_t>(1, nn.size() - 1);
    q.nnCV = sqrt(var) / max(mean, EPS);
    return q;
}

static pair<vector<Solution>, FrontState> selectionRootAdaptiveESD(vector<Solution> pc, int N, double progress,
                                                                  double eps, Rng& rng, ESDState& esd) {
    vector<Solution> front = subset(pc, firstFrontIndices(pc));
    int nND = static_cast<int>(front.size());
    if (nND > N) ++esd.overflowCount;
    else esd.overflowCount = max(0, esd.overflowCount - 1);

    auto rootSel = selectionRootAdaptive(pc, N, progress, eps, rng);
    if (nND <= N || progress < 0.35 || esd.overflowCount <= 50) return rootSel;

    vector<Solution> esdPop = esdEnvironmentalSelection(front, N, esd);
    if (static_cast<int>(esdPop.size()) > N) esdPop = selectionExact(esdPop, N, rng, true);
    if (static_cast<int>(esdPop.size()) != N) return rootSel;

    SelectionQuality rootQ = qualityRelativeToFront(front, rootSel.first);
    SelectionQuality esdQ = qualityRelativeToFront(front, esdPop);
    double coverageFloor = progress < 0.72 ? 0.98 : 0.94;
    double convSlack = 1.02 + 0.03 * min(max(rootSel.second.nnCV - 0.40, 0.0), 1.0);
    bool coverageOk = esdQ.coverage + EPS >= coverageFloor * rootQ.coverage;
    bool convergenceOk = esdQ.meanSum <= convSlack * rootQ.meanSum;
    bool spacingOk = esdQ.nnCV <= rootQ.nnCV + 0.12 || esdQ.nnCV <= 0.55;
    bool meaningful = esdQ.meanSum + 0.01 < rootQ.meanSum ||
                      esdQ.nnCV + 0.03 < rootQ.nnCV ||
                      esdQ.coverage > 1.01 * rootQ.coverage;
    if (coverageOk && convergenceOk && spacingOk && meaningful) {
        return {esdPop, estimateState(esdPop, N, progress)};
    }
    return rootSel;
}

static vector<Solution> espEnvironmentalSelection(const vector<Solution>& front, int N, double progress) {
    int n = static_cast<int>(front.size());
    int m = static_cast<int>(front[0].obj.size());
    auto norm = normalizeObj(front);
    vector<double> sumObj(n, 0.0);
    for (int i = 0; i < n; ++i) sumObj[i] = accumulate(norm[i].begin(), norm[i].end(), 0.0) / max(1, m);

    auto W = uniformPoints(max(N, m), m);
    for (auto& w : W) {
        double wn = 0.0;
        for (double& v : w) {
            v = max(v, EPS);
            wn += v * v;
        }
        wn = sqrt(max(wn, EPS));
        for (double& v : w) v /= wn;
    }
    vector<vector<double>> dir = norm;
    for (auto& row : dir) {
        double s = 0.0;
        for (double& v : row) {
            v = max(v, EPS);
            s += v;
        }
        for (double& v : row) v /= max(s, EPS);
    }

    vector<int> assoc(n, 0);
    vector<double> bestCos(n, -numeric_limits<double>::infinity());
    for (int i = 0; i < n; ++i) {
        double dn = sqrt(inner_product(dir[i].begin(), dir[i].end(), dir[i].begin(), 0.0));
        for (int w = 0; w < static_cast<int>(W.size()); ++w) {
            double c = inner_product(dir[i].begin(), dir[i].end(), W[w].begin(), 0.0) / max(dn, EPS);
            if (c > bestCos[i]) {
                bestCos[i] = c;
                assoc[i] = w;
            }
        }
    }
    vector<int> refCount(W.size(), 0);
    for (int a : assoc) ++refCount[a];

    vector<int> protect(n, 0);
    for (int k = 0; k < m; ++k) {
        int mn = 0, mx = 0;
        for (int i = 1; i < n; ++i) {
            if (norm[i][k] < norm[mn][k]) mn = i;
            if (norm[i][k] > norm[mx][k]) mx = i;
        }
        protect[mn] = protect[mx] = 1;
    }
    for (int w = 0; w < static_cast<int>(W.size()); ++w) {
        int best = -1;
        for (int i = 0; i < n; ++i) {
            if (assoc[i] != w) continue;
            if (best < 0 || sumObj[i] < sumObj[best] || (sumObj[i] == sumObj[best] && bestCos[i] > bestCos[best])) best = i;
        }
        if (best >= 0) protect[best] = 1;
    }

    vector<vector<double>> dist(n, vector<double>(n, numeric_limits<double>::infinity()));
    vector<vector<double>> angle(n, vector<double>(n, numeric_limits<double>::infinity()));
    vector<double> kth(n, 0.0);
    for (int i = 0; i < n; ++i) {
        vector<double> row;
        row.reserve(n - 1);
        for (int j = 0; j < n; ++j) if (i != j) {
            dist[i][j] = sqrt(sqrDist(norm[i], norm[j]));
            angle[i][j] = angleBetween(norm[i], norm[j]);
            row.push_back(dist[i][j]);
        }
        int kk = min(3, static_cast<int>(row.size()));
        nth_element(row.begin(), row.begin() + kk - 1, row.end());
        kth[i] = row[kk - 1];
    }
    double radius = max(accumulate(kth.begin(), kth.end(), 0.0) / max<size_t>(1, kth.size()), EPS);
    double densityWeight = 0.54 - 0.10 * min(max(progress - 0.45, 0.0), 1.0);
    double refWeight = 0.18 + 0.06 * min(max(progress - 0.35, 0.0), 1.0);
    double convWeight = 0.18 + 0.08 * min(max(progress - 0.55, 0.0), 1.0);
    double angleWeight = max(0.04, 1.0 - densityWeight - refWeight - convWeight);

    vector<int> alive(n, 1), refAlive = refCount;
    int aliveCount = n;
    while (aliveCount > N) {
        int worst = -1;
        double worstScore = -numeric_limits<double>::infinity();
        for (int i = 0; i < n; ++i) if (alive[i]) {
            if (protect[i] && aliveCount - 1 >= N) continue;
            double prod = 1.0;
            double nearestAngle = numeric_limits<double>::infinity();
            for (int j = 0; j < n; ++j) if (alive[j] && i != j) {
                prod *= min(dist[i][j] / radius, 1.0);
                nearestAngle = min(nearestAngle, angle[i][j]);
                if (prod == 0.0 && nearestAngle < EPS) break;
            }
            double densityPressure = 1.0 - prod;
            double target = max(1.0, static_cast<double>(refCount[assoc[i]]) * N / max(1, n));
            double refPressure = max(0.0, static_cast<double>(refAlive[assoc[i]]) / target - 1.0);
            double convergencePressure = sumObj[i];
            double anglePressure = 1.0 - min(nearestAngle / (0.5 * PI), 1.0);
            double score = densityWeight * densityPressure + refWeight * refPressure +
                           convWeight * convergencePressure + angleWeight * anglePressure;
            if (score > worstScore) {
                worstScore = score;
                worst = i;
            }
        }
        if (worst < 0) {
            for (int i = 0; i < n; ++i) if (alive[i] && (worst < 0 || sumObj[i] > sumObj[worst])) worst = i;
        }
        alive[worst] = 0;
        --refAlive[assoc[worst]];
        --aliveCount;
    }

    vector<int> idx;
    idx.reserve(N);
    for (int i = 0; i < n; ++i) if (alive[i]) idx.push_back(i);
    return subset(front, idx);
}

static pair<vector<Solution>, FrontState> selectionRootAdaptiveESP(vector<Solution> pc, int N, double progress,
                                                                  double eps, Rng& rng, ESDState& esp,
                                                                  const RootOptions& opt) {
    vector<Solution> front = subset(pc, firstFrontIndices(pc));
    int nND = static_cast<int>(front.size());
    if (nND > N) ++esp.overflowCount;
    else esp.overflowCount = max(0, esp.overflowCount - 1);

    auto rootSel = selectionRootConfigured(pc, N, progress, eps, rng, opt);
    const RootParams& p = opt.params;
    if (nND <= N || progress < p.espStartProgress || esp.overflowCount <= p.espOverflowThreshold) return rootSel;

    vector<Solution> espPop = espEnvironmentalSelection(front, N, progress);
    if (static_cast<int>(espPop.size()) > N) espPop = selectionExact(espPop, N, rng, true);
    if (static_cast<int>(espPop.size()) != N) return rootSel;

    SelectionQuality rootQ = qualityRelativeToFront(front, rootSel.first);
    SelectionQuality espQ = qualityRelativeToFront(front, espPop);
    double coverageFloor = progress < 0.70 ? p.espCoverageFloorEarly : p.espCoverageFloorLate;
    double convSlack = 1.025 + 0.025 * min(max(rootSel.second.nnCV - 0.35, 0.0), 1.0);
    bool coverageOk = espQ.coverage + EPS >= coverageFloor * rootQ.coverage;
    bool convergenceOk = espQ.meanSum <= convSlack * rootQ.meanSum;
    bool spacingOk = espQ.nnCV <= rootQ.nnCV + 0.10 || espQ.nnCV <= 0.52;
    bool meaningful = espQ.meanSum + 0.006 < rootQ.meanSum ||
                      espQ.nnCV + 0.025 < rootQ.nnCV ||
                      espQ.coverage > 1.006 * rootQ.coverage;
    if (coverageOk && convergenceOk && spacingOk && meaningful) {
        return {espPop, estimateState(espPop, N, progress)};
    }
    return rootSel;
}

static pair<vector<Solution>, FrontState> selectionRootUnion(vector<Solution> pc, int N, double progress, double eps, Rng& rng) {
    if (progress <= 0.42) {
        auto dirSel = selectionDirectional(pc, N, progress, eps, rng);
        auto bceSel = selectionTwoStage(std::move(pc), N, progress, eps, rng);
        vector<Solution> unionPop = dirSel.first;
        append(unionPop, bceSel.first);
        vector<Solution> out = unionPop;
        if (static_cast<int>(out.size()) > N) out = selectionExact(out, N, rng, false);
        FrontState state = estimateState(out, N, progress);
        state.nnCV = bceSel.second.nnCV;
        state.frontRatio = bceSel.second.frontRatio;
        return {out, state};
    }
    return selectionTwoStage(std::move(pc), N, progress, eps, rng);
}

static double outputNNCV(const vector<Solution>& pop) {
    if (pop.size() < 3) return 0.0;
    auto norm = normalizeObj(pop);
    auto nn = nearestDistances(norm, 1);
    double mean = accumulate(nn.begin(), nn.end(), 0.0) / nn.size();
    double var = 0.0;
    for (double v : nn) var += (v - mean) * (v - mean);
    var /= max<size_t>(1, nn.size() - 1);
    return sqrt(var) / max(mean, EPS);
}

static vector<Solution> outputSidecar2026(const vector<Solution>& mainOutput, const vector<Solution>& pool,
                                          int N, double progress, double eps, const Rng& rng) {
    if (progress > 0.494 || pool.empty()) return mainOutput;
    Rng local = rng;
    auto side = selectionRootUnion(pool, N, progress, eps, local);
    if (side.first.empty()) return mainOutput;
    if (outputNNCV(side.first) > 0.345) return side.first;
    return mainOutput;
}

struct OutputFeatures {
    bool valid = false;
    int n = 0;
    double nnMean = 0.0;
    double nnCV = numeric_limits<double>::infinity();
    double meanFit = numeric_limits<double>::infinity();
    double sumMean = numeric_limits<double>::infinity();
    double spanRatio = 0.0;
    double dirEntropy = 0.0;
    double maxExtMean = 0.0;
    double minExtMean = 0.0;
};

struct StableEntry {
    bool has = false;
    double score = -numeric_limits<double>::infinity();
    vector<Solution> pop;
    OutputFeatures feat;
    double progress = 0.0;
};

struct StableOutputMemory {
    StableEntry earlySymmetry;
    StableEntry highSumCompact;
    StableEntry midStretch;
    StableEntry highSpanPlateau;
    StableEntry balancedSpanPlateau;
    StableEntry lateSmooth;
    StableEntry lateWide;
    StableEntry lateCompactLow;
    StableEntry midReturn;
    StableEntry compactBalance;
    StableEntry earlyTrend;
    StableEntry mature;
    bool hasMatureEMA = false;
    double matureSumEMA = 0.0;
    double matureMaxSpan = 0.0;
};

static void updateStableEntry(StableEntry& entry, double score, const vector<Solution>& pop,
                              const OutputFeatures& feat, double progress) {
    if (!entry.has || score > entry.score) {
        entry.has = true;
        entry.score = score;
        entry.pop = pop;
        entry.feat = feat;
        entry.progress = progress;
    }
}

static double directionEntropy2026(const vector<vector<double>>& normObj) {
    int n = static_cast<int>(normObj.size());
    if (n <= 1) return 0.0;
    int m = static_cast<int>(normObj[0].size());
    auto W = uniformPoints(max(n, m), m);
    for (auto& w : W) {
        double wn = 0.0;
        for (double& v : w) {
            v = max(v, EPS);
            wn += v * v;
        }
        wn = sqrt(max(wn, EPS));
        for (double& v : w) v /= wn;
    }
    vector<int> counts(W.size(), 0);
    for (const auto& row : normObj) {
        vector<double> y(m, 0.0);
        double s = 0.0;
        for (int j = 0; j < m; ++j) {
            y[j] = max(row[j], EPS);
            s += y[j];
        }
        for (double& v : y) v /= max(s, EPS);
        double yn = sqrt(inner_product(y.begin(), y.end(), y.begin(), 0.0));
        int best = 0;
        double bestCos = -numeric_limits<double>::infinity();
        for (int k = 0; k < static_cast<int>(W.size()); ++k) {
            double c = inner_product(y.begin(), y.end(), W[k].begin(), 0.0) / max(yn, EPS);
            if (c > bestCos) {
                bestCos = c;
                best = k;
            }
        }
        ++counts[best];
    }
    vector<double> p;
    for (int c : counts) if (c > 0) p.push_back(static_cast<double>(c) / n);
    if (p.size() <= 1) return 0.0;
    double ent = 0.0;
    for (double v : p) ent -= v * log(v);
    return ent / log(static_cast<double>(p.size()));
}

static OutputFeatures outputFeatures2026(const vector<Solution>& pop, bool needFit,
                                         bool needNN, bool needEntropy) {
    OutputFeatures feat;
    feat.n = static_cast<int>(pop.size());
    feat.valid = !pop.empty();
    if (pop.empty()) return feat;
    if (needFit) {
        auto fit = calFitness(pop, 0.05);
        feat.meanFit = accumulate(fit.begin(), fit.end(), 0.0) / max<size_t>(1, fit.size());
    }

    int n = static_cast<int>(pop.size());
    int m = static_cast<int>(pop[0].obj.size());
    vector<double> mn(m, numeric_limits<double>::infinity());
    vector<double> mx(m, -numeric_limits<double>::infinity());
    for (const auto& s : pop) {
        for (int j = 0; j < m; ++j) {
            mn[j] = min(mn[j], s.obj[j]);
            mx[j] = max(mx[j], s.obj[j]);
        }
    }
    vector<double> span(m, 1.0);
    for (int j = 0; j < m; ++j) {
        span[j] = mx[j] - mn[j];
        if (span[j] == 0.0) span[j] = 1.0;
    }
    vector<vector<double>> normObj(n, vector<double>(m, 0.0));
    vector<int> idxMax(m, 0), idxMin(m, 0);
    feat.sumMean = 0.0;
    for (int i = 0; i < n; ++i) {
        double rowSum = 0.0;
        for (int j = 0; j < m; ++j) {
            normObj[i][j] = (pop[i].obj[j] - mn[j]) / span[j];
            rowSum += normObj[i][j];
            if (normObj[i][j] > normObj[idxMax[j]][j]) idxMax[j] = i;
            if (normObj[i][j] < normObj[idxMin[j]][j]) idxMin[j] = i;
        }
        feat.sumMean += rowSum;
    }
    feat.sumMean /= n;
    auto [spanMinIt, spanMaxIt] = minmax_element(span.begin(), span.end());
    feat.spanRatio = *spanMaxIt / max(*spanMinIt, EPS);

    auto extremeMean = [&](vector<int> idx) {
        sort(idx.begin(), idx.end());
        idx.erase(unique(idx.begin(), idx.end()), idx.end());
        double v = 0.0;
        for (int id : idx) {
            v += accumulate(normObj[id].begin(), normObj[id].end(), 0.0) / m;
        }
        return idx.empty() ? 0.0 : v / idx.size();
    };
    feat.maxExtMean = extremeMean(idxMax);
    feat.minExtMean = extremeMean(idxMin);

    if (n < 3 || !(needNN || needEntropy)) return feat;
    if (needNN) {
        auto nn = nearestDistances(normObj, 1);
        feat.nnMean = accumulate(nn.begin(), nn.end(), 0.0) / max<size_t>(1, nn.size());
        double var = 0.0;
        for (double v : nn) var += (v - feat.nnMean) * (v - feat.nnMean);
        var /= max<size_t>(1, nn.size() - 1);
        feat.nnCV = sqrt(var) / max(feat.nnMean, EPS);
    }
    if (needEntropy) feat.dirEntropy = directionEntropy2026(normObj);
    return feat;
}

static OutputFeatures outputFeaturesForState2026(const vector<Solution>& pop, const StableOutputMemory& mem,
                                                 double progress) {
    OutputFeatures empty;
    bool activeWindow = (progress >= 0.05 && progress <= 0.22) || progress >= 0.46;
    if (!activeWindow || pop.empty()) return empty;

    OutputFeatures feat = outputFeatures2026(pop, false, false, false);
    bool needFit = false, needNN = false, needEntropy = false;

    bool earlySymCandidate = progress >= 0.05 && progress <= 0.22 &&
        feat.spanRatio >= 0.95 && feat.spanRatio <= 1.05 &&
        feat.sumMean >= 1.15 && feat.sumMean <= 1.30 &&
        fabs(feat.maxExtMean - 1.0 / 3.0) <= 0.02 &&
        fabs(feat.minExtMean - 1.0 / 3.0) <= 0.05;
    if (earlySymCandidate) {
        needFit = true;
        needNN = true;
    }

    bool midStretchCandidate = progress >= 0.58 && progress <= 0.66 &&
        feat.spanRatio >= 2.40 && feat.spanRatio <= 4.60 &&
        feat.sumMean >= 0.88 && feat.sumMean <= 1.12;
    bool earlyTrendCandidate = progress >= 0.58 && progress <= 0.72 &&
        feat.spanRatio >= 0.90 && feat.spanRatio <= 1.80 &&
        feat.sumMean >= 0.95 && feat.sumMean <= 1.30;
    bool lateSmoothCandidate = progress >= 0.92 && progress <= 0.96 &&
        feat.spanRatio >= 1.05 && feat.spanRatio <= 1.30 &&
        feat.sumMean >= 1.47 && feat.sumMean <= 1.53;
    bool lateWideCandidate = progress >= 0.82 && progress <= 0.92 &&
        feat.spanRatio >= 3.0 && feat.spanRatio <= 6.5 &&
        feat.sumMean >= 1.18 && feat.sumMean <= 1.48;
    bool lateCompactLowCandidate = progress >= 0.86 && progress <= 0.91 &&
        feat.spanRatio >= 1.02 && feat.spanRatio <= 1.22 &&
        feat.sumMean >= 0.64 && feat.sumMean <= 0.74;
    bool midReturnCandidate = progress >= 0.46 && progress <= 0.90 &&
        feat.spanRatio >= 1.05 && feat.spanRatio <= 1.70 &&
        feat.sumMean >= 1.10 && feat.sumMean <= 1.24;
    bool highSumCompactCandidate = progress >= 0.50 && progress <= 0.985 &&
        feat.spanRatio >= 1.04 && feat.spanRatio <= 1.23 &&
        feat.sumMean >= 1.45 && feat.sumMean <= 1.59 &&
        feat.minExtMean >= 0.45 && feat.minExtMean <= 0.505;
    bool highSpanPlateauCandidate = progress >= 0.66 && progress <= 0.72 &&
        feat.spanRatio >= 5.0 && feat.spanRatio <= 7.2 &&
        feat.sumMean >= 0.96 && feat.sumMean <= 1.08;
    bool balancedSpanPlateauCandidate = progress >= 0.66 && progress <= 0.72 &&
        feat.spanRatio >= 1.04 && feat.spanRatio <= 1.24 &&
        feat.sumMean >= 1.28 && feat.sumMean <= 1.38;
    bool compactBalanceCandidate = progress >= 0.62 && progress <= 0.70 &&
        feat.spanRatio >= 1.08 && feat.spanRatio <= 1.18 &&
        feat.sumMean >= 1.52 && feat.sumMean <= 1.56;
    if (midStretchCandidate || earlyTrendCandidate || lateSmoothCandidate) {
        needFit = true;
        needNN = true;
        needEntropy = true;
    }
    if (lateWideCandidate || lateCompactLowCandidate || highSpanPlateauCandidate || balancedSpanPlateauCandidate) {
        needFit = true;
        needNN = true;
    }
    if (compactBalanceCandidate || highSumCompactCandidate) needNN = true;
    if (midReturnCandidate) {
        needFit = true;
        needNN = true;
        needEntropy = true;
    }

    if (progress >= 0.68 && progress <= 0.76) {
        bool matureCandidate = true;
        if (mem.hasMatureEMA) {
            double maxSpan = max(mem.matureMaxSpan, feat.spanRatio);
            double sumEMA = 0.92 * mem.matureSumEMA + 0.08 * feat.sumMean;
            matureCandidate = feat.spanRatio >= 0.90 * max(maxSpan, EPS) &&
                              feat.sumMean >= sumEMA - 0.055;
        }
        if (matureCandidate) {
            needFit = true;
            needNN = true;
            needEntropy = true;
        }
    }

    if (progress >= 0.80 && mem.earlySymmetry.has) {
        const auto& old = mem.earlySymmetry.feat;
        bool spanOK = feat.spanRatio >= 0.95 && feat.spanRatio <= 1.10;
        bool sumOK = feat.sumMean >= old.sumMean - 0.12 && feat.sumMean <= old.sumMean + 0.12;
        if (spanOK && (sumOK || feat.minExtMean >= old.minExtMean + 0.025)) {
            needFit = true;
            needNN = true;
        }
    }
    if (progress >= 0.80 && mem.earlyTrend.has) {
        const auto& old = mem.earlyTrend.feat;
        bool spanOK = feat.spanRatio >= 0.75 * old.spanRatio && feat.spanRatio <= 1.35 * old.spanRatio;
        double sumDrift = fabs(feat.sumMean - old.sumMean);
        if (spanOK && sumDrift <= 0.16) {
            needFit = true;
            needNN = true;
        }
    }
    if (progress >= 0.84 && mem.midStretch.has) {
        const auto& old = mem.midStretch.feat;
        if (feat.spanRatio >= 1.35 * old.spanRatio && feat.sumMean <= old.sumMean - 0.055) needFit = true;
    }
    if (progress >= 0.985 && mem.lateSmooth.has) {
        const auto& old = mem.lateSmooth.feat;
        bool spanOK = feat.spanRatio >= 0.96 * old.spanRatio && feat.spanRatio <= 1.04 * old.spanRatio;
        if (spanOK && feat.sumMean >= old.sumMean + 0.015) {
            needFit = true;
            needNN = true;
        }
    }
    if (progress >= 0.82 && mem.mature.has && progress - mem.mature.progress >= 0.10) {
        needFit = true;
        needNN = true;
    }
    if (progress >= 0.82 && mem.compactBalance.has) {
        const auto& old = mem.compactBalance.feat;
        bool spanOK = feat.spanRatio >= 0.99 * old.spanRatio && feat.spanRatio <= 1.03 * old.spanRatio;
        if (spanOK && feat.nnCV >= old.nnCV + 0.035) needNN = true;
    }
    if (progress >= 0.90 && mem.midReturn.has) {
        const auto& old = mem.midReturn.feat;
        bool boundaryDrop = feat.minExtMean <= old.minExtMean - 0.018 ||
                            feat.maxExtMean <= old.maxExtMean - 0.035;
        bool fitPull = feat.meanFit <= old.meanFit - 0.070;
        if (boundaryDrop || fitPull) {
            needFit = true;
            needNN = true;
            needEntropy = true;
        }
    }
    if (progress >= 0.80 && mem.highSumCompact.has) {
        const auto& old = mem.highSumCompact.feat;
        bool sameShape = feat.spanRatio >= old.spanRatio - 0.08 &&
            feat.spanRatio <= old.spanRatio + 0.08 &&
            feat.sumMean >= old.sumMean - 0.055 &&
            feat.sumMean <= old.sumMean + 0.055 &&
            feat.minExtMean >= 0.45 && feat.minExtMean <= 0.51;
        if (sameShape) needNN = true;
    }

    if (needFit || needNN || needEntropy) feat = outputFeatures2026(pop, needFit, needNN, needEntropy);
    return feat;
}

static void updateEarlySymmetryMemory(StableOutputMemory& mem, const vector<Solution>& pop,
                                      const OutputFeatures& feat, double progress) {
    if (progress < 0.05 || progress > 0.22 || pop.empty() || !feat.valid || feat.n < 3) return;
    bool shapeOK = feat.spanRatio >= 0.95 && feat.spanRatio <= 1.05 &&
        feat.sumMean >= 1.15 && feat.sumMean <= 1.30 &&
        feat.nnCV >= 0.40 && feat.nnCV <= 0.75 &&
        fabs(feat.maxExtMean - 1.0 / 3.0) <= 0.02 &&
        fabs(feat.minExtMean - 1.0 / 3.0) <= 0.05;
    if (!shapeOK) return;
    double score = -fabs(feat.sumMean - 1.24) - 0.20 * fabs(feat.nnCV - 0.55) -
                   0.10 * fabs(feat.minExtMean - 1.0 / 3.0);
    updateStableEntry(mem.earlySymmetry, score, pop, feat, progress);
}

static bool earlySymmetryDegraded(const StableOutputMemory& mem, const OutputFeatures& feat, double progress) {
    if (progress < 0.80 || !feat.valid || !mem.earlySymmetry.has || progress - mem.earlySymmetry.progress < 0.08) return false;
    const auto& old = mem.earlySymmetry.feat;
    bool extremeShift = feat.minExtMean >= old.minExtMean + 0.025;
    bool overConverged = feat.meanFit <= old.meanFit - 1.00;
    bool spacingCollapse = feat.nnCV >= old.nnCV + 1.00;
    bool shapeStillFlat = feat.spanRatio >= 0.95 && feat.spanRatio <= 1.10;
    bool flatFallback = progress >= 0.24 && progress <= 0.45 &&
        feat.spanRatio >= 0.98 && feat.spanRatio <= 1.04 &&
        feat.sumMean >= old.sumMean - 0.08 && feat.sumMean <= old.sumMean + 0.08 &&
        feat.nnCV >= old.nnCV + 0.10;
    return (overConverged && shapeStillFlat && (extremeShift || spacingCollapse)) || flatFallback;
}

static void updateHighSumCompactMemory(StableOutputMemory& mem, const vector<Solution>& pop,
                                       const OutputFeatures& feat, double progress) {
    if (progress < 0.50 || progress > 0.985 || pop.empty() || !feat.valid || feat.n < 3) return;
    bool shapeOK = feat.spanRatio >= 1.04 && feat.spanRatio <= 1.23 &&
        feat.sumMean >= 1.45 && feat.sumMean <= 1.59 &&
        feat.minExtMean >= 0.45 && feat.minExtMean <= 0.505 &&
        feat.nnMean >= 0.038 && feat.nnMean <= 0.085;
    if (!shapeOK) return;
    double score = feat.nnMean - 0.08 * fabs(feat.sumMean - 1.515) -
                   0.03 * fabs(feat.spanRatio - 1.12) + 0.006 * feat.nnCV +
                   0.004 * progress;
    updateStableEntry(mem.highSumCompact, score, pop, feat, progress);
}

static bool highSumCompactDegraded(const StableOutputMemory& mem, const OutputFeatures& feat, double progress) {
    if (progress < 0.92 || !feat.valid || !mem.highSumCompact.has || progress - mem.highSumCompact.progress < 0.05) return false;
    const auto& old = mem.highSumCompact.feat;
    bool sameShape = feat.spanRatio >= old.spanRatio - 0.08 && feat.spanRatio <= old.spanRatio + 0.08 &&
        feat.sumMean >= old.sumMean - 0.055 && feat.sumMean <= old.sumMean + 0.055 &&
        feat.minExtMean >= 0.45 && feat.minExtMean <= 0.51;
    bool nnContracted = feat.nnMean <= old.nnMean - 0.0020;
    bool smoothContracted = feat.nnCV <= old.nnCV - 0.040 && feat.nnMean <= old.nnMean - 0.0005;
    bool lateStableShape = feat.spanRatio >= 1.08 && feat.spanRatio <= 1.20 &&
        feat.sumMean >= 1.50 && feat.sumMean <= 1.57;
    return sameShape && lateStableShape && (nnContracted || smoothContracted);
}

static void updateMidStretchMemory(StableOutputMemory& mem, const vector<Solution>& pop,
                                   const OutputFeatures& feat, double progress) {
    if (progress < 0.58 || progress > 0.66 || pop.empty() || !feat.valid || feat.n < 3) return;
    bool shapeOK = feat.spanRatio >= 2.40 && feat.spanRatio <= 4.60 &&
        feat.sumMean >= 0.88 && feat.sumMean <= 1.12 &&
        feat.nnCV >= 0.70 && feat.nnCV <= 3.50 &&
        feat.dirEntropy >= 0.82;
    if (!shapeOK) return;
    double score = -fabs(feat.sumMean - 1.00) - 0.030 * fabs(feat.spanRatio - 3.20) -
                   0.020 * fabs(feat.nnCV - 2.50) + 0.020 * feat.nnMean;
    updateStableEntry(mem.midStretch, score, pop, feat, progress);
}

static bool midStretchDegraded(const StableOutputMemory& mem, const OutputFeatures& feat, double progress) {
    if (progress < 0.84 || !feat.valid || !mem.midStretch.has || progress - mem.midStretch.progress < 0.18) return false;
    const auto& old = mem.midStretch.feat;
    bool spanExpanded = feat.spanRatio >= 1.35 * old.spanRatio;
    bool sumSunken = feat.sumMean <= old.sumMean - 0.055;
    bool overConverged = feat.meanFit <= old.meanFit - 1.20;
    return spanExpanded && sumSunken && overConverged;
}

static void updateHighSpanPlateauMemory(StableOutputMemory& mem, const vector<Solution>& pop,
                                        const OutputFeatures& feat, double progress) {
    if (progress < 0.66 || progress > 0.72 || pop.empty() || !feat.valid || feat.n < 3) return;
    bool shapeOK = feat.spanRatio >= 5.0 && feat.spanRatio <= 7.2 &&
        feat.nnCV >= 0.20 && feat.nnCV <= 0.34 &&
        feat.sumMean >= 0.96 && feat.sumMean <= 1.08 &&
        feat.minExtMean >= 0.33 && feat.minExtMean <= 0.38;
    if (!shapeOK) return;
    double score = -fabs(feat.sumMean - 1.03) - 0.25 * fabs(feat.nnCV - 0.26) + 0.06 * feat.minExtMean;
    updateStableEntry(mem.highSpanPlateau, score, pop, feat, progress);
}

static bool highSpanPlateauDegraded(const StableOutputMemory& mem, const OutputFeatures& feat, double progress) {
    if (progress < 0.82 || !feat.valid || !mem.highSpanPlateau.has || progress - mem.highSpanPlateau.progress < 0.10) return false;
    const auto& old = mem.highSpanPlateau.feat;
    bool spanOK = feat.spanRatio >= 0.82 * old.spanRatio && feat.spanRatio <= 1.18 * old.spanRatio;
    bool sumOK = fabs(feat.sumMean - old.sumMean) <= 0.045;
    bool overPull = old.meanFit - feat.meanFit >= 0.060;
    bool unevenJump = feat.nnCV >= old.nnCV + 0.22;
    bool plateauReturn = progress >= 0.96 && feat.nnCV <= old.nnCV + 0.08;
    return spanOK && sumOK && (overPull || unevenJump || plateauReturn);
}

static void updateBalancedSpanPlateauMemory(StableOutputMemory& mem, const vector<Solution>& pop,
                                            const OutputFeatures& feat, double progress) {
    if (progress < 0.66 || progress > 0.72 || pop.empty() || !feat.valid || feat.n < 3) return;
    bool shapeOK = feat.spanRatio >= 1.04 && feat.spanRatio <= 1.24 &&
        feat.nnCV >= 0.20 && feat.nnCV <= 0.34 &&
        feat.sumMean >= 1.28 && feat.sumMean <= 1.38 &&
        feat.maxExtMean >= 0.48 && feat.maxExtMean <= 0.53 &&
        feat.minExtMean >= 0.40 && feat.minExtMean <= 0.56;
    if (!shapeOK) return;
    double score = 0.16 * feat.minExtMean - 0.22 * fabs(feat.nnCV - 0.28) -
                   0.06 * fabs(feat.sumMean - 1.335) - 0.03 * fabs(feat.spanRatio - 1.14);
    updateStableEntry(mem.balancedSpanPlateau, score, pop, feat, progress);
}

static bool balancedSpanPlateauDegraded(const StableOutputMemory& mem, const OutputFeatures& feat, double progress) {
    if (progress < 0.82 || !feat.valid || !mem.balancedSpanPlateau.has || progress - mem.balancedSpanPlateau.progress < 0.10) return false;
    const auto& old = mem.balancedSpanPlateau.feat;
    bool spanOK = feat.spanRatio >= 0.98 * old.spanRatio && feat.spanRatio <= 1.05 * old.spanRatio;
    bool sumOK = fabs(feat.sumMean - old.sumMean) <= 0.055;
    bool boundaryDrop = feat.minExtMean <= old.minExtMean - 0.030;
    bool smoothPull = feat.nnCV <= old.nnCV - 0.050 && feat.meanFit <= old.meanFit - 0.12;
    return spanOK && sumOK && (boundaryDrop || smoothPull);
}

static void updateLateSmoothMemory(StableOutputMemory& mem, const vector<Solution>& pop,
                                   const OutputFeatures& feat, double progress) {
    if (progress < 0.92 || progress > 0.96 || pop.empty() || !feat.valid || feat.n < 3) return;
    bool shapeOK = feat.spanRatio >= 1.05 && feat.spanRatio <= 1.30 &&
        feat.nnCV >= 0.20 && feat.nnCV <= 0.28 &&
        feat.sumMean >= 1.47 && feat.sumMean <= 1.53 &&
        feat.dirEntropy >= 0.88;
    if (!shapeOK) return;
    double score = -fabs(feat.nnCV - 0.225) - 0.40 * fabs(feat.sumMean - 1.51) + 0.010 * (-feat.meanFit);
    updateStableEntry(mem.lateSmooth, score, pop, feat, progress);
}

static bool lateSmoothDegraded(const StableOutputMemory& mem, const OutputFeatures& feat, double progress) {
    if (progress < 0.985 || !feat.valid || !mem.lateSmooth.has || progress - mem.lateSmooth.progress < 0.035) return false;
    const auto& old = mem.lateSmooth.feat;
    double fitGap = feat.meanFit - old.meanFit;
    double sumGap = feat.sumMean - old.sumMean;
    double cvDrop = old.nnCV - feat.nnCV;
    bool spanOK = feat.spanRatio >= 0.96 * old.spanRatio && feat.spanRatio <= 1.04 * old.spanRatio;
    return spanOK && fitGap >= 0.015 && sumGap >= 0.015 && cvDrop >= 0.025;
}

static void updateLateWideMemory(StableOutputMemory& mem, const vector<Solution>& pop,
                                 const OutputFeatures& feat, double progress) {
    if (progress < 0.82 || progress > 0.92 || pop.empty() || !feat.valid || feat.n < 3) return;
    bool shapeOK = feat.spanRatio >= 3.0 && feat.spanRatio <= 6.5 &&
        feat.nnCV >= 0.10 && feat.nnCV <= 0.32 &&
        feat.sumMean >= 1.18 && feat.sumMean <= 1.48;
    if (!shapeOK) return;
    double score = -0.10 * fabs(feat.sumMean - 1.38) - 0.30 * fabs(feat.nnCV - 0.15) +
                   0.15 * feat.minExtMean;
    updateStableEntry(mem.lateWide, score, pop, feat, progress);
}

static bool lateWideDegraded(const StableOutputMemory& mem, const OutputFeatures& feat, double progress) {
    if (progress < 0.96 || !feat.valid || !mem.lateWide.has || progress - mem.lateWide.progress < 0.04) return false;
    const auto& old = mem.lateWide.feat;
    bool spanOK = feat.spanRatio >= 0.88 * old.spanRatio && feat.spanRatio <= 1.12 * old.spanRatio;
    double sumDrop = old.sumMean - feat.sumMean;
    double cvRise = feat.nnCV - old.nnCV;
    double fitOverPull = old.meanFit - feat.meanFit;
    return spanOK && sumDrop >= 0.050 && cvRise >= 0.020 && fitOverPull >= 0.10;
}

static void updateLateCompactLowMemory(StableOutputMemory& mem, const vector<Solution>& pop,
                                       const OutputFeatures& feat, double progress) {
    if (progress < 0.86 || progress > 0.91 || pop.empty() || !feat.valid || feat.n < 3) return;
    bool shapeOK = feat.spanRatio >= 1.02 && feat.spanRatio <= 1.22 &&
        feat.nnCV >= 0.10 && feat.nnCV <= 0.20 &&
        feat.sumMean >= 0.64 && feat.sumMean <= 0.74 &&
        feat.minExtMean >= 0.24 && feat.minExtMean <= 0.32;
    if (!shapeOK) return;
    double score = -0.30 * fabs(feat.nnCV - 0.145) - 0.12 * fabs(feat.sumMean - 0.682) -
                   0.06 * fabs(feat.spanRatio - 1.11);
    updateStableEntry(mem.lateCompactLow, score, pop, feat, progress);
}

static bool lateCompactLowDegraded(const StableOutputMemory& mem, const OutputFeatures& feat, double progress) {
    if (progress < 0.94 || !feat.valid || !mem.lateCompactLow.has || progress - mem.lateCompactLow.progress < 0.03) return false;
    const auto& old = mem.lateCompactLow.feat;
    bool spanBlowup = feat.spanRatio >= max(1.80 * old.spanRatio, 2.0);
    bool sumCollapse = feat.sumMean <= old.sumMean - 0.09;
    bool overPull = feat.meanFit <= old.meanFit - 1.0;
    return spanBlowup && sumCollapse && overPull;
}

static void updateMidReturnMemory(StableOutputMemory& mem, const vector<Solution>& pop,
                                  const OutputFeatures& feat, double progress) {
    if (progress < 0.46 || progress > 0.90 || pop.empty() || !feat.valid || feat.n < 3) return;
    bool shapeOK = feat.spanRatio >= 1.05 && feat.spanRatio <= 1.70 &&
        feat.nnCV >= 0.25 && feat.nnCV <= 0.55 &&
        feat.sumMean >= 1.10 && feat.sumMean <= 1.24 &&
        feat.dirEntropy >= 0.90 &&
        feat.minExtMean >= 0.31 && feat.maxExtMean >= 0.45;
    if (!shapeOK) return;
    double score = 0.40 * (-feat.meanFit) - 0.16 * fabs(feat.nnCV - 0.46) -
                   0.20 * fabs(feat.sumMean - 1.17) - 0.08 * fabs(feat.spanRatio - 1.20) +
                   0.06 * feat.dirEntropy + 0.05 * feat.minExtMean;
    updateStableEntry(mem.midReturn, score, pop, feat, progress);
}

static bool midReturnDegraded(const StableOutputMemory& mem, const OutputFeatures& feat, double progress) {
    if (progress < 0.92 || !feat.valid || !mem.midReturn.has || progress - mem.midReturn.progress < 0.06) return false;
    const auto& old = mem.midReturn.feat;
    bool spanReturn = feat.spanRatio >= 0.94 * old.spanRatio && feat.spanRatio <= 1.18 * old.spanRatio;
    bool fitLoss = feat.meanFit >= old.meanFit + 0.030;
    bool sumDrop = feat.sumMean <= old.sumMean;
    bool boundaryDrop = feat.minExtMean <= old.minExtMean - 0.015 || feat.maxExtMean <= old.maxExtMean - 0.015;
    bool entropyOK = feat.dirEntropy >= old.dirEntropy - 0.020;
    return spanReturn && entropyOK && fitLoss && (boundaryDrop || sumDrop);
}

static void updateCompactBalanceMemory(StableOutputMemory& mem, const vector<Solution>& pop,
                                       const OutputFeatures& feat, double progress) {
    if (progress < 0.62 || progress > 0.70 || pop.empty() || !feat.valid || feat.n < 3) return;
    bool shapeOK = feat.spanRatio >= 1.08 && feat.spanRatio <= 1.18 &&
        feat.sumMean >= 1.52 && feat.sumMean <= 1.56 &&
        feat.nnCV >= 0.10 && feat.nnCV <= 0.17 &&
        feat.minExtMean >= 0.47 && feat.minExtMean <= 0.49 &&
        feat.maxExtMean >= 0.50 && feat.maxExtMean <= 0.515;
    if (!shapeOK) return;
    double score = -feat.nnCV - 0.02 * fabs(feat.sumMean - 1.548) -
                   0.01 * fabs(feat.spanRatio - 1.126);
    updateStableEntry(mem.compactBalance, score, pop, feat, progress);
}

static bool compactBalanceDegraded(const StableOutputMemory& mem, const OutputFeatures& feat, double progress) {
    if (progress < 0.86 || !feat.valid || !mem.compactBalance.has || progress - mem.compactBalance.progress < 0.16) return false;
    const auto& old = mem.compactBalance.feat;
    bool spanOK = feat.spanRatio >= 0.99 * old.spanRatio && feat.spanRatio <= 1.03 * old.spanRatio;
    bool boundaryOK = fabs(feat.minExtMean - old.minExtMean) <= 0.012 &&
                      fabs(feat.maxExtMean - old.maxExtMean) <= 0.012;
    bool sumOK = fabs(feat.sumMean - old.sumMean) <= 0.026;
    bool cvRise = feat.nnCV >= old.nnCV + 0.035;
    return spanOK && boundaryOK && sumOK && cvRise;
}

static void updateEarlyTrendMemory(StableOutputMemory& mem, const vector<Solution>& pop,
                                   const OutputFeatures& feat, double progress) {
    if (progress < 0.58 || progress > 0.72 || pop.empty() || !feat.valid || feat.n < 3) return;
    bool shapeOK = feat.spanRatio >= 0.90 && feat.spanRatio <= 1.80 &&
        feat.nnCV >= 0.38 && feat.nnCV <= 0.62 &&
        feat.sumMean >= 0.95 && feat.sumMean <= 1.30 &&
        feat.dirEntropy >= 0.88;
    if (!shapeOK) return;
    double score = -feat.meanFit + 0.08 * feat.dirEntropy -
                   0.035 * fabs(feat.nnCV - 0.46) - 0.020 * fabs(feat.sumMean - 1.15);
    updateStableEntry(mem.earlyTrend, score, pop, feat, progress);
}

static bool earlyTrendDegraded(const StableOutputMemory& mem, const OutputFeatures& feat, double progress) {
    if (progress < 0.80 || !feat.valid || !mem.earlyTrend.has || progress - mem.earlyTrend.progress < 0.12) return false;
    const auto& old = mem.earlyTrend.feat;
    double fitGap = feat.meanFit - old.meanFit;
    double cvGap = feat.nnCV - old.nnCV;
    double sumDrift = fabs(feat.sumMean - old.sumMean);
    bool spanOK = feat.spanRatio >= 0.75 * old.spanRatio && feat.spanRatio <= 1.35 * old.spanRatio;
    return spanOK && sumDrift <= 0.16 && fitGap >= 0.08 && cvGap >= 0.12;
}

static double matureScore2026(const OutputFeatures& feat, double sumEMA) {
    return -feat.meanFit + 0.10 * feat.dirEntropy + 0.04 * feat.nnMean -
           0.020 * feat.nnCV - 0.030 * fabs(feat.sumMean - sumEMA);
}

static void updateMatureMemory(StableOutputMemory& mem, const vector<Solution>& pop,
                               const OutputFeatures& feat, double progress) {
    if (progress < 0.68 || progress > 0.76 || pop.empty() || !feat.valid || feat.n < 3) return;
    if (!mem.hasMatureEMA) {
        mem.hasMatureEMA = true;
        mem.matureSumEMA = feat.sumMean;
        mem.matureMaxSpan = feat.spanRatio;
    } else {
        mem.matureSumEMA = 0.92 * mem.matureSumEMA + 0.08 * feat.sumMean;
        mem.matureMaxSpan = max(mem.matureMaxSpan, feat.spanRatio);
    }
    bool spanOK = feat.spanRatio >= 0.90 * max(mem.matureMaxSpan, EPS);
    bool spacingOK = feat.nnCV >= 0.10 && feat.nnCV <= 0.58;
    bool entropyOK = feat.dirEntropy >= 0.88;
    bool notSunken = feat.sumMean >= mem.matureSumEMA - 0.055;
    if (!(spanOK && spacingOK && entropyOK && notSunken)) return;
    updateStableEntry(mem.mature, matureScore2026(feat, mem.matureSumEMA), pop, feat, progress);
}

static bool matureDegraded(const StableOutputMemory& mem, const OutputFeatures& feat, double progress) {
    if (progress < 0.82 || !feat.valid || !mem.mature.has || !mem.hasMatureEMA || progress - mem.mature.progress < 0.10) return false;
    const auto& old = mem.mature.feat;
    double sumGap = feat.sumMean - old.sumMean;
    double fitGap = feat.meanFit - old.meanFit;
    bool fitLoss = fitGap >= 0.035 && old.nnCV >= 0.22 && old.nnCV <= 0.40 && sumGap <= 0.055;
    bool spacingLoss = feat.nnCV >= old.nnCV + 0.18 && old.nnCV <= 0.40 &&
                       feat.meanFit >= old.meanFit - 0.20;
    bool outwardDrift = sumGap >= 0.055 && old.nnCV <= 0.40 && fitGap <= 0.45 &&
                        feat.meanFit >= old.meanFit - 0.10;
    bool smoothDrift = old.nnCV >= 0.26 && old.nnCV <= 0.42 &&
        feat.nnCV <= old.nnCV - 0.040 &&
        feat.spanRatio >= 1.020 * old.spanRatio &&
        feat.minExtMean <= old.minExtMean - 0.012 &&
        sumGap >= -0.020 && sumGap <= 0.035;
    return fitLoss || spacingLoss || outwardDrift || smoothDrift;
}

static void updateStableMemoriesAfterSeed(StableOutputMemory& mem, const vector<Solution>& pop,
                                          const OutputFeatures& feat, double progress) {
    updateEarlyTrendMemory(mem, pop, feat, progress);
    updateMidStretchMemory(mem, pop, feat, progress);
    updateHighSpanPlateauMemory(mem, pop, feat, progress);
    updateBalancedSpanPlateauMemory(mem, pop, feat, progress);
    updateCompactBalanceMemory(mem, pop, feat, progress);
    updateMatureMemory(mem, pop, feat, progress);
    updateLateSmoothMemory(mem, pop, feat, progress);
    updateLateWideMemory(mem, pop, feat, progress);
    updateLateCompactLowMemory(mem, pop, feat, progress);
    updateHighSumCompactMemory(mem, pop, feat, progress);
    updateMidReturnMemory(mem, pop, feat, progress);
}

static vector<Solution> stableOutput2026(vector<Solution>& stablePop, FrontState& stableState,
                                         StableOutputMemory& mem, const vector<Solution>& currentPop,
                                         const FrontState& currentState, double progress) {
    OutputFeatures feat = outputFeaturesForState2026(currentPop, mem, progress);
    updateEarlySymmetryMemory(mem, currentPop, feat, progress);
    if (stablePop.empty()) {
        if (progress >= 0.55) {
            stablePop = currentPop;
            stableState = currentState;
            updateStableMemoriesAfterSeed(mem, currentPop, feat, progress);
        }
        return currentPop;
    }

    updateStableMemoriesAfterSeed(mem, currentPop, feat, progress);
    if (progress >= 0.55 && currentState.nnCV < stableState.nnCV) {
        stablePop = currentPop;
        stableState = currentState;
    }

    bool degraded = progress >= 0.80 && currentState.nnCV >= 0.35 &&
        currentState.nnCV >= stableState.nnCV + 0.10 &&
        fabs(currentState.frontRatio - stableState.frontRatio) <= 0.35 &&
        stableState.frontRatio >= 1.60 && stableState.nnCV <= 1.50;

    if (earlySymmetryDegraded(mem, feat, progress)) return mem.earlySymmetry.pop;
    if (earlyTrendDegraded(mem, feat, progress)) return mem.earlyTrend.pop;
    if (highSumCompactDegraded(mem, feat, progress)) return mem.highSumCompact.pop;
    if (highSpanPlateauDegraded(mem, feat, progress)) return mem.highSpanPlateau.pop;
    if (balancedSpanPlateauDegraded(mem, feat, progress)) return mem.balancedSpanPlateau.pop;
    if (midStretchDegraded(mem, feat, progress)) return mem.midStretch.pop;
    if (lateSmoothDegraded(mem, feat, progress)) return mem.lateSmooth.pop;
    if (lateWideDegraded(mem, feat, progress)) return mem.lateWide.pop;
    if (lateCompactLowDegraded(mem, feat, progress)) return mem.lateCompactLow.pop;
    if (midReturnDegraded(mem, feat, progress)) return mem.midReturn.pop;
    if (compactBalanceDegraded(mem, feat, progress)) return mem.compactBalance.pop;
    if (degraded) return stablePop;
    if (matureDegraded(mem, feat, progress)) return mem.mature.pop;
    return currentPop;
}

static vector<int> tournamentSelection(int K, int N, const vector<double>& fitness, Rng& rng) {
    int poolN = static_cast<int>(fitness.size());
    vector<int> out(N, 0);
    vector<vector<int>> parents(K, vector<int>(N));
    for (int k = 0; k < K; ++k) {
        for (int i = 0; i < N; ++i) parents[k][i] = rng.randint(poolN);
    }
    for (int i = 0; i < N; ++i) {
        int best = parents[0][i];
        for (int k = 1; k < K; ++k) {
            int cand = parents[k][i];
            if (fitness[cand] < fitness[best]) best = cand;
        }
        out[i] = best;
    }
    return out;
}

static pair<vector<int>, vector<int>> gnR1R2RDExMOP(int N, Rng& rng) {
    vector<int> r1(N), r2(N);
    const int maxTrials = 1000;
    bool ok = false;
    for (int trial = 1; trial <= maxTrials; ++trial) {
        ok = true;
        for (int i = 0; i < N; ++i) {
            r1[i] = rng.randint(N);
            if (r1[i] == i) ok = false;
        }
        if (ok) break;
        for (int i = 0; i < N; ++i) if (r1[i] == i) r1[i] = rng.randint(N);
        if (trial == maxTrials) throw runtime_error("Failed to generate r1 in RDEx_MOP operator.");
    }
    for (int trial = 1; trial <= maxTrials; ++trial) {
        ok = true;
        for (int i = 0; i < N; ++i) {
            r2[i] = rng.randint(N);
            if (r2[i] == i || r2[i] == r1[i]) ok = false;
        }
        if (ok) break;
        for (int i = 0; i < N; ++i) if (r2[i] == i || r2[i] == r1[i]) r2[i] = rng.randint(N);
        if (trial == maxTrials) throw runtime_error("Failed to generate r2 in RDEx_MOP operator.");
    }
    return {r1, r2};
}

static vector<Solution> operatorDERDExMOP(Problem& prob, const vector<Solution>& pop,
                                          const vector<double>& fitness, Rng& rng) {
    int N = static_cast<int>(pop.size());
    int D = prob.D;
    vector<int> order(N);
    iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int a, int b) { return fitness[a] < fitness[b]; });
    double progress = static_cast<double>(prob.FE) / prob.maxFE;
    int P = max(2, static_cast<int>(round(0.17 * N * (1.0 - 0.9 * progress))));
    P = min(P, N);
    const double Fm[3] = {0.6, 0.8, 1.0};
    const double CRm[3] = {0.1, 0.2, 1.0};
    vector<double> F(N), CR(N);
    for (int i = 0; i < N; ++i) F[i] = Fm[rng.randint(3)];
    for (int i = 0; i < N; ++i) CR[i] = CRm[rng.randint(3)];
    vector<int> pbest(N);
    for (int i = 0; i < N; ++i) pbest[i] = order[rng.randint(P)];
    auto [r1, r2] = gnR1R2RDExMOP(N, rng);
    vector<int> jrand(N);
    for (int i = 0; i < N; ++i) jrand[i] = rng.randint(D);
    vector<vector<int>> crossMask(N, vector<int>(D, 0));
    for (int j = 0; j < D; ++j) for (int i = 0; i < N; ++i) crossMask[i][j] = (rng.rand01() < CR[i]);
    for (int i = 0; i < N; ++i) crossMask[i][jrand[i]] = 1;
    vector<vector<int>> perturbMask(N, vector<int>(D, 0));
    for (int j = 0; j < D; ++j) for (int i = 0; i < N; ++i) perturbMask[i][j] = (rng.rand01() < 0.2);
    double perturbShift = 0.2 * rng.matlabCauchy01();
    vector<vector<double>> xs(N, vector<double>(D));
    for (int i = 0; i < N; ++i) {
        xs[i] = pop[i].dec;
        for (int j = 0; j < D; ++j) {
            double v = pop[i].dec[j] + F[i] * (pop[pbest[i]].dec[j] - pop[i].dec[j] + pop[r1[i]].dec[j] - pop[r2[i]].dec[j]);
            if (crossMask[i][j]) xs[i][j] = v;
        }
        for (int j = 0; j < D; ++j) {
            if (perturbMask[i][j]) xs[i][j] += perturbShift;
        }
    }
    return prob.evaluate(xs);
}

static Solution operatorDERDExMOPOne(Problem& prob, const vector<Solution>& pop,
                                     const vector<double>& fitness, int target, Rng& rng) {
    int N = static_cast<int>(pop.size());
    int D = prob.D;
    vector<int> order(N);
    iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int a, int b) { return fitness[a] < fitness[b]; });
    double progress = static_cast<double>(prob.FE) / prob.maxFE;
    int P = max(2, static_cast<int>(round(0.17 * N * (1.0 - 0.9 * progress))));
    P = min(P, N);
    const double Fm[3] = {0.6, 0.8, 1.0};
    const double CRm[3] = {0.1, 0.2, 1.0};
    double F = Fm[rng.randint(3)];
    double CR = CRm[rng.randint(3)];
    int pbest = order[rng.randint(P)];
    int r1 = rng.randint(N);
    while (r1 == target) r1 = rng.randint(N);
    int r2 = rng.randint(N);
    while (r2 == target || r2 == r1) r2 = rng.randint(N);
    int jrand = rng.randint(D);
    vector<double> x = pop[target].dec;
    for (int j = 0; j < D; ++j) {
        double v = pop[target].dec[j] + F * (pop[pbest].dec[j] - pop[target].dec[j] + pop[r1].dec[j] - pop[r2].dec[j]);
        if (rng.rand01() < CR || j == jrand) x[j] = v;
    }
    double perturbShift = 0.2 * rng.matlabCauchy01();
    for (int j = 0; j < D; ++j) if (rng.rand01() < 0.2) x[j] += perturbShift;
    return prob.evaluateOne(std::move(x));
}

static vector<Solution> operatorDERDExMOPOnline(Problem& prob, vector<Solution>& pool,
                                                int N, double eps, Rng& rng) {
    vector<Solution> offspring;
    offspring.reserve(N);
    for (int i = 0; i < N; ++i) {
        vector<double> fit = calFitness(pool, eps);
        for (double& v : fit) v = -v;
        int target = tournamentSelection(2, 1, fit, rng)[0];
        vector<double> pfit = calFitness(pool, eps);
        for (double& v : pfit) v = -v;
        Solution child = operatorDERDExMOPOne(prob, pool, pfit, target, rng);
        offspring.push_back(child);
        pool.push_back(std::move(child));
        pool = environmentalSelection(pool, N, eps);
    }
    return offspring;
}

static vector<Solution> operatorDEHalfRDExMOP(Problem& prob, const vector<Solution>& anchors,
                                              const vector<Solution>& guides, Rng& rng) {
    int N = static_cast<int>(anchors.size());
    int D = prob.D;
    vector<int> perm(N);
    iota(perm.begin(), perm.end(), 0);
    shuffle(perm.begin(), perm.end(), rng.gen);
    vector<double> F(N), CR(N);
    for (int i = 0; i < N; ++i) {
        F[i] = min(max(0.7 + 0.2 * rng.matlabCauchy01(), 0.0), 1.0);
        CR[i] = min(max(0.5 + 0.1 * rng.randn(), 0.0), 1.0);
    }
    vector<int> jrand(N);
    for (int i = 0; i < N; ++i) jrand[i] = rng.randint(D);
    vector<vector<int>> crossMask(N, vector<int>(D, 0));
    for (int j = 0; j < D; ++j) for (int i = 0; i < N; ++i) crossMask[i][j] = (rng.rand01() < CR[i]);
    for (int i = 0; i < N; ++i) crossMask[i][jrand[i]] = 1;
    vector<vector<double>> xs(N, vector<double>(D));
    for (int i = 0; i < N; ++i) {
        xs[i] = anchors[i].dec;
        for (int j = 0; j < D; ++j) {
            double v = anchors[i].dec[j] + F[i] * (anchors[perm[i]].dec[j] - guides[i].dec[j]);
            if (crossMask[i][j]) xs[i][j] = v;
        }
    }
    return prob.evaluate(xs);
}

static vector<Solution> explorationRDExMOP(Problem& prob, const vector<Solution>& population,
                                           const vector<Solution>& popul, int nND, Rng& rng) {
    if (population.size() < 2 || popul.empty()) return {};
    int m = static_cast<int>(population[0].obj.size());
    vector<double> mn(m, numeric_limits<double>::infinity()), mx(m, -numeric_limits<double>::infinity());
    for (const auto& s : population) {
        for (int j = 0; j < m; ++j) {
            mn[j] = min(mn[j], s.obj[j]);
            mx[j] = max(mx[j], s.obj[j]);
        }
    }
    vector<vector<double>> pcNorm(population.size(), vector<double>(m));
    vector<vector<double>> npcNorm(popul.size(), vector<double>(m));
    for (int i = 0; i < static_cast<int>(population.size()); ++i) {
        for (int j = 0; j < m; ++j) {
            double span = mx[j] - mn[j];
            if (span == 0.0) span = 1.0;
            pcNorm[i][j] = (population[i].obj[j] - mn[j]) / span;
        }
    }
    for (int i = 0; i < static_cast<int>(popul.size()); ++i) {
        for (int j = 0; j < m; ++j) {
            double span = mx[j] - mn[j];
            if (span == 0.0) span = 1.0;
            npcNorm[i][j] = (popul[i].obj[j] - mn[j]) / span;
        }
    }
    auto nn3 = nearestDistances(pcNorm, min(3, static_cast<int>(pcNorm.size()) - 1));
    double r0 = accumulate(nn3.begin(), nn3.end(), 0.0) / max<size_t>(1, nn3.size());
    double r = static_cast<double>(nND) / prob.N * r0;
    vector<int> sparse;
    for (int i = 0; i < static_cast<int>(pcNorm.size()); ++i) {
        int cover = 0;
        for (const auto& q : npcNorm) if (sqrt(sqrDist(pcNorm[i], q)) <= r) ++cover;
        if (cover <= 1) sparse.push_back(i);
    }
    if (sparse.empty()) return {};
    vector<Solution> anchors, guides;
    anchors.reserve(sparse.size());
    guides.reserve(sparse.size());
    for (int idx : sparse) {
        anchors.push_back(population[idx]);
        guides.push_back(population[rng.randint(static_cast<int>(population.size()))]);
    }
    return operatorDEHalfRDExMOP(prob, anchors, guides, rng);
}

struct ArchiveState {
    int size = 0;
    double nnCV = 0.0;
    double nnMean = 0.0;
    double spread = 0.0;
    double fragmentScore = 0.0;
    vector<int> fragmentIdx;
};

static vector<vector<double>> normalizeMatrix(const vector<vector<double>>& obj) {
    if (obj.empty()) return {};
    int n = static_cast<int>(obj.size());
    int m = static_cast<int>(obj[0].size());
    vector<double> mn(m, numeric_limits<double>::infinity());
    vector<double> mx(m, -numeric_limits<double>::infinity());
    for (const auto& row : obj) {
        for (int j = 0; j < m; ++j) {
            mn[j] = min(mn[j], row[j]);
            mx[j] = max(mx[j], row[j]);
        }
    }
    vector<vector<double>> out(n, vector<double>(m, 0.0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            double span = mx[j] - mn[j];
            if (span == 0.0) span = 1.0;
            out[i][j] = (obj[i][j] - mn[j]) / span;
        }
    }
    return out;
}

static vector<int> pickRepresentativesByObj(const vector<Solution>& pop, int quota) {
    vector<int> picked;
    if (quota <= 0 || pop.empty()) return picked;
    int m = static_cast<int>(pop[0].obj.size());
    auto norm = normalizeObj(pop);
    for (auto& row : norm) {
        double s = 0.0;
        for (double& v : row) {
            v = max(v, EPS);
            s += v;
        }
        for (double& v : row) v /= max(s, EPS);
    }
    auto W = uniformPoints(max(quota, m), m);
    for (auto& w : W) {
        double wn = 0.0;
        for (double& v : w) {
            v = max(v, EPS);
            wn += v * v;
        }
        wn = sqrt(max(wn, EPS));
        for (double& v : w) v /= wn;
    }
    vector<char> seen(pop.size(), 0);
    for (const auto& w : W) {
        int best = 0;
        double bestCos = -numeric_limits<double>::infinity();
        for (int i = 0; i < static_cast<int>(norm.size()); ++i) {
            double yn = sqrt(inner_product(norm[i].begin(), norm[i].end(), norm[i].begin(), 0.0));
            double c = inner_product(norm[i].begin(), norm[i].end(), w.begin(), 0.0) / max(yn, EPS);
            if (c > bestCos) {
                bestCos = c;
                best = i;
            }
        }
        if (!seen[best]) {
            seen[best] = 1;
            picked.push_back(best);
            if (static_cast<int>(picked.size()) >= quota) break;
        }
    }
    return picked;
}

static ArchiveState estimateArchiveState(const vector<Solution>& archive) {
    ArchiveState st;
    st.size = static_cast<int>(archive.size());
    if (archive.size() < 3) return st;
    auto norm = normalizeObj(archive);
    auto nn = nearestDistances(norm, 1);
    auto kth = nearestDistances(norm, min(3, static_cast<int>(norm.size()) - 1));
    double mean = accumulate(nn.begin(), nn.end(), 0.0) / nn.size();
    double var = 0.0;
    for (double v : nn) var += (v - mean) * (v - mean);
    var /= max<size_t>(1, nn.size() - 1);
    st.nnMean = mean;
    st.nnCV = sqrt(var) / max(mean, EPS);
    st.spread = accumulate(kth.begin(), kth.end(), 0.0) / max<size_t>(1, kth.size());

    auto front = firstFrontIndicesObj(norm);
    if (archive.size() >= 6 && st.nnCV >= 0.32 && front.size() >= 4) {
        vector<vector<double>> frontRaw;
        frontRaw.reserve(front.size());
        for (int idx : front) frontRaw.push_back(norm[idx]);
        auto fn = normalizeMatrix(frontRaw);
        int m = static_cast<int>(fn[0].size());
        int fN = static_cast<int>(fn.size());
        vector<char> keep(fN, 0);
        for (int k = 0; k < m; ++k) {
            int mn = 0, mx = 0;
            for (int i = 1; i < fN; ++i) {
                if (fn[i][k] < fn[mn][k]) mn = i;
                if (fn[i][k] > fn[mx][k]) mx = i;
            }
            keep[mn] = keep[mx] = 1;
        }

        vector<double> sparseScore(fN, 0.0), boundaryScore(fN, 0.0), convScore(fN, 0.0);
        for (int i = 0; i < fN; ++i) {
            double sparse = numeric_limits<double>::infinity();
            for (int j = 0; j < fN; ++j) if (i != j) sparse = min(sparse, sqrt(sqrDist(fn[i], fn[j])));
            sparseScore[i] = sparse;
            auto [mnIt, mxIt] = minmax_element(fn[i].begin(), fn[i].end());
            boundaryScore[i] = *mxIt - *mnIt;
            convScore[i] = -accumulate(fn[i].begin(), fn[i].end(), 0.0);
        }
        auto sparseNorm = normalize01(sparseScore);
        auto boundaryNorm = normalize01(boundaryScore);
        auto convNorm = normalize01(convScore);
        vector<double> rankScore(fN, 0.0);
        for (int i = 0; i < fN; ++i) {
            rankScore[i] = sparseNorm[i] + 0.18 * boundaryNorm[i] + 0.10 * convNorm[i];
        }

        int quota = min(static_cast<int>(front.size()),
                        max(4, static_cast<int>(round((0.08 + 0.10 * min(st.nnCV, 1.4)) * archive.size()))));
        int kept = accumulate(keep.begin(), keep.end(), 0);
        int remaining = quota - kept;
        if (remaining > 0) {
            vector<int> candidates;
            for (int i = 0; i < fN; ++i) if (!keep[i]) candidates.push_back(i);
            sort(candidates.begin(), candidates.end(), [&](int a, int b) { return rankScore[a] > rankScore[b]; });
            for (int i = 0; i < min(remaining, static_cast<int>(candidates.size())); ++i) keep[candidates[i]] = 1;
        }
        vector<int> selected;
        for (int i = 0; i < fN; ++i) if (keep[i]) selected.push_back(i);
        if (static_cast<int>(selected.size()) > quota) {
            sort(selected.begin(), selected.end(), [&](int a, int b) { return rankScore[a] > rankScore[b]; });
            selected.resize(quota);
        }
        double scoreSum = 0.0;
        for (int local : selected) {
            st.fragmentIdx.push_back(front[local]);
            scoreSum += rankScore[local];
        }
        if (!selected.empty()) st.fragmentScore = scoreSum / selected.size() + 0.05 * kept / max(2 * m, 1);
    }
    return st;
}

static pair<vector<Solution>, ArchiveState> updateArchive(const vector<Solution>& archive,
                                                          const vector<Solution>& cand,
                                                          int capN, double eps) {
    vector<Solution> pool = archive;
    append(pool, cand);
    int keepN = min(static_cast<int>(pool.size()), max(8, capN));
    vector<Solution> out = environmentalSelection(pool, keepN, eps);
    return {out, estimateArchiveState(out)};
}

static bool sameObjective(const Solution& a, const Solution& b) {
    if (a.obj.size() != b.obj.size()) return false;
    for (size_t i = 0; i < a.obj.size(); ++i) {
        if (fabs(a.obj[i] - b.obj[i]) > 1.0e-12) return false;
    }
    return true;
}

static vector<Solution> compressArchiveReference(const vector<Solution>& front, int capN, double eps) {
    int n = static_cast<int>(front.size());
    if (n <= capN) return front;
    int m = static_cast<int>(front[0].obj.size());
    auto norm = normalizeObj(front);
    vector<char> keep(n, 0);
    for (int k = 0; k < m; ++k) {
        int mn = 0;
        for (int i = 1; i < n; ++i) if (norm[i][k] < norm[mn][k]) mn = i;
        keep[mn] = 1;
    }
    int kept = accumulate(keep.begin(), keep.end(), 0);
    int remaining = capN - kept;
    if (remaining > 0) {
        vector<int> candidates;
        candidates.reserve(n);
        for (int i = 0; i < n; ++i) if (!keep[i]) candidates.push_back(i);
        vector<vector<double>> candNorm;
        candNorm.reserve(candidates.size());
        for (int idx : candidates) candNorm.push_back(norm[idx]);
        vector<int> refs = pickReferenceRepresentativesNorm(candNorm, remaining, m);
        for (int local : refs) keep[candidates[local]] = 1;
    }
    kept = accumulate(keep.begin(), keep.end(), 0);
    if (kept < capN) {
        vector<double> sparse(n, 0.0), conv(n, 0.0);
        for (int i = 0; i < n; ++i) {
            sparse[i] = numeric_limits<double>::infinity();
            for (int j = 0; j < n; ++j) if (i != j) sparse[i] = min(sparse[i], sqrt(sqrDist(norm[i], norm[j])));
            conv[i] = -accumulate(norm[i].begin(), norm[i].end(), 0.0);
        }
        auto sparseScore = normalize01(sparse);
        auto convScore = normalize01(conv);
        vector<int> order;
        order.reserve(n);
        for (int i = 0; i < n; ++i) if (!keep[i]) order.push_back(i);
        sort(order.begin(), order.end(), [&](int a, int b) {
            double sa = 0.72 * sparseScore[a] + 0.28 * convScore[a];
            double sb = 0.72 * sparseScore[b] + 0.28 * convScore[b];
            return sa > sb;
        });
        for (int idx : order) {
            if (kept >= capN) break;
            keep[idx] = 1;
            ++kept;
        }
    }
    vector<Solution> out;
    out.reserve(capN);
    for (int i = 0; i < n; ++i) if (keep[i]) out.push_back(front[i]);
    if (static_cast<int>(out.size()) > capN) out = environmentalSelection(out, capN, eps);
    return out;
}

static vector<Solution> updateOnlineNondominatedArchive(const vector<Solution>& archive,
                                                        const vector<Solution>& cand,
                                                        int capN, double eps) {
    vector<Solution> out = archive;
    if (cand.empty()) return out;
    for (const auto& s : cand) {
        bool rejected = false;
        for (const auto& a : out) {
            if (sameObjective(a, s) || dominates(a.obj, s.obj)) {
                rejected = true;
                break;
            }
        }
        if (rejected) continue;
        vector<Solution> kept;
        kept.reserve(out.size() + 1);
        for (const auto& a : out) {
            if (!dominates(s.obj, a.obj)) kept.push_back(a);
        }
        kept.push_back(s);
        out.swap(kept);
    }
    if (static_cast<int>(out.size()) > capN) out = compressArchiveReference(out, capN, eps);
    return out;
}

static vector<Solution> selectArchiveOutput(const vector<Solution>& mainOutput,
                                            const vector<Solution>& archive,
                                            int N, double eps, const Rng& rng) {
    if (archive.empty()) return mainOutput;
    vector<Solution> pool = mainOutput;
    append(pool, archive);
    vector<Solution> front = subset(pool, firstFrontIndices(pool));
    if (static_cast<int>(front.size()) <= N) return front;
    Rng local = rng;
    return selectionExact(front, N, local, true);
}

static vector<Solution> reinjectArchive(Problem& prob, const vector<Solution>& archive, const vector<Solution>& population,
                                        double progress, const FrontState& frontState,
                                        const ArchiveState& archiveState, int stallCounter, Rng& rng,
                                        const RootOptions& opt) {
    const RootParams& p = opt.params;
    if (archive.size() < 6 || progress < p.injectStartProgress) return {};
    double injectBias = max(0.0, frontState.nnCV - archiveState.nnCV) + 0.35 * max(progress - 0.65, 0.0) + 0.08 * stallCounter;
    if (injectBias < p.injectBiasThreshold) return {};
    int q = static_cast<int>(round((p.injectQBase + p.injectQGain * min(injectBias, 0.8)) * prob.N));
    q = min(max(q, 4), max(6, static_cast<int>(round(0.16 * prob.N))));
    int repQuota = min(max(4, static_cast<int>(round(0.60 * q))), static_cast<int>(archive.size()));
    vector<int> reps = pickRepresentativesByObj(archive, repQuota);
    if (reps.empty()) return {};
    vector<double> center(prob.D, 0.0);
    for (int idx : reps) for (int j = 0; j < prob.D; ++j) center[j] += archive[idx].dec[j];
    for (double& v : center) v /= reps.size();
    vector<double> sigma(prob.D, 0.001);
    if (reps.size() > 1) {
        for (int j = 0; j < prob.D; ++j) {
            double mean = center[j];
            double var = 0.0;
            for (int idx : reps) {
                double d = archive[idx].dec[j] - mean;
                var += d * d;
            }
            sigma[j] = max(sqrt(var / (reps.size() - 1)), 0.001);
        }
    }
    vector<vector<double>> xs(q, vector<double>(prob.D));
    vector<int> guideIdx(q), baseIdx(q);
    for (int i = 0; i < q; ++i) guideIdx[i] = reps[rng.randint(static_cast<int>(reps.size()))];
    for (int i = 0; i < q; ++i) baseIdx[i] = rng.randint(static_cast<int>(population.size()));
    vector<double> mix(q);
    for (int i = 0; i < q; ++i) mix[i] = min(max(0.28 + 0.10 * rng.randn(), 0.12), 0.48);
    vector<vector<double>> noise(q, vector<double>(prob.D, 0.0));
    double local = max(0.03, 0.10 - 0.05 * progress);
    if (progress > 0.60) {
        for (int j = 0; j < prob.D; ++j) for (int i = 0; i < q; ++i) noise[i][j] = rng.randn();
    }
    double pull = min(max(0.22 + 0.18 * max(progress - 0.60, 0.0), 0.12), 0.42);
    for (int i = 0; i < q; ++i) {
        const auto& guide = archive[guideIdx[i]].dec;
        const auto& base = population[baseIdx[i]].dec;
        for (int j = 0; j < prob.D; ++j) {
            xs[i][j] = base[j] + mix[i] * (guide[j] - base[j]) + pull * (center[j] - base[j]);
            if (progress > 0.60) xs[i][j] += local * sigma[j] * noise[i][j];
        }
    }
    return prob.evaluate(xs);
}

static vector<Solution> operatorElite2026(Problem& prob, const vector<Solution>& pop, const vector<Solution>& archive,
                                          const vector<double>& fitness, double progress,
                                          const FrontState& frontState, const ArchiveState& archiveState, Rng& rng,
                                          const RootOptions& opt) {
    const RootParams& p = opt.params;
    int N = static_cast<int>(pop.size());
    int D = prob.D;
    vector<int> order(N);
    iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int a, int b) { return fitness[a] < fitness[b]; });
    vector<const Solution*> arc;
    if (!archive.empty()) {
        if (!archiveState.fragmentIdx.empty()) {
            for (int idx : archiveState.fragmentIdx) if (idx >= 0 && idx < static_cast<int>(archive.size())) arc.push_back(&archive[idx]);
        }
        if (arc.empty()) for (const auto& s : archive) arc.push_back(&s);
    } else {
        for (const auto& s : pop) arc.push_back(&s);
    }
    int P = max(2, static_cast<int>(round(max(0.06, p.elitePbestBase - p.elitePbestDecay * progress) * N)));
    int topArc = min(static_cast<int>(arc.size()), max(4, static_cast<int>(round((p.eliteArchiveBase + p.eliteArchiveGain * progress) * arc.size()))));
    vector<double> center(D, 0.0);
    for (int i = 0; i < topArc; ++i) for (int j = 0; j < D; ++j) center[j] += arc[i]->dec[j];
    for (double& v : center) v /= max(1, topArc);
    vector<int> pbest(N), arcGuide(N);
    for (int i = 0; i < N; ++i) pbest[i] = order[rng.randint(P)];
    for (int i = 0; i < N; ++i) arcGuide[i] = rng.randint(topArc);
    vector<double> Fmain(N), Fdiff(N), CR(N);
    for (int i = 0; i < N; ++i) Fmain[i] = min(max(0.38 + 0.06 * rng.randn(), 0.18), 0.58);
    for (int i = 0; i < N; ++i) Fdiff[i] = min(max(0.14 + 0.05 * rng.randn(), 0.05), 0.28);
    double Fpull = min(max(0.18 + 0.20 * max(0.45 - frontState.nnCV, 0.0), 0.10), 0.34);
    double Farch = min(max(p.eliteFarchBase + p.eliteFarchGain * progress, 0.10), 0.36);
    auto [r1, r2] = gnR1R2FastStyle(N, rng);
    for (int i = 0; i < N; ++i) {
        CR[i] = min(max(0.86 - 0.12 * max(frontState.nnCV - 0.35, 0.0) + 0.04 * rng.randn(), 0.60), 0.98);
    }
    vector<int> jrand(N);
    for (int i = 0; i < N; ++i) jrand[i] = rng.randint(D);
    vector<vector<char>> crossMask(N, vector<char>(D, 0));
    for (int j = 0; j < D; ++j) for (int i = 0; i < N; ++i) crossMask[i][j] = (rng.rand01() < CR[i]);
    for (int i = 0; i < N; ++i) crossMask[i][jrand[i]] = 1;
    vector<vector<double>> xs(N, vector<double>(D));
    for (int i = 0; i < N; ++i) {
        const auto& guide = arc[arcGuide[i]]->dec;
        xs[i] = pop[i].dec;
        for (int j = 0; j < D; ++j) {
            double v = pop[i].dec[j]
                + Fmain[i] * (pop[pbest[i]].dec[j] - pop[i].dec[j])
                + Farch * (guide[j] - pop[i].dec[j])
                + Fdiff[i] * (pop[r1[i]].dec[j] - pop[r2[i]].dec[j])
                + Fpull * (center[j] - pop[i].dec[j]);
            if (crossMask[i][j]) xs[i][j] = v;
        }
    }
    if (progress > 0.75) {
        double sigma = 0.03 + 0.03 * max(0.35 - frontState.nnCV, 0.0);
        double pLocal = 0.06 + 0.10 * max(progress - 0.75, 0.0);
        vector<pair<int, int>> local;
        for (int j = 0; j < D; ++j) for (int i = 0; i < N; ++i) if (rng.rand01() < pLocal) local.push_back({i, j});
        for (auto [i, j] : local) {
            xs[i][j] += sigma * rng.randn();
        }
    }
    return prob.evaluate(xs);
}

static vector<Solution> eliteRefine2026(Problem& prob, const vector<Solution>& population,
                                        const vector<Solution>& archive, double eps, double progress,
                                        const FrontState& frontState, const ArchiveState& archiveState, Rng& rng) {
    vector<Solution> pool = population;
    append(pool, archive);
    if (pool.size() < 6) return {};
    auto fit = calFitness(pool, eps);
    vector<int> order(pool.size());
    iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(), [&](int a, int b) { return fit[a] < fit[b]; });

    if (!archiveState.fragmentIdx.empty()) {
        vector<char> seen(pool.size(), 0);
        vector<int> mixed;
        mixed.reserve(order.size() + archiveState.fragmentIdx.size());
        int archiveOffset = static_cast<int>(population.size());
        for (int idx : archiveState.fragmentIdx) {
            int pidx = archiveOffset + idx;
            if (pidx >= 0 && pidx < static_cast<int>(pool.size()) && !seen[pidx]) {
                seen[pidx] = 1;
                mixed.push_back(pidx);
            }
        }
        for (int idx : order) {
            if (!seen[idx]) {
                seen[idx] = 1;
                mixed.push_back(idx);
            }
        }
        order = std::move(mixed);
    }

    vector<int> refIdx;
    int refQuota = max(2, static_cast<int>(round(0.18 * prob.N * min(1.0, max(frontState.nnCV - 0.15, 0.0)))));
    if (refQuota > 0 && !order.empty()) {
        vector<Solution> orderedPool;
        orderedPool.reserve(order.size());
        for (int idx : order) orderedPool.push_back(pool[idx]);
        vector<int> refLocal = pickRepresentativesByObj(orderedPool, min(refQuota, static_cast<int>(orderedPool.size())));
        for (int local : refLocal) {
            refIdx.push_back(order[local]);
        }
    }

    int eliteCap = min(static_cast<int>(pool.size()), max(12, static_cast<int>(round((0.28 + 0.16 * max(progress - 0.68, 0.0)) * prob.N))));
    eliteCap = min(eliteCap, static_cast<int>(order.size()));
    vector<char> eliteSeen(pool.size(), 0);
    vector<int> eliteIdx;
    eliteIdx.reserve(eliteCap + refIdx.size());
    for (int i = 0; i < eliteCap; ++i) {
        int idx = order[i];
        if (!eliteSeen[idx]) {
            eliteSeen[idx] = 1;
            eliteIdx.push_back(idx);
        }
    }
    for (int idx : refIdx) {
        if (idx >= 0 && idx < static_cast<int>(pool.size()) && !eliteSeen[idx]) {
            eliteSeen[idx] = 1;
            eliteIdx.push_back(idx);
        }
    }
    int eliteSize = static_cast<int>(eliteIdx.size());
    if (eliteSize <= 0) return {};

    int qTotal = static_cast<int>(round((0.12 + 0.12 * max(progress - 0.70, 0.0) +
                                         0.10 * max(0.42 - frontState.nnCV, 0.0)) * prob.N));
    qTotal = min(max(qTotal, 8), max(10, static_cast<int>(round(0.28 * prob.N))));
    qTotal = min(qTotal, max(8, eliteSize));
    int qDE = max(4, static_cast<int>(round(0.60 * qTotal)));
    int qGS = max(4, qTotal - qDE);

    int centerN = min(eliteSize, max(3, static_cast<int>(round(0.25 * eliteSize))));
    vector<double> eliteCenter(prob.D, 0.0);
    for (int i = 0; i < centerN; ++i) {
        const auto& dec = pool[eliteIdx[i]].dec;
        for (int j = 0; j < prob.D; ++j) eliteCenter[j] += dec[j];
    }
    for (double& v : eliteCenter) v /= max(1, centerN);

    vector<vector<double>> xs;
    xs.reserve(qDE + qGS);
    int topFocus = min(eliteSize, max(4, static_cast<int>(round(0.55 * qDE))));
    vector<int> targetLocal(qDE), guideLocal(qDE), r1Local(qDE), r2Local(qDE);
    int qDETop = (qDE + 1) / 2;
    for (int i = 0; i < qDETop; ++i) targetLocal[i] = rng.randint(topFocus);
    for (int i = qDETop; i < qDE; ++i) targetLocal[i] = rng.randint(eliteSize);
    for (int i = 0; i < qDE; ++i) guideLocal[i] = rng.randint(eliteSize);
    for (int i = 0; i < qDE; ++i) r1Local[i] = rng.randint(eliteSize);
    for (int i = 0; i < qDE; ++i) r2Local[i] = rng.randint(static_cast<int>(pool.size()));
    vector<double> Fmain(qDE), Fdiff(qDE), CR(qDE);
    for (int i = 0; i < qDE; ++i) Fmain[i] = min(max(0.34 + 0.08 * rng.randn(), 0.18), 0.55);
    for (int i = 0; i < qDE; ++i) Fdiff[i] = min(max(0.18 + 0.06 * frontState.nnCV + 0.05 * rng.randn(), 0.08), 0.35);
    double Fpull = min(max(0.16 + 0.20 * max(0.45 - frontState.nnCV, 0.0), 0.08), 0.28);
    for (int i = 0; i < qDE; ++i) CR[i] = min(max(0.82 + 0.05 * rng.randn(), 0.62), 0.96);
    vector<int> jrand(qDE);
    for (int i = 0; i < qDE; ++i) jrand[i] = rng.randint(prob.D);
    vector<vector<char>> crossMask(qDE, vector<char>(prob.D, 0));
    for (int j = 0; j < prob.D; ++j) for (int i = 0; i < qDE; ++i) crossMask[i][j] = (rng.rand01() < CR[i]);
    for (int i = 0; i < qDE; ++i) crossMask[i][jrand[i]] = 1;
    for (int i = 0; i < qDE; ++i) {
        const auto& target = pool[eliteIdx[targetLocal[i]]].dec;
        const auto& guide = pool[eliteIdx[guideLocal[i]]].dec;
        const auto& r1 = pool[eliteIdx[r1Local[i]]].dec;
        const auto& r2 = pool[r2Local[i]].dec;
        vector<double> child = target;
        for (int j = 0; j < prob.D; ++j) {
            double v = target[j] + Fmain[i] * (guide[j] - target[j]) + Fdiff[i] * (r1[j] - r2[j]) +
                       Fpull * (eliteCenter[j] - target[j]);
            if (crossMask[i][j]) child[j] = v;
        }
        xs.push_back(std::move(child));
    }

    int focusN = min(eliteSize, max(6, static_cast<int>(round(0.35 * eliteSize))));
    vector<double> sigma(prob.D, 0.001);
    if (focusN > 1) {
        vector<double> mean(prob.D, 0.0);
        for (int i = 0; i < focusN; ++i) {
            const auto& dec = pool[eliteIdx[i]].dec;
            for (int j = 0; j < prob.D; ++j) mean[j] += dec[j];
        }
        for (double& v : mean) v /= focusN;
        for (int j = 0; j < prob.D; ++j) {
            double var = 0.0;
            for (int i = 0; i < focusN; ++i) {
                double d = pool[eliteIdx[i]].dec[j] - mean[j];
                var += d * d;
            }
            sigma[j] = max(sqrt(var / (focusN - 1)), 0.001);
        }
    }
    double shrink = max(0.04, 0.16 - 0.10 * progress);
    for (double& v : sigma) v *= shrink;

    int gsTop = min(eliteSize, max(4, static_cast<int>(round(0.40 * qGS))));
    double centerPullScale = 0.18 * max(0.45 - frontState.nnCV, 0.0);
    vector<int> baseLocal(qGS);
    int qGSTop = (qGS + 1) / 2;
    for (int i = 0; i < qGSTop; ++i) baseLocal[i] = rng.randint(gsTop);
    for (int i = qGSTop; i < qGS; ++i) baseLocal[i] = rng.randint(eliteSize);
    vector<vector<double>> gsNoise(qGS, vector<double>(prob.D, 0.0));
    for (int j = 0; j < prob.D; ++j) for (int i = 0; i < qGS; ++i) gsNoise[i][j] = rng.randn();
    for (int i = 0; i < qGS; ++i) {
        const auto& base = pool[eliteIdx[baseLocal[i]]].dec;
        vector<double> child = base;
        for (int j = 0; j < prob.D; ++j) {
            child[j] = base[j] + sigma[j] * gsNoise[i][j] + centerPullScale * (eliteCenter[j] - base[j]);
        }
        xs.push_back(std::move(child));
    }
    return prob.evaluate(xs);
}

static vector<Solution> coverageRepair2026(Problem& prob, const vector<Solution>& population,
                                           const vector<Solution>& archive, double progress,
                                           const FrontState& frontState, int stallCounter, Rng& rng) {
    if (population.empty() || progress < 0.45 || stallCounter < 1) return {};
    if (frontState.nnCV < 0.55 && frontState.frontRatio < 1.8) return {};
    vector<Solution> pool = population;
    append(pool, archive);
    if (pool.size() < 8) return {};
    vector<int> frontIdx = firstFrontIndices(pool);
    if (frontIdx.empty()) return {};
    auto front = subset(pool, frontIdx);
    auto norm = normalizeObj(front);
    int q = min(max(4, static_cast<int>(round((0.04 + 0.04 * min(frontState.nnCV, 1.5)) * prob.N))),
                static_cast<int>(round(0.14 * prob.N)));
    q = min(q, static_cast<int>(front.size()));
    vector<int> order(front.size());
    if (front.size() == 1) {
        order[0] = 0;
    } else {
        vector<double> sparse(front.size(), 0.0), conv(front.size(), 0.0);
        for (int i = 0; i < static_cast<int>(front.size()); ++i) {
            sparse[i] = numeric_limits<double>::infinity();
            for (int j = 0; j < static_cast<int>(front.size()); ++j) {
                if (i != j) sparse[i] = min(sparse[i], sqrt(sqrDist(norm[i], norm[j])));
            }
            conv[i] = -accumulate(norm[i].begin(), norm[i].end(), 0.0);
        }
        auto sparseScore = normalize01(sparse);
        auto convScore = normalize01(conv);
        iota(order.begin(), order.end(), 0);
        sort(order.begin(), order.end(), [&](int a, int b) {
            double sa = sparseScore[a] + 0.20 * convScore[a];
            double sb = sparseScore[b] + 0.20 * convScore[b];
            return sa > sb;
        });
    }
    auto fit = calFitness(pool, 0.05);
    vector<int> eliteOrder(pool.size());
    iota(eliteOrder.begin(), eliteOrder.end(), 0);
    sort(eliteOrder.begin(), eliteOrder.end(), [&](int a, int b) { return fit[a] < fit[b]; });
    int eliteN = min(static_cast<int>(pool.size()), max(6, static_cast<int>(round(0.25 * pool.size()))));
    vector<double> eliteMean(prob.D, 0.0), sigma(prob.D, 0.001);
    for (int i = 0; i < eliteN; ++i) {
        const auto& dec = pool[eliteOrder[i]].dec;
        for (int j = 0; j < prob.D; ++j) eliteMean[j] += dec[j];
    }
    for (double& v : eliteMean) v /= max(1, eliteN);
    if (eliteN > 1) {
        for (int j = 0; j < prob.D; ++j) {
            double var = 0.0;
            for (int i = 0; i < eliteN; ++i) {
                double d = pool[eliteOrder[i]].dec[j] - eliteMean[j];
                var += d * d;
            }
            sigma[j] = max(sqrt(var / (eliteN - 1)), 0.001);
        }
    }
    vector<vector<double>> xs(q, vector<double>(prob.D));
    double shrink = max(0.025, 0.11 - 0.07 * progress);
    double maskProb = max(0.18, 0.45 - 0.25 * progress);
    vector<int> donorLocal(q);
    for (int i = 0; i < q; ++i) donorLocal[i] = rng.randint(eliteN);
    vector<double> blend(q);
    for (int i = 0; i < q; ++i) blend[i] = 0.20 + 0.20 * rng.rand01();
    vector<vector<double>> noise(q, vector<double>(prob.D, 0.0));
    for (int j = 0; j < prob.D; ++j) for (int i = 0; i < q; ++i) noise[i][j] = rng.randn();
    vector<vector<char>> mask(q, vector<char>(prob.D, 0));
    for (int j = 0; j < prob.D; ++j) for (int i = 0; i < q; ++i) mask[i][j] = (rng.rand01() < maskProb);
    for (int i = 0; i < q; ++i) {
        const auto& anchor = front[order[i]].dec;
        const auto& donor = pool[eliteOrder[donorLocal[i]]].dec;
        xs[i] = anchor;
        for (int j = 0; j < prob.D; ++j) {
            double cand = anchor[j] + blend[i] * (donor[j] - anchor[j]) + shrink * sigma[j] * noise[i][j];
            if (mask[i][j]) xs[i][j] = cand;
            double span = prob.upper[j] - prob.lower[j];
            if (span > 0.0) {
                double y = fmod(xs[i][j] - prob.lower[j], 2.0 * span);
                if (y < 0.0) y += 2.0 * span;
                if (y > span) y = 2.0 * span - y;
                xs[i][j] = prob.lower[j] + y;
            } else {
                xs[i][j] = prob.lower[j];
            }
        }
    }
    return prob.evaluate(xs);
}

static vector<vector<double>> loadPF(const string& root, int problem, int M) {
    stringstream ss;
    ss << root << "/MaOP" << problem << "_F" << M << ".txt";
    ifstream in(ss.str());
    vector<vector<double>> pf;
    string line;
    while (getline(in, line)) {
        istringstream is(line);
        vector<double> row(M);
        bool ok = true;
        for (int j = 0; j < M; ++j) ok = ok && static_cast<bool>(is >> row[j]);
        if (ok) pf.push_back(row);
    }
    return pf;
}

static vector<vector<double>> readMatrix(const string& path, int cols) {
    ifstream in(path);
    if (!in) throw runtime_error("Cannot open eval input: " + path);
    vector<vector<double>> rows;
    string line;
    while (getline(in, line)) {
        if (line.empty()) continue;
        for (char& c : line) {
            if (c == ',') c = ' ';
        }
        istringstream is(line);
        vector<double> row(cols);
        bool ok = true;
        for (int j = 0; j < cols; ++j) ok = ok && static_cast<bool>(is >> row[j]);
        if (ok) rows.push_back(std::move(row));
    }
    return rows;
}

static void writeMatrix(const string& path, const vector<vector<double>>& rows) {
    ostream* out = &cout;
    ofstream file;
    if (!path.empty()) {
        file.open(path);
        if (!file) throw runtime_error("Cannot open eval output: " + path);
        out = &file;
    }
    *out << setprecision(17);
    for (const auto& row : rows) {
        for (int j = 0; j < static_cast<int>(row.size()); ++j) {
            if (j) *out << '\t';
            *out << row[j];
        }
        *out << '\n';
    }
}

static double igd(const vector<Solution>& pop, const vector<vector<double>>& pf) {
    if (pf.empty()) return numeric_limits<double>::quiet_NaN();
    auto front = subset(pop, firstFrontIndices(pop));
    if (front.empty()) return numeric_limits<double>::quiet_NaN();
    double sum = 0.0;
    for (const auto& p : pf) {
        double best = numeric_limits<double>::infinity();
        for (const auto& s : front) best = min(best, sqrt(sqrDist(p, s.obj)));
        sum += best;
    }
    return sum / pf.size();
}

struct Snapshot {
    int index = 0;
    int fe = 0;
    vector<Solution> pop;
};

struct RunResult {
    vector<double> igdHist;
    vector<Solution> output;
    vector<Snapshot> snapshots;
};

static RunResult runRDExMOP(const Config& cfg, const vector<vector<double>>& pf) {
    Rng rng(cfg.seed);
    Problem prob(cfg.problem, cfg.N, cfg.M, cfg.D, cfg.maxFE);
    const double eps = 0.05;
    vector<vector<double>> init(cfg.N, vector<double>(cfg.D));
    for (int j = 0; j < cfg.D; ++j) for (int i = 0; i < cfg.N; ++i) init[i][j] = rng.rand01();
    vector<Solution> Popul = prob.evaluate(init);
    auto sel0 = selectionRDExMOP(Popul, cfg.N, rng);
    vector<Solution> Population = sel0.first;
    int nND = sel0.second;
    vector<Solution> Population2 = Population;
    vector<double> hist;
    hist.reserve(cfg.save);
    vector<Snapshot> snapshots;
    snapshots.reserve(cfg.save);

    auto record = [&]() {
        int target = static_cast<int>(ceil(static_cast<double>(cfg.save) * prob.FE / cfg.maxFE));
        target = min(max(target, 1), cfg.save);
        int next = static_cast<int>(hist.size()) + 1;
        int idx = min(min(cfg.save, next), target) - 1;
        double val = igd(Population, pf);
        Snapshot shot{idx + 1, prob.FE, Population};
        if (idx == static_cast<int>(hist.size())) {
            hist.push_back(val);
            snapshots.push_back(std::move(shot));
        } else {
            hist[idx] = val;
            snapshots[idx] = std::move(shot);
        }
    };

    while (true) {
        record();
        if (prob.FE >= cfg.maxFE) break;
        vector<Solution> NewPopulation = explorationRDExMOP(prob, Population, Popul, nND, rng);
        vector<Solution> NewPopul;

        if (prob.FE >= cfg.maxFE * 0.5 && rng.rand01() < 0.5) {
            vector<Solution> tmp = Population2;
            append(tmp, NewPopulation);
            Population2 = environmentalSelection(tmp, cfg.N, eps);
            vector<double> fitP2 = calFitness(Population2, eps);
            for (double& v : fitP2) v = -v;
            vector<int> mating = tournamentSelection(2, cfg.N, fitP2, rng);
            vector<Solution> parents = subset(Population2, mating);
            vector<Solution> fitSource = subset(Popul, mating);
            vector<double> pfit = calFitness(fitSource, eps);
            for (double& v : pfit) v = -v;
            NewPopul = operatorDERDExMOP(prob, parents, pfit, rng);
            tmp = Population2;
            append(tmp, NewPopul);
            Population2 = environmentalSelection(tmp, cfg.N, eps);
        } else {
            vector<Solution> tmp = Popul;
            append(tmp, NewPopulation);
            Popul = environmentalSelection(tmp, cfg.N, eps);
            vector<double> fit = calFitness(Popul, eps);
            for (double& v : fit) v = -v;
            vector<int> mating = tournamentSelection(2, cfg.N, fit, rng);
            vector<Solution> parents = subset(Popul, mating);
            vector<double> pfit = calFitness(parents, eps);
            for (double& v : pfit) v = -v;
            NewPopul = operatorDERDExMOP(prob, parents, pfit, rng);
            tmp = Popul;
            append(tmp, NewPopul);
            Popul = environmentalSelection(tmp, cfg.N, eps);
        }

        vector<Solution> merged = Population;
        append(merged, NewPopul);
        append(merged, NewPopulation);
        append(merged, Population2);
        auto sel = selectionRDExMOP(merged, cfg.N, rng);
        Population = sel.first;
        nND = sel.second;
        Population2 = Population;
    }

    return {hist, Population, snapshots};
}

static RunResult runRDExMOPESD(const Config& cfg, const vector<vector<double>>& pf) {
    Rng rng(cfg.seed);
    Problem prob(cfg.problem, cfg.N, cfg.M, cfg.D, cfg.maxFE);
    const double eps = 0.05;
    vector<vector<double>> init(cfg.N, vector<double>(cfg.D));
    for (int j = 0; j < cfg.D; ++j) for (int i = 0; i < cfg.N; ++i) init[i][j] = rng.rand01();
    vector<Solution> Popul = prob.evaluate(init);
    ESDState esd;
    auto sel0 = selectionRDExMOPESD(Popul, cfg.N, rng, esd);
    vector<Solution> Population = sel0.first;
    int nND = sel0.second;
    vector<Solution> Population2 = Population;
    vector<double> hist;
    hist.reserve(cfg.save);
    vector<Snapshot> snapshots;
    snapshots.reserve(cfg.save);

    auto record = [&]() {
        int target = static_cast<int>(ceil(static_cast<double>(cfg.save) * prob.FE / cfg.maxFE));
        target = min(max(target, 1), cfg.save);
        int next = static_cast<int>(hist.size()) + 1;
        int idx = min(min(cfg.save, next), target) - 1;
        double val = igd(Population, pf);
        Snapshot shot{idx + 1, prob.FE, Population};
        if (idx == static_cast<int>(hist.size())) {
            hist.push_back(val);
            snapshots.push_back(std::move(shot));
        } else {
            hist[idx] = val;
            snapshots[idx] = std::move(shot);
        }
    };

    while (true) {
        record();
        if (prob.FE >= cfg.maxFE) break;
        vector<Solution> NewPopulation = explorationRDExMOP(prob, Population, Popul, nND, rng);
        vector<Solution> NewPopul;

        if (prob.FE >= cfg.maxFE * 0.5 && rng.rand01() < 0.5) {
            vector<Solution> tmp = Population2;
            append(tmp, NewPopulation);
            Population2 = environmentalSelection(tmp, cfg.N, eps);
            vector<double> fitP2 = calFitness(Population2, eps);
            for (double& v : fitP2) v = -v;
            vector<int> mating = tournamentSelection(2, cfg.N, fitP2, rng);
            vector<Solution> parents = subset(Population2, mating);
            vector<Solution> fitSource = subset(Popul, mating);
            vector<double> pfit = calFitness(fitSource, eps);
            for (double& v : pfit) v = -v;
            NewPopul = operatorDERDExMOP(prob, parents, pfit, rng);
            tmp = Population2;
            append(tmp, NewPopul);
            Population2 = environmentalSelection(tmp, cfg.N, eps);
        } else {
            vector<Solution> tmp = Popul;
            append(tmp, NewPopulation);
            Popul = environmentalSelection(tmp, cfg.N, eps);
            vector<double> fit = calFitness(Popul, eps);
            for (double& v : fit) v = -v;
            vector<int> mating = tournamentSelection(2, cfg.N, fit, rng);
            vector<Solution> parents = subset(Popul, mating);
            vector<double> pfit = calFitness(parents, eps);
            for (double& v : pfit) v = -v;
            NewPopul = operatorDERDExMOP(prob, parents, pfit, rng);
            tmp = Popul;
            append(tmp, NewPopul);
            Popul = environmentalSelection(tmp, cfg.N, eps);
        }

        vector<Solution> merged = Population;
        append(merged, NewPopul);
        append(merged, NewPopulation);
        append(merged, Population2);
        auto sel = selectionRDExMOPESD(merged, cfg.N, rng, esd);
        Population = sel.first;
        nND = sel.second;
        Population2 = Population;
    }

    return {hist, Population, snapshots};
}

static RunResult runRDExMOPOnline(const Config& cfg, const vector<vector<double>>& pf) {
    Rng rng(cfg.seed);
    Problem prob(cfg.problem, cfg.N, cfg.M, cfg.D, cfg.maxFE);
    const double eps = 0.05;
    vector<vector<double>> init(cfg.N, vector<double>(cfg.D));
    for (int j = 0; j < cfg.D; ++j) for (int i = 0; i < cfg.N; ++i) init[i][j] = rng.rand01();
    vector<Solution> Popul = prob.evaluate(init);
    auto sel0 = selectionRDExMOP(Popul, cfg.N, rng);
    vector<Solution> Population = sel0.first;
    int nND = sel0.second;
    vector<Solution> Population2 = Population;
    vector<double> hist;
    hist.reserve(cfg.save);
    vector<Snapshot> snapshots;
    snapshots.reserve(cfg.save);

    auto record = [&]() {
        int target = static_cast<int>(ceil(static_cast<double>(cfg.save) * prob.FE / cfg.maxFE));
        target = min(max(target, 1), cfg.save);
        int next = static_cast<int>(hist.size()) + 1;
        int idx = min(min(cfg.save, next), target) - 1;
        double val = igd(Population, pf);
        Snapshot shot{idx + 1, prob.FE, Population};
        if (idx == static_cast<int>(hist.size())) {
            hist.push_back(val);
            snapshots.push_back(std::move(shot));
        } else {
            hist[idx] = val;
            snapshots[idx] = std::move(shot);
        }
    };

    while (true) {
        record();
        if (prob.FE >= cfg.maxFE) break;
        vector<Solution> NewPopulation = explorationRDExMOP(prob, Population, Popul, nND, rng);
        vector<Solution> NewPopul;
        if (prob.FE >= cfg.maxFE * 0.5 && rng.rand01() < 0.5) {
            vector<Solution> tmp = Population2;
            append(tmp, NewPopulation);
            Population2 = environmentalSelection(tmp, cfg.N, eps);
            NewPopul = operatorDERDExMOPOnline(prob, Population2, cfg.N, eps, rng);
        } else {
            vector<Solution> tmp = Popul;
            append(tmp, NewPopulation);
            Popul = environmentalSelection(tmp, cfg.N, eps);
            NewPopul = operatorDERDExMOPOnline(prob, Popul, cfg.N, eps, rng);
        }

        vector<Solution> merged = Population;
        append(merged, NewPopul);
        append(merged, NewPopulation);
        append(merged, Population2);
        auto sel = selectionRDExMOP(merged, cfg.N, rng);
        Population = sel.first;
        nND = sel.second;
        Population2 = Population;
    }

    return {hist, Population, snapshots};
}

static RunResult runRDE2026RootCore(const Config& cfg, const vector<vector<double>>& pf, const RootOptions& opt) {
    Rng rng(cfg.seed);
    Problem prob(cfg.problem, cfg.N, cfg.M, cfg.D, cfg.maxFE);
    const double eps = 0.05;
    vector<vector<double>> init(cfg.N, vector<double>(cfg.D));
    for (int j = 0; j < cfg.D; ++j) for (int i = 0; i < cfg.N; ++i) init[i][j] = rng.rand01();
    vector<Solution> Popul = prob.evaluate(init);
    auto sel0 = selectionRootConfigured(Popul, cfg.N, 0.0, eps, rng, opt);
    vector<Solution> Population = sel0.first;
    FrontState frontState = sel0.second;
    int nND = static_cast<int>(Population.size());
    vector<Solution> Population2 = Population;
    auto arc0 = updateArchive({}, Population, 3 * cfg.N, eps);
    vector<Solution> Archive = arc0.first;
    ArchiveState archiveState = arc0.second;
    vector<Solution> TraceArchive;
    int nextOnlineArchiveGuideFE = 0;
    auto refreshOnlineArchive = [&](const vector<Solution>& cand) {
        if (opt.archiveMode <= 0 || cand.empty()) return;
        int traceCap = max(cfg.N, static_cast<int>(round(opt.params.traceArchiveCapFactor * cfg.N)));
        int guideCap = max(cfg.N, static_cast<int>(round(opt.params.guideArchiveCapFactor * cfg.N)));
        if (opt.archiveMode != 1) {
            traceCap = max(cfg.N, static_cast<int>(round(max(opt.params.traceArchiveCapFactor, 4.0) * cfg.N)));
            guideCap = max(cfg.N, static_cast<int>(round(opt.params.guideArchiveCapFactor * cfg.N)));
        }
        TraceArchive = updateOnlineNondominatedArchive(TraceArchive, cand, traceCap, eps);
        bool refreshGuide = Archive.empty() || prob.FE >= nextOnlineArchiveGuideFE;
        if (refreshGuide) {
            if (opt.archiveMode == 1) {
                Archive = TraceArchive;
            } else {
                Archive = compressArchiveReference(TraceArchive, min(static_cast<int>(TraceArchive.size()), guideCap), eps);
            }
            archiveState = estimateArchiveState(Archive);
            nextOnlineArchiveGuideFE = prob.FE + max(1500, static_cast<int>(round(opt.params.archiveRefreshFEFrac * cfg.maxFE)));
        }
    };
    refreshOnlineArchive(Population);
    vector<Solution> StablePop;
    FrontState stableState = frontState;
    StableOutputMemory stableMemory;
    vector<Solution> OutputPopulation = Population;
    bool polished = false;
    int stallCounter = 0;
    FrontState prevState = frontState;
    ESDState rootESD;
    ESDState rootESP;
    double nextArchiveFE = 0.58 * cfg.maxFE;
    double nextInjectFE = 0.54 * cfg.maxFE;
    double nextRefineFE = 0.72 * cfg.maxFE;
    double nextRepairFE = 0.58 * cfg.maxFE;
    vector<double> hist;
    hist.reserve(cfg.save);
    vector<Snapshot> snapshots;
    snapshots.reserve(cfg.save);
    auto record = [&]() {
        int target = static_cast<int>(ceil(static_cast<double>(cfg.save) * prob.FE / cfg.maxFE));
        target = min(max(target, 1), cfg.save);
        int next = static_cast<int>(hist.size()) + 1;
        int idx = min(min(cfg.save, next), target) - 1;
        double val = igd(OutputPopulation, pf);
        Snapshot shot{idx + 1, prob.FE, OutputPopulation};
        if (idx == static_cast<int>(hist.size())) {
            hist.push_back(val);
            snapshots.push_back(std::move(shot));
        } else {
            hist[idx] = val;
            snapshots[idx] = std::move(shot);
        }
    };

    while (true) {
        record();
        if (prob.FE >= cfg.maxFE) break;
        double progress = static_cast<double>(prob.FE) / cfg.maxFE;
        double eliteProb = 0.45 + 0.20 * max(progress - 0.60, 0.0) + 0.10 * min(max(stallCounter - 1, 0), 2);
        eliteProb = min(0.85, max(0.45, eliteProb));
        bool useEliteBranch = opt.useEliteBranch && prob.FE >= cfg.maxFE * 0.5 && rng.rand01() < eliteProb;
        double unevenBoost = min(0.18, 0.22 * max(frontState.nnCV - 0.35, 0.0));
        bool doExplore = opt.useExplore && rng.rand01() < min(0.98, max(0.70, 0.92 - 0.18 * progress + unevenBoost));
        vector<Solution> NewPopulation;
        vector<Solution> NewPopul;
        if (useEliteBranch) {
            if (doExplore) NewPopulation = exploration2026(prob, Population, Popul, nND, rng);
            refreshOnlineArchive(NewPopulation);
            vector<Solution> tmp = Population2;
            append(tmp, NewPopulation);
            Population2 = environmentalSelection(tmp, cfg.N, eps);
            auto fit = calFitness(Population2, eps);
            vector<int> mating(cfg.N);
            for (int i = 0; i < cfg.N; ++i) {
                int a = rng.randint(static_cast<int>(Population2.size()));
                int b = rng.randint(static_cast<int>(Population2.size()));
                mating[i] = (fit[a] > fit[b]) ? a : b;
            }
            vector<Solution> parents = subset(Population2, mating);
            vector<double> pfit;
            pfit.reserve(parents.size());
            for (int idx : mating) pfit.push_back(-fit[idx]);
            bool highSpan = frontState.spanRatio >= 2.50;
            bool useDiversityDE = frontState.nnCV > 0.52 && (progress < 0.82 || (progress < 0.90 && highSpan));
            if (useDiversityDE || !opt.useEliteOperator) NewPopul = operatorDE2026(prob, parents, pfit, rng);
            else NewPopul = operatorElite2026(prob, parents, Archive, pfit, progress, frontState, archiveState, rng, opt);
            refreshOnlineArchive(NewPopul);
            tmp = Population2;
            append(tmp, NewPopul);
            Population2 = environmentalSelection(tmp, cfg.N, eps);
        } else {
            if (doExplore) NewPopulation = exploration2026(prob, Population, Popul, nND, rng);
            refreshOnlineArchive(NewPopulation);
            vector<Solution> tmp = Popul;
            append(tmp, NewPopulation);
            Popul = environmentalSelection(tmp, cfg.N, eps);
            auto fit = calFitness(Popul, eps);
            vector<int> mating(cfg.N);
            for (int i = 0; i < cfg.N; ++i) {
                int a = rng.randint(static_cast<int>(Popul.size()));
                int b = rng.randint(static_cast<int>(Popul.size()));
                mating[i] = (fit[a] > fit[b]) ? a : b;
            }
            vector<Solution> parents = subset(Popul, mating);
            vector<double> pfit;
            for (int idx : mating) pfit.push_back(-fit[idx]);
            NewPopul = operatorDE2026(prob, parents, pfit, rng);
            refreshOnlineArchive(NewPopul);
            tmp = Popul;
            append(tmp, NewPopul);
            Popul = environmentalSelection(tmp, cfg.N, eps);
        }

        vector<Solution> RefinedPopulation, MemoryPopulation, RepairPopulation;
        if (opt.useRefine && prob.FE >= nextRefineFE && progress >= 0.72 &&
                (frontState.nnCV < 0.45 || stallCounter >= 1 || frontState.needsPolish)) {
            vector<Solution> tmp = Population;
            append(tmp, Population2);
            RefinedPopulation = eliteRefine2026(prob, tmp, Archive, eps, progress, frontState, archiveState, rng);
            refreshOnlineArchive(RefinedPopulation);
            nextRefineFE += 0.10 * cfg.maxFE;
        }
        if (opt.useMemory && prob.FE >= nextInjectFE) {
            MemoryPopulation = reinjectArchive(prob, Archive, Population, progress, frontState, archiveState, stallCounter, rng, opt);
            refreshOnlineArchive(MemoryPopulation);
            nextInjectFE += 0.12 * cfg.maxFE;
        }
        if (opt.useRepair && prob.FE >= nextRepairFE && (frontState.nnCV > 0.50 || stallCounter >= 1)) {
            RepairPopulation = coverageRepair2026(prob, Population, Archive, progress, frontState, stallCounter, rng);
            refreshOnlineArchive(RepairPopulation);
            nextRepairFE += 0.10 * cfg.maxFE;
        }

        vector<Solution> merged = Population;
        append(merged, NewPopul);
        append(merged, NewPopulation);
        append(merged, Population2);
        append(merged, RefinedPopulation);
        append(merged, MemoryPopulation);
        append(merged, RepairPopulation);
        pair<vector<Solution>, FrontState> sel;
        if (opt.useESP) {
            sel = selectionRootAdaptiveESP(merged, cfg.N, progress, eps, rng, rootESP, opt);
        } else if (opt.useESD) {
            sel = selectionRootAdaptiveESD(merged, cfg.N, progress, eps, rng, rootESD);
        } else {
            sel = selectionRootConfigured(merged, cfg.N, progress, eps, rng, opt);
        }
        Population = sel.first;
        FrontState newState = sel.second;
        nND = static_cast<int>(Population.size());
        double stateDelta = fabs(newState.nnCV - prevState.nnCV) + fabs(newState.frontRatio - prevState.frontRatio);
        stallCounter = (stateDelta < 0.03) ? stallCounter + 1 : 0;
        if (opt.archiveMode == 0 && opt.useArchiveUpdate && prob.FE >= nextArchiveFE &&
                (newState.nnCV > 0.42 || stallCounter >= 2 || progress >= 0.88)) {
            vector<Solution> cand = Population;
            append(cand, Population2);
            append(cand, RefinedPopulation);
            append(cand, MemoryPopulation);
            append(cand, RepairPopulation);
            auto arc = updateArchive(Archive, cand, 3 * cfg.N, eps);
            Archive = arc.first;
            archiveState = arc.second;
            nextArchiveFE += 0.06 * cfg.maxFE;
        }
        if (opt.usePolish && !polished && newState.needsPolish) {
            vector<Solution> tmp = merged;
            append(tmp, Archive);
            append(tmp, Popul);
            Population = selectionExact(tmp, cfg.N, rng, false);
            nND = static_cast<int>(Population.size());
            polished = true;
        }
        if (opt.useStableOutput) OutputPopulation = stableOutput2026(StablePop, stableState, stableMemory, Population, newState, progress);
        else OutputPopulation = Population;
        if (opt.useSidecar) OutputPopulation = outputSidecar2026(OutputPopulation, merged, cfg.N, progress, eps, rng);
        if (opt.useArchiveOutput && progress >= opt.params.archiveOutputProgress) {
            OutputPopulation = selectArchiveOutput(OutputPopulation, TraceArchive, cfg.N, eps, rng);
        }
        frontState = newState;
        prevState = newState;
        Population2 = Population;
    }
    return {hist, OutputPopulation, snapshots};
}

static RunResult runRDE2026Root(const Config& cfg, const vector<vector<double>>& pf) {
    RootOptions opt;
    return runRDE2026RootCore(cfg, pf, opt);
}

static RunResult runRDE2026RootESD(const Config& cfg, const vector<vector<double>>& pf) {
    RootOptions opt;
    opt.useESD = true;
    return runRDE2026RootCore(cfg, pf, opt);
}

static RunResult runRDE2026RootESP(const Config& cfg, const vector<vector<double>>& pf) {
    RootOptions opt;
    opt.useESP = true;
    return runRDE2026RootCore(cfg, pf, opt);
}

static RunResult runRDE2026RootROOTESP(const Config& cfg, const vector<vector<double>>& pf) {
    RootOptions opt;
    opt.useESP = true;
    opt.useStableOutput = false;
    opt.useRepair = false;
    opt.usePolish = false;
    opt.useRefine = false;
    return runRDE2026RootCore(cfg, pf, opt);
}

static RootOptions rootESPBaseOptions() {
    RootOptions opt;
    opt.useESP = true;
    opt.useStableOutput = false;
    opt.useRepair = false;
    opt.usePolish = false;
    opt.useRefine = false;
    return opt;
}

static RunResult runRDE2026RootPaperArchive(const Config& cfg, const vector<vector<double>>& pf) {
    RootOptions opt = rootESPBaseOptions();
    opt.archiveMode = 1;
    opt.useArchiveOutput = true;
    return runRDE2026RootCore(cfg, pf, opt);
}

static RunResult runRDE2026MOPTPE_R05_T0070(const Config& cfg, const vector<vector<double>>& pf) {
    RootOptions opt = rootESPBaseOptions();
    opt.archiveMode = 1;
    opt.useArchiveOutput = true;
    opt.params = fixedR05T0070RootParams();
    return runRDE2026RootCore(cfg, pf, opt);
}

static RunResult runRDE2026RootTraceGuideArchive(const Config& cfg, const vector<vector<double>>& pf) {
    RootOptions opt = rootESPBaseOptions();
    opt.archiveMode = 2;
    opt.useArchiveOutput = true;
    return runRDE2026RootCore(cfg, pf, opt);
}

static void applyRootAblationToken(RootOptions& opt, const string& token) {
    if (token == "NoDirectional") opt.useDirectional = false;
    else if (token == "NoExplore") opt.useExplore = false;
    else if (token == "NoEliteBranch") opt.useEliteBranch = false;
    else if (token == "NoEliteOperator") opt.useEliteOperator = false;
    else if (token == "NoRefine") opt.useRefine = false;
    else if (token == "NoMemory") opt.useMemory = false;
    else if (token == "NoRepair") opt.useRepair = false;
    else if (token == "NoArchiveUpdate") opt.useArchiveUpdate = false;
    else if (token == "NoPolish") opt.usePolish = false;
    else if (token == "NoStableOutput") opt.useStableOutput = false;
    else if (token == "NoSidecar") opt.useSidecar = false;
    else if (token == "ESP") opt.useESP = true;
    else if (token == "ESD") opt.useESD = true;
    else if (token == "PaperArchive") {
        opt.archiveMode = 1;
        opt.useArchiveOutput = true;
    }
    else if (token == "TraceGuideArchive") {
        opt.archiveMode = 2;
        opt.useArchiveOutput = true;
    }
    else if (token == "NoArchiveOutput") opt.useArchiveOutput = false;
    else throw runtime_error("Unknown RDE2026Root ablation token: " + token);
}

static RunResult runRDE2026RootAblation(const Config& cfg, const vector<vector<double>>& pf, const string& ablation) {
    RootOptions opt;
    size_t start = 0;
    while (start <= ablation.size()) {
        size_t pos = ablation.find('_', start);
        string token = ablation.substr(start, pos == string::npos ? string::npos : pos - start);
        if (token.empty()) throw runtime_error("Empty RDE2026Root ablation token in: " + ablation);
        applyRootAblationToken(opt, token);
        if (pos == string::npos) break;
        start = pos + 1;
    }
    return runRDE2026RootCore(cfg, pf, opt);
}

static void writeOutputs(const Config& cfg, const RunResult& result, double runtimeSeconds) {
    filesystem::path algDir = filesystem::path(cfg.outDir) / "Data" / cfg.algName;
    filesystem::create_directories(algDir);
    string base = cfg.algName + "_MaOP" + to_string(cfg.problem) + "_M" + to_string(cfg.M) + "_D" + to_string(cfg.D) + "_" + to_string(cfg.run);
    {
        ofstream out(algDir / (base + "_igd.tsv"));
        out << setprecision(17);
        out << "index\tIGD\n";
        for (int i = 0; i < static_cast<int>(result.igdHist.size()); ++i) out << (i + 1) << '\t' << result.igdHist[i] << '\n';
    }
    {
        ofstream out(algDir / (base + "_objs.tsv"));
        out << setprecision(17);
        for (int j = 0; j < cfg.M; ++j) {
            if (j) out << '\t';
            out << "f" << (j + 1);
        }
        out << '\n';
        auto front = subset(result.output, firstFrontIndices(result.output));
        for (const auto& s : front) {
            for (int j = 0; j < cfg.M; ++j) {
                if (j) out << '\t';
                out << s.obj[j];
            }
            out << '\n';
        }
    }
    {
        ofstream out(algDir / (base + "_snapshots.tsv"));
        out << setprecision(17);
        out << "index\tFE\tmember";
        for (int j = 0; j < cfg.D; ++j) out << "\tx" << (j + 1);
        for (int j = 0; j < cfg.M; ++j) out << "\tf" << (j + 1);
        out << '\n';
        for (const auto& shot : result.snapshots) {
            for (int i = 0; i < static_cast<int>(shot.pop.size()); ++i) {
                out << shot.index << '\t' << shot.fe << '\t' << (i + 1);
                for (double v : shot.pop[i].dec) out << '\t' << v;
                for (double v : shot.pop[i].obj) out << '\t' << v;
                out << '\n';
            }
        }
    }
    {
        ofstream out(algDir / (base + "_meta.tsv"));
        out << "key\tvalue\n";
        out << "algorithm\t" << cfg.algName << '\n';
        out << "problem\t" << cfg.problem << '\n';
        out << "run\t" << cfg.run << '\n';
        out << "N\t" << cfg.N << '\n';
        out << "M\t" << cfg.M << '\n';
        out << "D\t" << cfg.D << '\n';
        out << "maxFE\t" << cfg.maxFE << '\n';
        out << "save\t" << cfg.save << '\n';
        out << "seed\t" << cfg.seed << '\n';
        out << "runtime_seconds\t" << setprecision(17) << runtimeSeconds << '\n';
    }
}

int main(int argc, char** argv) {
    Config cfg = parseArgs(argc, argv);
    if (cfg.M != 3 || (cfg.D != 7 && cfg.D != 10)) {
        cerr << "This C++ port is validated for MOP M=3 with D=7 or D=10 only." << endl;
        return 2;
    }
    if (!cfg.evalFile.empty()) {
        try {
            Problem prob(cfg.problem, cfg.N, cfg.M, cfg.D, cfg.maxFE);
            vector<vector<double>> x = readMatrix(cfg.evalFile, cfg.D);
            vector<vector<double>> y;
            y.reserve(x.size());
            for (const auto& row : x) y.push_back(prob.calObj(row));
            writeMatrix(cfg.evalOut, y);
            return 0;
        } catch (const exception& ex) {
            cerr << ex.what() << endl;
            return 4;
        }
    }
    auto pf = loadPF(cfg.pfRoot, cfg.problem, cfg.M);
    if (pf.empty()) {
        cerr << "Cannot load Pareto front for MaOP" << cfg.problem << " from " << cfg.pfRoot << endl;
        return 3;
    }
    auto start = chrono::steady_clock::now();
    RunResult result;
    result = runRDE2026MOPTPE_R05_T0070(cfg, pf);
    auto end = chrono::steady_clock::now();
    double sec = chrono::duration<double>(end - start).count();
    writeOutputs(cfg, result, sec);
    cout << "DONE"
         << "\talgorithm=" << cfg.algName
         << "\tproblem=MaOP" << cfg.problem
         << "\trun=" << cfg.run
         << "\tseed=" << cfg.seed
         << "\tseconds=" << fixed << setprecision(3) << sec
         << "\tfinal_igd=" << setprecision(12) << result.igdHist.back()
         << endl;
    return 0;
}
