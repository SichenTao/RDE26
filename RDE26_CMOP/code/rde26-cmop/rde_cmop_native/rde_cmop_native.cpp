/*
  Native C++ RDE runner for the CEC2026 CMOP SDC benchmark.

  This file intentionally avoids PlatEMO, MATLAB objects, and MEX.  It keeps
  the algorithm state in contiguous vectors and calls the C++ CMOP evaluator
  directly.
*/

#include "../../../benchmark/cpp/cmop_cec2026/cmop_cec2026_test.h"

#include <algorithm>
#include <chrono>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <limits>
#include <numeric>
#include <random>
#include <string>
#include <vector>

struct Options
{
    int problem = 1;
    int run = 1;
    int seed = 20260608;
    int popsize = 100;
    int dim = 30;
    int maxFE = 200000;
    int save = 1000;
    bool toprace = true;
    bool igd = true;
    bool stateTrace = false;
    int stateTraceLimit = 0;
    std::string pfRoot;
};

struct MatlabRng
{
    static const int N = 624;
    static const int M = 397;
    unsigned long mt[N];
    int mti = N + 1;

    explicit MatlabRng(unsigned long seed = 5489UL) { init(seed); }

    void init(unsigned long seed)
    {
        mt[0] = seed & 0xffffffffUL;
        for (mti = 1; mti < N; ++mti)
        {
            mt[mti] = (1812433253UL * (mt[mti - 1] ^ (mt[mti - 1] >> 30)) + (unsigned long)mti);
            mt[mti] &= 0xffffffffUL;
        }
    }

    unsigned long uint32()
    {
        static const unsigned long mag01[2] = {0x0UL, 0x9908b0dfUL};
        unsigned long y;
        if (mti >= N)
        {
            int kk = 0;
            for (; kk < N - M; ++kk)
            {
                y = (mt[kk] & 0x80000000UL) | (mt[kk + 1] & 0x7fffffffUL);
                mt[kk] = mt[kk + M] ^ (y >> 1) ^ mag01[y & 0x1UL];
            }
            for (; kk < N - 1; ++kk)
            {
                y = (mt[kk] & 0x80000000UL) | (mt[kk + 1] & 0x7fffffffUL);
                mt[kk] = mt[kk + (M - N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
            }
            y = (mt[N - 1] & 0x80000000UL) | (mt[0] & 0x7fffffffUL);
            mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ mag01[y & 0x1UL];
            mti = 0;
        }
        y = mt[mti++];
        y ^= (y >> 11);
        y ^= (y << 7) & 0x9d2c5680UL;
        y ^= (y << 15) & 0xefc60000UL;
        y ^= (y >> 18);
        return y;
    }

    double rand()
    {
        const unsigned long a = uint32() >> 5;
        const unsigned long b = uint32() >> 6;
        return ((double)a * 67108864.0 + (double)b) * (1.0 / 9007199254740992.0);
    }

    int randi(int highInclusive)
    {
        int value = (int)std::floor(rand() * (double)highInclusive) + 1;
        return value > highInclusive ? highInclusive : value;
    }
};

struct Population
{
    int n = 0;
    int dim = 30;
    int nobj = 2;
    int ncon = 3;
    std::vector<double> dec;
    std::vector<double> obj;
    std::vector<double> con;
    std::vector<double> fitness;

    Population() {}
    Population(int n_, int dim_, int nobj_, int ncon_)
        : n(n_), dim(dim_), nobj(nobj_), ncon(ncon_),
          dec((size_t)n_ * dim_), obj((size_t)n_ * nobj_),
          con((size_t)n_ * ncon_), fitness(n_, 0.0)
    {
    }
};

struct ReferenceFront
{
    int n = 0;
    int nobj = 0;
    bool available = false;
    std::vector<double> obj;
};

static double now_seconds()
{
    using clock = std::chrono::steady_clock;
    static const auto t0 = clock::now();
    auto dt = clock::now() - t0;
    return std::chrono::duration<double>(dt).count();
}

static Options parse_options(int argc, char **argv)
{
    Options opt;
    for (int i = 1; i < argc; ++i)
    {
        auto next_int = [&](int &value) {
            if (i + 1 >= argc)
            {
                std::fprintf(stderr, "Missing value for %s\n", argv[i]);
                std::exit(2);
            }
            value = std::atoi(argv[++i]);
        };
        auto next_string = [&](std::string &value) {
            if (i + 1 >= argc)
            {
                std::fprintf(stderr, "Missing value for %s\n", argv[i]);
                std::exit(2);
            }
            value = argv[++i];
        };
        if (!std::strcmp(argv[i], "--problem")) next_int(opt.problem);
        else if (!std::strcmp(argv[i], "--run")) next_int(opt.run);
        else if (!std::strcmp(argv[i], "--seed")) next_int(opt.seed);
        else if (!std::strcmp(argv[i], "--N")) next_int(opt.popsize);
        else if (!std::strcmp(argv[i], "--D")) next_int(opt.dim);
        else if (!std::strcmp(argv[i], "--maxFE")) next_int(opt.maxFE);
        else if (!std::strcmp(argv[i], "--save")) next_int(opt.save);
        else if (!std::strcmp(argv[i], "--pfRoot")) next_string(opt.pfRoot);
        else if (!std::strcmp(argv[i], "--noIGD")) opt.igd = false;
        else if (!std::strcmp(argv[i], "--toprace")) opt.toprace = true;
        else if (!std::strcmp(argv[i], "--noTopRace")) opt.toprace = false;
        else if (!std::strcmp(argv[i], "--stateTrace")) opt.stateTrace = true;
        else if (!std::strcmp(argv[i], "--stateTraceLimit")) next_int(opt.stateTraceLimit);
        else if (!std::strcmp(argv[i], "--help"))
        {
            std::printf("Usage: %s --problem 1 --run 1 --seed 20261611 --N 100 --maxFE 200000 --pfRoot references/cmop_sdc_pf\n", argv[0]);
            std::exit(0);
        }
        else
        {
            std::fprintf(stderr, "Unknown option: %s\n", argv[i]);
            std::exit(2);
        }
    }
    if (opt.problem < 1 || opt.problem > 15)
    {
        std::fprintf(stderr, "problem must be in 1..15\n");
        std::exit(2);
    }
    return opt;
}

static std::string default_pf_root()
{
    const char *env = std::getenv("RDE_CPP_PF_ROOT");
    if (env && *env) return std::string(env);
    return "references/cmop_sdc_pf";
}

static ReferenceFront load_reference_front(int problem, int nobj, const std::string &root)
{
    ReferenceFront ref;
    ref.nobj = nobj;
    char fileName[256];
    std::snprintf(fileName, sizeof(fileName), "F%d_M%d_N10000.txt", problem, nobj);
    std::string path = root;
    if (!path.empty() && path[path.size() - 1] != '/') path += "/";
    path += fileName;

    std::ifstream in(path.c_str());
    if (!in)
    {
        std::fprintf(stderr, "WARNING reference front not found: %s\n", path.c_str());
        return ref;
    }
    std::vector<double> values;
    double x = 0.0;
    while (in >> x) values.push_back(x);
    if (values.empty() || values.size() % (size_t)nobj != 0)
    {
        std::fprintf(stderr, "WARNING invalid reference front: %s values=%zu nobj=%d\n",
                     path.c_str(), values.size(), nobj);
        return ref;
    }
    ref.n = (int)(values.size() / (size_t)nobj);
    ref.obj.swap(values);
    ref.available = true;
    return ref;
}

static double constraint_violation(const Population &pop, int i)
{
    double cv = 0.0;
    const size_t base = (size_t)i * pop.ncon;
    for (int j = 0; j < pop.ncon; ++j)
    {
        if (pop.con[base + j] > 0.0) cv += pop.con[base + j];
    }
    return cv;
}

static bool dominates_by_cv_obj(const Population &pop, int a, int b, double epsilon)
{
    double cva = constraint_violation(pop, a);
    double cvb = constraint_violation(pop, b);
    if (cva <= epsilon) cva = 0.0;
    if (cvb <= epsilon) cvb = 0.0;
    if (cva < cvb) return true;
    if (cva > cvb) return false;

    bool less = false;
    bool greater = false;
    for (int m = 0; m < pop.nobj; ++m)
    {
        double oa = pop.obj[(size_t)a * pop.nobj + m];
        double ob = pop.obj[(size_t)b * pop.nobj + m];
        if (oa < ob) less = true;
        else if (oa > ob) greater = true;
    }
    return less && !greater;
}

static bool feasible_all_constraints(const Population &pop, int i)
{
    const size_t base = (size_t)i * pop.ncon;
    for (int j = 0; j < pop.ncon; ++j)
    {
        if (pop.con[base + j] > 0.0) return false;
    }
    return true;
}

static bool dominates_objectives(const Population &pop, int a, int b)
{
    bool strictlyBetter = false;
    for (int m = 0; m < pop.nobj; ++m)
    {
        const double oa = pop.obj[(size_t)a * pop.nobj + m];
        const double ob = pop.obj[(size_t)b * pop.nobj + m];
        if (oa > ob) return false;
        if (oa < ob) strictlyBetter = true;
    }
    return strictlyBetter;
}

static std::vector<int> feasible_nondominated_indices(const Population &pop)
{
    std::vector<int> feasible;
    feasible.reserve(pop.n);
    for (int i = 0; i < pop.n; ++i)
    {
        if (feasible_all_constraints(pop, i)) feasible.push_back(i);
    }

    std::vector<int> firstFront;
    firstFront.reserve(feasible.size());
    for (int idx : feasible)
    {
        bool dominated = false;
        for (int other : feasible)
        {
            if (other != idx && dominates_objectives(pop, other, idx))
            {
                dominated = true;
                break;
            }
        }
        if (!dominated) firstFront.push_back(idx);
    }
    return firstFront;
}

static double igd_value(const Population &pop, const ReferenceFront &ref)
{
    if (!ref.available || ref.nobj != pop.nobj)
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
    std::vector<int> front = feasible_nondominated_indices(pop);
    if (front.empty())
    {
        return std::numeric_limits<double>::quiet_NaN();
    }

    double sum = 0.0;
    for (int r = 0; r < ref.n; ++r)
    {
        double bestDist2 = std::numeric_limits<double>::infinity();
        const size_t refBase = (size_t)r * ref.nobj;
        for (int idx : front)
        {
            double dist2 = 0.0;
            const size_t popBase = (size_t)idx * pop.nobj;
            for (int m = 0; m < pop.nobj; ++m)
            {
                const double diff = ref.obj[refBase + m] - pop.obj[popBase + m];
                dist2 += diff * diff;
            }
            if (dist2 < bestDist2) bestDist2 = dist2;
        }
        sum += std::sqrt(bestDist2);
    }
    return sum / (double)ref.n;
}

static std::vector<double> pairwise_density(const Population &pop, const std::vector<int> &idx)
{
    const int n = (int)idx.size();
    const int kth = std::max(0, (int)std::floor(std::sqrt((double)n)) - 1);
    std::vector<double> d(n, 0.0);
    std::vector<double> distances(n);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            if (i == j)
            {
                distances[j] = std::numeric_limits<double>::infinity();
                continue;
            }
            double s = 0.0;
            for (int m = 0; m < pop.nobj; ++m)
            {
                const double diff = pop.obj[(size_t)idx[i] * pop.nobj + m] - pop.obj[(size_t)idx[j] * pop.nobj + m];
                s += diff * diff;
            }
            distances[j] = std::sqrt(s);
        }
        std::nth_element(distances.begin(), distances.begin() + kth, distances.end());
        d[i] = 1.0 / (distances[kth] + 2.0);
    }
    return d;
}

static std::vector<double> pairwise_density_matrix(const std::vector<double> &obj, int n, int nobj)
{
    const int kth = std::max(0, (int)std::floor(std::sqrt((double)n)) - 1);
    std::vector<double> d(n, 0.0);
    std::vector<double> distances(n);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            if (i == j)
            {
                distances[j] = std::numeric_limits<double>::infinity();
                continue;
            }
            double s = 0.0;
            for (int m = 0; m < nobj; ++m)
            {
                const double diff = obj[(size_t)i * nobj + m] - obj[(size_t)j * nobj + m];
                s += diff * diff;
            }
            distances[j] = std::sqrt(s);
        }
        std::sort(distances.begin(), distances.end());
        d[i] = 1.0 / (distances[(size_t)kth] + 2.0);
    }
    return d;
}

static bool dominates_matrix(const std::vector<double> &obj, const std::vector<double> &cv,
                             int nobj, int a, int b)
{
    if (cv[a] < cv[b]) return true;
    if (cv[a] > cv[b]) return false;
    bool anyLess = false;
    bool anyGreater = false;
    for (int m = 0; m < nobj; ++m)
    {
        const double oa = obj[(size_t)a * nobj + m];
        const double ob = obj[(size_t)b * nobj + m];
        if (oa < ob) anyLess = true;
        else if (oa > ob) anyGreater = true;
    }
    return anyLess && !anyGreater;
}

static std::vector<double> cal_fitness_matrix(const std::vector<double> &obj, int n, int nobj,
                                              const std::vector<double> *cvInput,
                                              double epsilon, bool applyEpsilon)
{
    if (n <= 0) return std::vector<double>();
    std::vector<unsigned char> dom((size_t)n * n, 0);
    std::vector<double> s(n, 0.0), r(n, 0.0);
    std::vector<double> cv(n, 0.0);
    if (cvInput)
    {
        cv = *cvInput;
        if (applyEpsilon)
        {
            for (double &v : cv) if (v <= epsilon) v = 0.0;
        }
    }

    for (int i = 0; i < n - 1; ++i)
    {
        for (int j = i + 1; j < n; ++j)
        {
            if (dominates_matrix(obj, cv, nobj, i, j))
            {
                dom[(size_t)i * n + j] = 1;
            }
            else if (dominates_matrix(obj, cv, nobj, j, i))
            {
                dom[(size_t)j * n + i] = 1;
            }
        }
    }
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j) s[i] += dom[(size_t)i * n + j] ? 1.0 : 0.0;
    }
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            if (dom[(size_t)j * n + i]) r[i] += s[j];
        }
    }
    std::vector<double> density = pairwise_density_matrix(obj, n, nobj);
    std::vector<double> fit(n);
    for (int i = 0; i < n; ++i) fit[i] = r[i] + density[i];
    return fit;
}

static std::vector<double> cal_fitness(const Population &pop, const std::vector<int> &idx, double epsilon)
{
    const int n = (int)idx.size();
    std::vector<double> obj((size_t)n * pop.nobj);
    std::vector<double> cv(n, 0.0);
    for (int i = 0; i < n; ++i)
    {
        int src = idx[(size_t)i];
        std::copy_n(&pop.obj[(size_t)src * pop.nobj], pop.nobj, &obj[(size_t)i * pop.nobj]);
        cv[(size_t)i] = constraint_violation(pop, src);
    }
    return cal_fitness_matrix(obj, n, pop.nobj, &cv, epsilon, epsilon >= 0.0);
}

static std::vector<double> cal_fitness_no_epsilon(const Population &pop, const std::vector<int> &idx)
{
    return cal_fitness(pop, idx, -1.0);
}

static std::vector<double> cal_fitness_objectives_only(const std::vector<double> &obj, int n, int nobj)
{
    return cal_fitness_matrix(obj, n, nobj, nullptr, 0.0, false);
}

static void evaluate(Population &pop, int func_num)
{
    static const double boundarySnapTol = []() {
        const char *raw = std::getenv("RDE_CPP_BOUNDARY_SNAP_TOL");
        return (raw && *raw) ? std::max(0.0, std::atof(raw)) : 0.0;
    }();
    for (double &x : pop.dec)
    {
        if (x < 0.0) x = 0.0;
        else if (x > 1.0) x = 1.0;
        else if (std::fabs(x) <= boundarySnapTol) x = 0.0;
        else if (std::fabs(x - 1.0) <= boundarySnapTol) x = 1.0;
    }
    cmop_cec2026_test(pop.dec.data(), pop.obj.data(), pop.con.data(), nullptr,
                      pop.dim, pop.n, func_num);
}

static void initialize(Population &pop, int func_num, MatlabRng &rng)
{
    for (int d = 0; d < pop.dim; ++d)
    {
        for (int i = 0; i < pop.n; ++i)
        {
            pop.dec[(size_t)i * pop.dim + d] = rng.rand();
        }
    }
    evaluate(pop, func_num);
}

static int sample_from_pool(const std::vector<int> &pool, MatlabRng &rng)
{
    if (pool.empty())
    {
        std::fprintf(stderr, "empty random pool\n");
        std::exit(2);
    }
    return pool[(size_t)rng.randi((int)pool.size()) - 1];
}

static std::vector<int> setdiff_one_based(int n, const std::vector<int> &excludeZeroBased)
{
    std::vector<int> pool;
    pool.reserve(n);
    for (int i = 0; i < n; ++i)
    {
        bool blocked = false;
        for (int e : excludeZeroBased)
        {
            if (i == e)
            {
                blocked = true;
                break;
            }
        }
        if (!blocked) pool.push_back(i);
    }
    return pool;
}

struct OffspringResult
{
    Population pop;
    int evaluations = 0;
    OffspringResult(const Population &p, int e) : pop(p), evaluations(e) {}
};

struct TopRaceParams
{
    bool enabled = true;
    double targetTopFrac = 0.05;
    double sampleTopFrac = 0.20;
    double progressMin = 0.15;
    double progressMax = 0.75;
    double feasRateMin = 0.50;
    double altCrMin = 0.26;
    double altCrMax = 0.85;
    double pbestPullScale = 0.62;
    double diffScale = 1.16;
    int pbestRetry = 7;
    bool useFScale = false;
    std::string crSource = "random";
    std::string boundRepair = "clip";
    double boundRepairMix = 0.5;
    std::string insertMode = "replace";
};

struct ArchiveParams
{
    bool enabled = false;
    bool output = false;
    bool reinject = false;
    double capFrac = 3.0;
    double outputProgress = 0.82;
    double reinjectFrac = 0.08;
    double reinjectProgressMin = 0.45;
    double reinjectInterval = 0.12;
};

static Population merge_pop(const Population &a, const Population &b);

static double env_number(const char *name, double defaultValue)
{
    const char *raw = std::getenv(name);
    if (!raw || !*raw) return defaultValue;
    return std::atof(raw);
}

static bool env_flag(const char *name, bool defaultValue)
{
    const char *raw = std::getenv(name);
    if (!raw || !*raw) return defaultValue;
    std::string value(raw);
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) { return (char)std::tolower(c); });
    return value == "1" || value == "true" || value == "yes" || value == "y" || value == "on";
}

static std::string env_string(const char *name, const char *defaultValue)
{
    const char *raw = std::getenv(name);
    std::string value = raw && *raw ? std::string(raw) : std::string(defaultValue);
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) { return (char)std::tolower(c); });
    return value;
}

static TopRaceParams default_toprace_params()
{
    TopRaceParams p;
    p.targetTopFrac = 0.08;
    p.sampleTopFrac = 0.27;
    p.progressMin = 0.0;
    p.progressMax = 0.76;
    p.feasRateMin = 0.48;
    p.altCrMin = 0.26;
    p.altCrMax = 0.85;
    p.pbestPullScale = 0.62;
    p.diffScale = 1.16;
    p.pbestRetry = 7;
    p.useFScale = true;
    p.crSource = "base";
    p.boundRepair = "mix";
    p.boundRepairMix = 0.21;
    p.insertMode = "replace";
    return p;
}

static TopRaceParams read_toprace_params()
{
    TopRaceParams p = default_toprace_params();
    p.enabled = env_flag("TOPRACE_ENABLED", true);
    p.targetTopFrac = env_number("TOPRACE_TARGET_TOP_FRAC", p.targetTopFrac);
    p.sampleTopFrac = env_number("TOPRACE_SAMPLE_TOP_FRAC", p.sampleTopFrac);
    p.progressMin = env_number("TOPRACE_PROGRESS_MIN", p.progressMin);
    p.progressMax = env_number("TOPRACE_PROGRESS_MAX", p.progressMax);
    p.feasRateMin = env_number("TOPRACE_FEAS_RATE_MIN", p.feasRateMin);
    p.altCrMin = env_number("TOPRACE_ALT_CR_MIN", p.altCrMin);
    p.altCrMax = env_number("TOPRACE_ALT_CR_MAX", p.altCrMax);
    p.pbestPullScale = env_number("TOPRACE_PBEST_PULL_SCALE", p.pbestPullScale);
    p.diffScale = env_number("TOPRACE_DIFF_SCALE", p.diffScale);
    p.pbestRetry = std::max(1, (int)std::llround(env_number("TOPRACE_PBEST_RETRY", (double)p.pbestRetry)));
    p.useFScale = env_flag("TOPRACE_USE_F_SCALE", p.useFScale);
    p.crSource = env_string("TOPRACE_CR_SOURCE", p.crSource.c_str());
    p.boundRepair = env_string("TOPRACE_BOUND_REPAIR", p.boundRepair.c_str());
    p.boundRepairMix = env_number("TOPRACE_BOUND_REPAIR_MIX", p.boundRepairMix);
    p.insertMode = env_string("TOPRACE_INSERT_MODE", p.insertMode.c_str());
    return p;
}

static ArchiveParams read_archive_params()
{
    ArchiveParams p;
    p.enabled = env_flag("CMOP_ARCHIVE_ENABLED", p.enabled);
    p.output = env_flag("CMOP_ARCHIVE_OUTPUT", p.output);
    p.reinject = env_flag("CMOP_ARCHIVE_REINJECT", p.reinject);
    p.capFrac = env_number("CMOP_ARCHIVE_CAP_FRAC", p.capFrac);
    p.outputProgress = env_number("CMOP_ARCHIVE_OUTPUT_PROGRESS", p.outputProgress);
    p.reinjectFrac = env_number("CMOP_ARCHIVE_REINJECT_FRAC", p.reinjectFrac);
    p.reinjectProgressMin = env_number("CMOP_ARCHIVE_REINJECT_PROGRESS_MIN", p.reinjectProgressMin);
    p.reinjectInterval = env_number("CMOP_ARCHIVE_REINJECT_INTERVAL", p.reinjectInterval);
    if (p.output || p.reinject) p.enabled = true;
    if (p.capFrac < 0.0) p.capFrac = 0.0;
    if (p.outputProgress < 0.0) p.outputProgress = 0.0;
    if (p.outputProgress > 1.0) p.outputProgress = 1.0;
    if (p.reinjectFrac < 0.0) p.reinjectFrac = 0.0;
    if (p.reinjectProgressMin < 0.0) p.reinjectProgressMin = 0.0;
    if (p.reinjectProgressMin > 1.0) p.reinjectProgressMin = 1.0;
    if (p.reinjectInterval <= 0.0) p.reinjectInterval = 0.12;
    return p;
}

static double repair_toprace_bound(double value, double current, const TopRaceParams &params)
{
    if (params.boundRepair == "midpoint")
    {
        if (value < 0.0) value = 0.5 * current;
        if (value > 1.0) value = 0.5 * (current + 1.0);
    }
    else if (params.boundRepair == "mix")
    {
        const double mix = std::min(std::max(params.boundRepairMix, 0.0), 1.0);
        if (value < 0.0) value = (1.0 - mix) * current;
        if (value > 1.0) value = (1.0 - mix) * current + mix;
    }
    if (value < 0.0) value = 0.0;
    if (value > 1.0) value = 1.0;
    return value;
}

static bool toprace_better_by_cmop(const Population &alt, int altIndex,
                                   const Population &base, int baseIndex,
                                   const Population &context, double epsilon)
{
    const double altCV = constraint_violation(alt, altIndex);
    const double baseCV = constraint_violation(base, baseIndex);
    const bool altGood = altCV <= epsilon;
    const bool baseGood = baseCV <= epsilon;
    if (altGood && !baseGood) return true;
    if (!altGood && baseGood) return false;
    if (!altGood && !baseGood && altCV != baseCV) return altCV < baseCV;

    bool altDominates = true, altStrict = false;
    bool baseDominates = true, baseStrict = false;
    for (int m = 0; m < context.nobj; ++m)
    {
        const double ao = alt.obj[(size_t)altIndex * alt.nobj + m];
        const double bo = base.obj[(size_t)baseIndex * base.nobj + m];
        if (ao > bo) altDominates = false;
        if (ao < bo) altStrict = true;
        if (bo > ao) baseDominates = false;
        if (bo < ao) baseStrict = true;
    }
    altDominates = altDominates && altStrict;
    baseDominates = baseDominates && baseStrict;
    if (altDominates && !baseDominates) return true;
    if (baseDominates && !altDominates) return false;

    const int n = context.n + 2;
    std::vector<double> obj((size_t)n * context.nobj);
    std::vector<double> cv(n, 0.0);
    for (int i = 0; i < context.n; ++i)
    {
        std::copy_n(&context.obj[(size_t)i * context.nobj], context.nobj, &obj[(size_t)i * context.nobj]);
        cv[(size_t)i] = constraint_violation(context, i);
    }
    std::copy_n(&base.obj[(size_t)baseIndex * base.nobj], base.nobj, &obj[(size_t)context.n * context.nobj]);
    cv[(size_t)context.n] = baseCV;
    std::copy_n(&alt.obj[(size_t)altIndex * alt.nobj], alt.nobj, &obj[(size_t)(context.n + 1) * context.nobj]);
    cv[(size_t)context.n + 1] = altCV;
    std::vector<double> fit = cal_fitness_matrix(obj, n, context.nobj, &cv, epsilon, true);
    return fit[(size_t)n - 1] < fit[(size_t)n - 2];
}

static OffspringResult make_offspring(const Population &pop, const std::vector<double> &fitness,
                                      const Options &opt, int max_trials, int fe, int maxFE,
                                      double epsilon, const TopRaceParams &topraceParams,
                                      MatlabRng &rng)
{
    static const double Fm[3] = {0.6, 0.8, 1.0};
    static const double CRm[3] = {0.1, 0.2, 1.0};
    const int trialCount = std::min(pop.n, std::max(0, max_trials));
    Population off(trialCount, pop.dim, pop.nobj, pop.ncon);
    if (trialCount <= 0) return OffspringResult(off, 0);

    std::vector<int> sorted(pop.n);
    std::iota(sorted.begin(), sorted.end(), 0);
    std::stable_sort(sorted.begin(), sorted.end(), [&](int a, int b) {
        return fitness[a] < fitness[b];
    });

    const double progress = (double)fe / (double)maxFE;
    const int pNP = std::max((int)std::floor(pop.n * (1.0 - 0.99 * progress)), 2);

    std::vector<int> r1(trialCount), r2(trialCount);
    for (int i = 0; i < trialCount; ++i)
    {
        r1[(size_t)i] = sample_from_pool(setdiff_one_based(pop.n, {i}), rng);
    }
    for (int i = 0; i < trialCount; ++i)
    {
        r2[(size_t)i] = sample_from_pool(setdiff_one_based(pop.n, {i, r1[(size_t)i]}), rng);
    }

    std::vector<int> pbest(trialCount);
    for (int i = 0; i < trialCount; ++i)
    {
        int pbestPos = (int)std::ceil(rng.rand() * (double)pNP) - 1;
        if (pbestPos < 0) pbestPos = 0;
        if (pbestPos >= pNP) pbestPos = pNP - 1;
        pbest[(size_t)i] = sorted[(size_t)pbestPos];
    }

    std::vector<double> baseF(trialCount), baseCr(trialCount), signedF2(trialCount);
    for (int i = 0; i < trialCount; ++i)
    {
        baseF[(size_t)i] = Fm[(size_t)rng.randi(3) - 1];
        signedF2[(size_t)i] = fitness[(size_t)r1[(size_t)i]] > fitness[(size_t)r2[(size_t)i]] ? -baseF[(size_t)i] : baseF[(size_t)i];
    }
    for (int i = 0; i < trialCount; ++i)
    {
        baseCr[(size_t)i] = CRm[(size_t)rng.randi(3) - 1];
    }

    std::vector<unsigned char> site((size_t)trialCount * pop.dim, 0);
    for (int d = 0; d < pop.dim; ++d)
    {
        for (int i = 0; i < trialCount; ++i)
        {
            site[(size_t)i * pop.dim + d] = rng.rand() < baseCr[(size_t)i] ? 1 : 0;
        }
    }
    const bool perturbation = rng.rand() < 0.2;

    for (int i = 0; i < trialCount; ++i)
    {
        for (int d = 0; d < pop.dim; ++d)
        {
            double value = pop.dec[(size_t)i * pop.dim + d];
            if (site[(size_t)i * pop.dim + d])
            {
                value = pop.dec[(size_t)i * pop.dim + d]
                    + baseF[(size_t)i] * (pop.dec[(size_t)pbest[(size_t)i] * pop.dim + d] - pop.dec[(size_t)i * pop.dim + d])
                    + signedF2[(size_t)i] * (pop.dec[(size_t)r1[(size_t)i] * pop.dim + d] - pop.dec[(size_t)r2[(size_t)i] * pop.dim + d]);
                if (value < 0.0) value = 0.0;
                if (value > 1.0) value = 1.0;
            }
            off.dec[(size_t)i * off.dim + d] = value;
        }
    }
    if (perturbation)
    {
        const double u = rng.rand();
        const double delta = 0.1 * std::tan(3.14159265358979323846 * (u - 0.5));
        for (int i = 0; i < trialCount; ++i)
        {
            for (int d = 0; d < pop.dim; ++d)
            {
                if (!site[(size_t)i * pop.dim + d])
                {
                    off.dec[(size_t)i * off.dim + d] = pop.dec[(size_t)i * pop.dim + d] + delta;
                }
            }
        }
    }
    for (double &value : off.dec)
    {
        if (value < 0.0) value = 0.0;
        if (value > 1.0) value = 1.0;
    }
    evaluate(off, opt.problem);
    int evaluations = trialCount;

    if (opt.toprace)
    {
        TopRaceParams params = topraceParams;
        double feasRate = 0.0;
        for (int i = 0; i < pop.n; ++i)
        {
            if (constraint_violation(pop, i) <= epsilon) feasRate += 1.0;
        }
        feasRate /= (double)pop.n;
        const int remainingFE = std::max(0, maxFE - (fe + evaluations));
        const int targetCount = std::min(std::min((int)std::ceil(pop.n * params.targetTopFrac), pop.n), remainingFE);
        if (params.enabled && trialCount == pop.n && targetCount > 0 &&
            progress >= params.progressMin && progress < params.progressMax &&
            feasRate > params.feasRateMin)
        {
            const int sampleLimit = std::min(pop.n, std::max(2, (int)std::ceil(pop.n * params.sampleTopFrac)));
            std::vector<int> targetIndex(targetCount);
            for (int i = 0; i < targetCount; ++i) targetIndex[(size_t)i] = sorted[(size_t)i];
            Population alt(targetCount, pop.dim, pop.nobj, pop.ncon);
            for (int i = 0; i < targetCount; ++i)
            {
                const int target = targetIndex[(size_t)i];
                int pbestIndex = target;
                for (int retry = 0; retry < params.pbestRetry; ++retry)
                {
                    int candidate = sorted[(size_t)rng.randi(sampleLimit) - 1];
                    if (candidate != target)
                    {
                        pbestIndex = candidate;
                        break;
                    }
                }
                if (pbestIndex == target)
                {
                    pbestIndex = sample_from_pool(setdiff_one_based(pop.n, {target}), rng);
                }
                int randA = sample_from_pool(setdiff_one_based(pop.n, {target, pbestIndex}), rng);
                int randC = sample_from_pool(setdiff_one_based(pop.n, {target, pbestIndex, randA}), rng);
                double cr = 0.0;
                if (params.crSource == "base")
                {
                    cr = baseCr[(size_t)target];
                }
                else
                {
                    cr = CRm[(size_t)rng.randi(3) - 1];
                }
                cr = std::min(std::max(cr, params.altCrMin), params.altCrMax);
                std::vector<unsigned char> altSite(pop.dim, 0);
                for (int d = 0; d < pop.dim; ++d) altSite[(size_t)d] = rng.rand() < cr ? 1 : 0;
                int jrand = rng.randi(pop.dim) - 1;
                altSite[(size_t)jrand] = 1;
                double stepScale = params.useFScale ? baseF[(size_t)target] : 1.0;
                for (int d = 0; d < pop.dim; ++d)
                {
                    double value = pop.dec[(size_t)target * pop.dim + d];
                    if (altSite[(size_t)d])
                    {
                        value = pop.dec[(size_t)target * pop.dim + d]
                            + params.pbestPullScale * stepScale * (pop.dec[(size_t)pbestIndex * pop.dim + d] - pop.dec[(size_t)target * pop.dim + d])
                            + params.diffScale * stepScale * (pop.dec[(size_t)randA * pop.dim + d] - pop.dec[(size_t)randC * pop.dim + d]);
                    }
                    value = repair_toprace_bound(value, pop.dec[(size_t)target * pop.dim + d], params);
                    alt.dec[(size_t)i * alt.dim + d] = value;
                }
            }
            evaluate(alt, opt.problem);
            evaluations += targetCount;
            if (params.insertMode == "append_all")
            {
                off = merge_pop(off, alt);
            }
            else
            {
                std::vector<int> accepted;
                for (int i = 0; i < targetCount; ++i)
                {
                    if (toprace_better_by_cmop(alt, i, off, targetIndex[(size_t)i], pop, epsilon))
                    {
                        accepted.push_back(i);
                    }
                }
                if (params.insertMode == "append_better")
                {
                    Population selected((int)accepted.size(), pop.dim, pop.nobj, pop.ncon);
                    for (int i = 0; i < selected.n; ++i)
                    {
                        int src = accepted[(size_t)i];
                        std::copy_n(&alt.dec[(size_t)src * alt.dim], alt.dim, &selected.dec[(size_t)i * selected.dim]);
                        std::copy_n(&alt.obj[(size_t)src * alt.nobj], alt.nobj, &selected.obj[(size_t)i * selected.nobj]);
                        std::copy_n(&alt.con[(size_t)src * alt.ncon], alt.ncon, &selected.con[(size_t)i * selected.ncon]);
                    }
                    off = merge_pop(off, selected);
                }
                else
                {
                    for (int src : accepted)
                    {
                        int dst = targetIndex[(size_t)src];
                        std::copy_n(&alt.dec[(size_t)src * alt.dim], alt.dim, &off.dec[(size_t)dst * off.dim]);
                        std::copy_n(&alt.obj[(size_t)src * alt.nobj], alt.nobj, &off.obj[(size_t)dst * off.nobj]);
                        std::copy_n(&alt.con[(size_t)src * alt.ncon], alt.ncon, &off.con[(size_t)dst * off.ncon]);
                    }
                }
            }
        }
    }

    return OffspringResult(off, evaluations);
}

static Population merge_pop(const Population &a, const Population &b)
{
    Population c(a.n + b.n, a.dim, a.nobj, a.ncon);
    std::copy(a.dec.begin(), a.dec.end(), c.dec.begin());
    std::copy(b.dec.begin(), b.dec.end(), c.dec.begin() + a.dec.size());
    std::copy(a.obj.begin(), a.obj.end(), c.obj.begin());
    std::copy(b.obj.begin(), b.obj.end(), c.obj.begin() + a.obj.size());
    std::copy(a.con.begin(), a.con.end(), c.con.begin());
    std::copy(b.con.begin(), b.con.end(), c.con.begin() + a.con.size());
    return c;
}

static std::vector<int> truncation(const Population &pop, const std::vector<int> &idx, int removeCount)
{
    const int n = (int)idx.size();
    std::vector<double> dist((size_t)n * n);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            if (i == j)
            {
                dist[(size_t)i * n + j] = std::numeric_limits<double>::infinity();
            }
            else
            {
                double s = 0.0;
                for (int m = 0; m < pop.nobj; ++m)
                {
                    double diff = pop.obj[(size_t)idx[i] * pop.nobj + m] - pop.obj[(size_t)idx[j] * pop.nobj + m];
                    s += diff * diff;
                }
                dist[(size_t)i * n + j] = std::sqrt(s);
            }
        }
    }
    std::vector<unsigned char> del(n, 0);
    for (int removed = 0; removed < removeCount; ++removed)
    {
        int best = -1;
        std::vector<double> bestRow;
        for (int i = 0; i < n; ++i)
        {
            if (del[i]) continue;
            std::vector<double> row;
            row.reserve(n);
            for (int j = 0; j < n; ++j)
            {
                if (!del[j]) row.push_back(dist[(size_t)i * n + j]);
            }
            std::sort(row.begin(), row.end());
            if (best < 0 || std::lexicographical_compare(row.begin(), row.end(), bestRow.begin(), bestRow.end()))
            {
                best = i;
                bestRow.swap(row);
            }
        }
        if (best >= 0) del[best] = 1;
    }
    std::vector<int> keep;
    keep.reserve(n - removeCount);
    for (int i = 0; i < n; ++i)
    {
        if (!del[i]) keep.push_back(idx[i]);
    }
    return keep;
}

static Population select_by_indices(const Population &pop, const std::vector<int> &idx)
{
    Population out((int)idx.size(), pop.dim, pop.nobj, pop.ncon);
    for (int i = 0; i < out.n; ++i)
    {
        int src = idx[i];
        std::copy_n(&pop.dec[(size_t)src * pop.dim], pop.dim, &out.dec[(size_t)i * out.dim]);
        std::copy_n(&pop.obj[(size_t)src * pop.nobj], pop.nobj, &out.obj[(size_t)i * out.nobj]);
        std::copy_n(&pop.con[(size_t)src * pop.ncon], pop.ncon, &out.con[(size_t)i * out.ncon]);
    }
    return out;
}

static std::vector<double> feasible_augmented_fitness(const Population &pop, const std::vector<int> &idx)
{
    const int n = (int)idx.size();
    std::vector<double> obj((size_t)n * (pop.nobj + 1));
    for (int i = 0; i < n; ++i)
    {
        const int src = idx[(size_t)i];
        for (int m = 0; m < pop.nobj; ++m)
        {
            obj[(size_t)i * (pop.nobj + 1) + m] = pop.obj[(size_t)src * pop.nobj + m];
        }
        obj[(size_t)i * (pop.nobj + 1) + pop.nobj] = constraint_violation(pop, src);
    }
    return cal_fitness_objectives_only(obj, n, pop.nobj + 1);
}

static void environmental_selection(Population &pop, int targetN, double epsilon)
{
    std::vector<int> feas, infeas;
    feas.reserve(pop.n);
    infeas.reserve(pop.n);
    for (int i = 0; i < pop.n; ++i)
    {
        if (constraint_violation(pop, i) <= epsilon) feas.push_back(i);
        else infeas.push_back(i);
    }

    std::vector<int> chosen;
    std::vector<double> chosenFitness;
    auto append_sorted = [&](const std::vector<int> &source, std::vector<double> fit, int limit, double offset) {
        if (source.empty() || limit <= 0) return;
        std::vector<unsigned char> next(source.size(), 0);
        int count = 0;
        for (size_t i = 0; i < source.size(); ++i)
        {
            if (fit[i] < 1.0)
            {
                next[i] = 1;
                ++count;
            }
        }
        std::vector<int> order(source.size());
        std::iota(order.begin(), order.end(), 0);
        std::stable_sort(order.begin(), order.end(), [&](int a, int b) { return fit[a] < fit[b]; });
        for (int k = 0; count < limit && k < (int)order.size(); ++k)
        {
            if (!next[order[k]])
            {
                next[order[k]] = 1;
                ++count;
            }
        }
        if (count > limit)
        {
            std::vector<int> nextIdx;
            for (size_t i = 0; i < source.size(); ++i) if (next[i]) nextIdx.push_back(source[i]);
            std::vector<int> keep = truncation(pop, nextIdx, count - limit);
            std::fill(next.begin(), next.end(), 0);
            for (int id : keep)
            {
                auto it = std::find(source.begin(), source.end(), id);
                if (it != source.end()) next[(size_t)(it - source.begin())] = 1;
            }
        }
        std::vector<int> selectedLocal;
        for (size_t i = 0; i < source.size(); ++i) if (next[i]) selectedLocal.push_back((int)i);
        std::stable_sort(selectedLocal.begin(), selectedLocal.end(), [&](int a, int b) { return fit[a] < fit[b]; });
        for (int li : selectedLocal)
        {
            chosen.push_back(source[(size_t)li]);
            chosenFitness.push_back(fit[(size_t)li] + offset);
        }
    };

    if (feas.empty())
    {
        std::vector<double> fit = cal_fitness_no_epsilon(pop, infeas);
        append_sorted(infeas, fit, targetN, 0.0);
    }
    else if ((int)feas.size() <= targetN)
    {
        std::vector<double> fitF = feasible_augmented_fitness(pop, feas);
        append_sorted(feas, fitF, (int)feas.size(), 0.0);
        double offset = chosenFitness.empty() ? 0.0 : *std::max_element(chosenFitness.begin(), chosenFitness.end());
        std::vector<double> fitI = cal_fitness_no_epsilon(pop, infeas);
        append_sorted(infeas, fitI, targetN - (int)chosen.size(), offset);
    }
    else
    {
        std::vector<double> fitF = feasible_augmented_fitness(pop, feas);
        append_sorted(feas, fitF, targetN, 0.0);
    }

    pop = select_by_indices(pop, chosen);
    pop.fitness = chosenFitness;
    if ((int)pop.fitness.size() != pop.n)
    {
        std::vector<int> all(pop.n);
        std::iota(all.begin(), all.end(), 0);
        pop.fitness = cal_fitness(pop, all, epsilon);
    }
}


static void push_solution(Population &out, const Population &src, int idx)
{
    for (int d = 0; d < src.dim; ++d) out.dec.push_back(src.dec[(size_t)idx * src.dim + d]);
    for (int m = 0; m < src.nobj; ++m) out.obj.push_back(src.obj[(size_t)idx * src.nobj + m]);
    for (int j = 0; j < src.ncon; ++j) out.con.push_back(src.con[(size_t)idx * src.ncon + j]);
    out.fitness.push_back(0.0);
    ++out.n;
}

static bool same_objectives_between(const Population &a, int ia, const Population &b, int ib)
{
    for (int m = 0; m < a.nobj; ++m)
    {
        const double x = a.obj[(size_t)ia * a.nobj + m];
        const double y = b.obj[(size_t)ib * b.nobj + m];
        const double tol = 1.0e-12 * std::max(1.0, std::max(std::fabs(x), std::fabs(y)));
        if (std::fabs(x - y) > tol) return false;
    }
    return true;
}

static bool dominates_objectives_between(const Population &a, int ia, const Population &b, int ib)
{
    bool strictlyBetter = false;
    for (int m = 0; m < a.nobj; ++m)
    {
        const double x = a.obj[(size_t)ia * a.nobj + m];
        const double y = b.obj[(size_t)ib * b.nobj + m];
        if (x > y) return false;
        if (x < y) strictlyBetter = true;
    }
    return strictlyBetter;
}

static std::vector<double> normalize01_values(const std::vector<double> &values)
{
    std::vector<double> out(values.size(), 0.0);
    if (values.empty()) return out;
    const double lo = *std::min_element(values.begin(), values.end());
    const double hi = *std::max_element(values.begin(), values.end());
    const double span = hi - lo;
    if (span <= 1.0e-30) return out;
    for (size_t i = 0; i < values.size(); ++i) out[i] = (values[i] - lo) / span;
    return out;
}

static Population compress_archive_reference(const Population &archive, int capN)
{
    if (archive.n <= capN) return archive;
    if (capN <= 0) return Population(0, archive.dim, archive.nobj, archive.ncon);

    const int n = archive.n;
    const int m = archive.nobj;
    std::vector<double> minObj(m, std::numeric_limits<double>::infinity());
    std::vector<double> maxObj(m, -std::numeric_limits<double>::infinity());
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < m; ++j)
        {
            const double v = archive.obj[(size_t)i * m + j];
            minObj[j] = std::min(minObj[j], v);
            maxObj[j] = std::max(maxObj[j], v);
        }
    }

    std::vector<double> norm((size_t)n * m, 0.0);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < m; ++j)
        {
            const double span = maxObj[j] - minObj[j];
            if (span > 1.0e-30) norm[(size_t)i * m + j] = (archive.obj[(size_t)i * m + j] - minObj[j]) / span;
        }
    }

    std::vector<unsigned char> keep(n, 0);
    for (int j = 0; j < m; ++j)
    {
        int best = 0;
        for (int i = 1; i < n; ++i)
        {
            if (norm[(size_t)i * m + j] < norm[(size_t)best * m + j]) best = i;
        }
        keep[(size_t)best] = 1;
    }

    int kept = 0;
    for (unsigned char v : keep) kept += v ? 1 : 0;
    if (kept < capN)
    {
        std::vector<double> sparse(n, std::numeric_limits<double>::infinity());
        std::vector<double> convergence(n, 0.0);
        for (int i = 0; i < n; ++i)
        {
            double sumNorm = 0.0;
            for (int j = 0; j < m; ++j) sumNorm += norm[(size_t)i * m + j];
            convergence[(size_t)i] = -sumNorm;
            for (int k = 0; k < n; ++k)
            {
                if (i == k) continue;
                double dist2 = 0.0;
                for (int j = 0; j < m; ++j)
                {
                    const double diff = norm[(size_t)i * m + j] - norm[(size_t)k * m + j];
                    dist2 += diff * diff;
                }
                sparse[(size_t)i] = std::min(sparse[(size_t)i], std::sqrt(dist2));
            }
        }
        std::vector<double> sparseScore = normalize01_values(sparse);
        std::vector<double> convScore = normalize01_values(convergence);
        std::vector<int> order;
        order.reserve(n);
        for (int i = 0; i < n; ++i) if (!keep[(size_t)i]) order.push_back(i);
        std::stable_sort(order.begin(), order.end(), [&](int a, int b) {
            const double sa = 0.72 * sparseScore[(size_t)a] + 0.28 * convScore[(size_t)a];
            const double sb = 0.72 * sparseScore[(size_t)b] + 0.28 * convScore[(size_t)b];
            return sa > sb;
        });
        for (int idx : order)
        {
            if (kept >= capN) break;
            keep[(size_t)idx] = 1;
            ++kept;
        }
    }

    Population out(0, archive.dim, archive.nobj, archive.ncon);
    for (int i = 0; i < n; ++i)
    {
        if (keep[(size_t)i] && out.n < capN) push_solution(out, archive, i);
    }
    return out;
}

static Population update_online_archive(Population archive, const Population &cand, int capN)
{
    if (capN <= 0 || cand.n <= 0) return archive;
    for (int i = 0; i < cand.n; ++i)
    {
        if (!feasible_all_constraints(cand, i)) continue;
        bool rejected = false;
        for (int a = 0; a < archive.n; ++a)
        {
            if (same_objectives_between(archive, a, cand, i) || dominates_objectives_between(archive, a, cand, i))
            {
                rejected = true;
                break;
            }
        }
        if (rejected) continue;
        Population kept(0, archive.dim, archive.nobj, archive.ncon);
        for (int a = 0; a < archive.n; ++a)
        {
            if (!dominates_objectives_between(cand, i, archive, a)) push_solution(kept, archive, a);
        }
        push_solution(kept, cand, i);
        archive = std::move(kept);
        if (archive.n > capN) archive = compress_archive_reference(archive, capN);
    }
    return archive;
}

static Population archive_output_population(const Population &mainPop, const Population &archive, int targetN, double progress, const ArchiveParams &params)
{
    if (!params.enabled || !params.output || progress < params.outputProgress || archive.n <= 0) return mainPop;
    Population pool = merge_pop(mainPop, archive);
    std::vector<int> front = feasible_nondominated_indices(pool);
    if (front.empty()) return mainPop;
    if ((int)front.size() > targetN) front = truncation(pool, front, (int)front.size() - targetN);
    return select_by_indices(pool, front);
}

static Population archive_reinject(const Population &archive, const Population &pop, const ArchiveParams &params,
                                   int funcNum, double progress, int remainingFE, MatlabRng &rng)
{
    if (!params.enabled || !params.reinject || archive.n < 2 || remainingFE <= 0) return Population(0, pop.dim, pop.nobj, pop.ncon);
    if (progress < params.reinjectProgressMin) return Population(0, pop.dim, pop.nobj, pop.ncon);
    int q = (int)std::ceil((double)pop.n * params.reinjectFrac);
    q = std::min(std::max(q, 1), remainingFE);
    if (q <= 0) return Population(0, pop.dim, pop.nobj, pop.ncon);

    const int guideCap = std::min(archive.n, std::max(2, std::min(pop.n, q * 2)));
    Population guides = compress_archive_reference(archive, guideCap);
    Population injected(q, pop.dim, pop.nobj, pop.ncon);
    const double localScale = std::max(0.015, 0.08 * (1.0 - progress));
    for (int i = 0; i < q; ++i)
    {
        const int base = rng.randi(pop.n) - 1;
        const int guide = rng.randi(guides.n) - 1;
        const double mix = 0.22 + 0.22 * rng.rand();
        for (int d = 0; d < pop.dim; ++d)
        {
            const double b = pop.dec[(size_t)base * pop.dim + d];
            const double g = guides.dec[(size_t)guide * guides.dim + d];
            const double noise = (2.0 * rng.rand() - 1.0) * localScale;
            double value = (1.0 - mix) * b + mix * g + noise;
            if (value < 0.0) value = 0.0;
            if (value > 1.0) value = 1.0;
            injected.dec[(size_t)i * injected.dim + d] = value;
        }
    }
    evaluate(injected, funcNum);
    return injected;
}

static double mean_cv(const Population &pop)
{
    double s = 0.0;
    for (int i = 0; i < pop.n; ++i) s += constraint_violation(pop, i);
    return s / (double)pop.n;
}

static double score_cv(const Population &pop)
{
    double s = 0.0;
    for (int i = 0; i < pop.n; ++i)
    {
        const size_t base = (size_t)i * pop.ncon;
        for (int j = 0; j < pop.ncon; ++j) s += pop.con[base + j];
    }
    return s / (double)pop.n;
}

static double feasible_rate(const Population &pop)
{
    int cnt = 0;
    for (int i = 0; i < pop.n; ++i)
    {
        if (constraint_violation(pop, i) <= 0.0) ++cnt;
    }
    return (double)cnt / (double)pop.n;
}

static double best_cv(const Population &pop)
{
    double b = std::numeric_limits<double>::infinity();
    for (int i = 0; i < pop.n; ++i) b = std::min(b, constraint_violation(pop, i));
    return b;
}

static void print_state_trace(const char *source, int generation, int fe, double epsilon,
                              double nextRand, const Population &pop)
{
    double sumDec = 0.0;
    double sumObj = 0.0;
    double sumCon = 0.0;
    double sumPosCon = 0.0;
    double sumFit = 0.0;
    double rowWeightedDec = 0.0;
    double rowWeightedObj = 0.0;
    double rowWeightedCon = 0.0;
    double rowWeightedFit = 0.0;
    double minFit = std::numeric_limits<double>::infinity();
    double maxFit = -std::numeric_limits<double>::infinity();
    for (double x : pop.dec) sumDec += x;
    for (double x : pop.obj) sumObj += x;
    for (double x : pop.con)
    {
        sumCon += x;
        if (x > 0.0) sumPosCon += x;
    }
    for (double x : pop.fitness)
    {
        sumFit += x;
        if (x < minFit) minFit = x;
        if (x > maxFit) maxFit = x;
    }
    for (int i = 0; i < pop.n; ++i)
    {
        const double w = (double)(i + 1);
        for (int d = 0; d < pop.dim; ++d) rowWeightedDec += w * pop.dec[(size_t)i * pop.dim + d];
        for (int m = 0; m < pop.nobj; ++m) rowWeightedObj += w * pop.obj[(size_t)i * pop.nobj + m];
        for (int j = 0; j < pop.ncon; ++j) rowWeightedCon += w * pop.con[(size_t)i * pop.ncon + j];
        if (i < (int)pop.fitness.size()) rowWeightedFit += w * pop.fitness[(size_t)i];
    }
    if (pop.fitness.empty())
    {
        minFit = std::numeric_limits<double>::quiet_NaN();
        maxFit = std::numeric_limits<double>::quiet_NaN();
    }
    std::printf("STATE source=%s gen=%d fe=%d epsilon=%.17g next_rand=%.17g sum_dec=%.17g sum_obj=%.17g sum_con=%.17g sum_pos_con=%.17g sum_fit=%.17g row_weighted_dec=%.17g row_weighted_obj=%.17g row_weighted_con=%.17g row_weighted_fit=%.17g min_fit=%.17g max_fit=%.17g mean_cv=%.17g feasible_rate=%.17g\n",
                source, generation, fe, epsilon, nextRand, sumDec, sumObj, sumCon, sumPosCon,
                sumFit, rowWeightedDec, rowWeightedObj, rowWeightedCon, rowWeightedFit,
                minFit, maxFit, mean_cv(pop), feasible_rate(pop));
}

static void print_selection_probe(const char *source, int generation, double epsilon,
                                  const Population &pop, int targetN)
{
    const char *dumpRaw = std::getenv("RDE_TRACE_DUMP_GEN");
    const int dumpGeneration = dumpRaw && *dumpRaw ? std::atoi(dumpRaw) : -1;
    std::vector<int> feas, infeas;
    for (int i = 0; i < pop.n; ++i)
    {
        if (constraint_violation(pop, i) <= epsilon) feas.push_back(i);
        else infeas.push_back(i);
    }
    const char *branch = "mixed";
    std::vector<double> fit;
    int limit = targetN;
    if (feas.empty())
    {
        branch = "no_feasible";
        fit = cal_fitness_no_epsilon(pop, infeas);
        limit = targetN;
    }
    else if ((int)feas.size() <= targetN)
    {
        branch = "some_feasible";
        fit = feasible_augmented_fitness(pop, feas);
        limit = (int)feas.size();
    }
    else
    {
        branch = "too_many_feasible";
        fit = feasible_augmented_fitness(pop, feas);
        limit = targetN;
    }
    int lt1 = 0;
    double sumFit = 0.0;
    double weightedFit = 0.0;
    double minFit = std::numeric_limits<double>::infinity();
    double maxFit = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < (int)fit.size(); ++i)
    {
        if (fit[(size_t)i] < 1.0) ++lt1;
        sumFit += fit[(size_t)i];
        weightedFit += (double)(i + 1) * fit[(size_t)i];
        if (fit[(size_t)i] < minFit) minFit = fit[(size_t)i];
        if (fit[(size_t)i] > maxFit) maxFit = fit[(size_t)i];
    }
    if (fit.empty())
    {
        minFit = std::numeric_limits<double>::quiet_NaN();
        maxFit = std::numeric_limits<double>::quiet_NaN();
    }
    const int removeCount = lt1 > limit ? lt1 - limit : 0;
    std::printf("SELECT source=%s gen=%d epsilon=%.17g branch=%s feas=%d infeas=%d fit_n=%zu limit=%d lt1=%d remove=%d sum_fit=%.17g weighted_fit=%.17g min_fit=%.17g max_fit=%.17g\n",
                source, generation, epsilon, branch, (int)feas.size(), (int)infeas.size(),
                fit.size(), limit, lt1, removeCount, sumFit, weightedFit, minFit, maxFit);
    if (generation == dumpGeneration)
    {
        const std::vector<int> &dumpIdx = feas.empty() ? infeas : feas;
        for (int local = 0; local < (int)dumpIdx.size(); ++local)
        {
            const int src = dumpIdx[(size_t)local];
            std::printf("AUG source=%s gen=%d row=%d src=%d", source, generation, local + 1, src + 1);
            for (int m = 0; m < pop.nobj; ++m)
            {
                std::printf(" obj%d=%.17g", m + 1, pop.obj[(size_t)src * pop.nobj + m]);
            }
            std::printf(" cv=%.17g fit=%.17g\n", constraint_violation(pop, src),
                        local < (int)fit.size() ? fit[(size_t)local] : std::numeric_limits<double>::quiet_NaN());
            std::printf("DEC source=%s gen=%d row=%d src=%d", source, generation, local + 1, src + 1);
            for (int d = 0; d < pop.dim; ++d)
            {
                std::printf(" x%d=%.17g", d + 1, pop.dec[(size_t)src * pop.dim + d]);
            }
            std::printf("\n");
        }
    }
}

static void print_trace(const char *tag, const Options &opt, int fe, const Population &pop,
                        const ReferenceFront &ref, double elapsed)
{
    const double igd = opt.igd ? igd_value(pop, ref) : std::numeric_limits<double>::quiet_NaN();
    std::printf("%s problem=%d run=%d fe=%d igd=%.17g score_cv=%.17g mean_cv=%.17g feasible_rate=%.17g best_cv=%.17g elapsed=%.6f\n",
                tag, opt.problem, opt.run, fe, igd, score_cv(pop), mean_cv(pop), feasible_rate(pop),
                best_cv(pop), elapsed);
}

int main(int argc, char **argv)
{
    Options opt = parse_options(argc, argv);
    const int nobj = cmop_cec2026_nobj(opt.problem);
    const int ncon = cmop_cec2026_ng(opt.problem);
    if (opt.pfRoot.empty()) opt.pfRoot = default_pf_root();
    ReferenceFront ref = opt.igd ? load_reference_front(opt.problem, nobj, opt.pfRoot) : ReferenceFront();
    MatlabRng rng((unsigned long)opt.seed);
    const double t0 = now_seconds();
    TopRaceParams activeTopRace = read_toprace_params();
    ArchiveParams activeArchive = read_archive_params();
    activeTopRace.enabled = activeTopRace.enabled && opt.toprace;
    const int archiveCap = std::max(0, (int)std::llround(activeArchive.capFrac * (double)opt.popsize));
    std::printf("CONFIG problem=%d run=%d toprace=%d target_top_frac=%.17g sample_top_frac=%.17g progress_min=%.17g progress_max=%.17g feas_rate_min=%.17g use_f_scale=%d cr_source=%s bound_repair=%s bound_repair_mix=%.17g insert_mode=%s\n",
                opt.problem, opt.run, activeTopRace.enabled ? 1 : 0,
                activeTopRace.targetTopFrac, activeTopRace.sampleTopFrac,
                activeTopRace.progressMin, activeTopRace.progressMax, activeTopRace.feasRateMin,
                activeTopRace.useFScale ? 1 : 0, activeTopRace.crSource.c_str(),
                activeTopRace.boundRepair.c_str(), activeTopRace.boundRepairMix,
                activeTopRace.insertMode.c_str());
    std::printf("ARCHIVE_CONFIG enabled=%d output=%d reinject=%d cap=%d cap_frac=%.17g output_progress=%.17g reinject_frac=%.17g reinject_progress_min=%.17g reinject_interval=%.17g\n",
                activeArchive.enabled ? 1 : 0, activeArchive.output ? 1 : 0, activeArchive.reinject ? 1 : 0,
                archiveCap, activeArchive.capFrac, activeArchive.outputProgress, activeArchive.reinjectFrac,
                activeArchive.reinjectProgressMin, activeArchive.reinjectInterval);

    Population pop(opt.popsize, opt.dim, nobj, ncon);
    initialize(pop, opt.problem, rng);
    int fe = opt.popsize;

    std::vector<int> all(pop.n);
    std::iota(all.begin(), all.end(), 0);
    std::vector<double> cv0(pop.n);
    for (int i = 0; i < pop.n; ++i) cv0[i] = constraint_violation(pop, i);
    double epsilon0 = *std::max_element(cv0.begin(), cv0.end());
    if (epsilon0 == 0.0) epsilon0 = 1.0;
    pop.fitness = cal_fitness(pop, all, epsilon0);
    Population onlineArchive(0, opt.dim, nobj, ncon);
    if (activeArchive.enabled) onlineArchive = update_online_archive(onlineArchive, pop, archiveCap);
    int nextArchiveReinjectFE = (int)std::ceil(activeArchive.reinjectProgressMin * (double)opt.maxFE);
    int generation = 0;
    if (opt.stateTrace && (opt.stateTraceLimit <= 0 || generation <= opt.stateTraceLimit))
    {
        MatlabRng preview = rng;
        print_state_trace("cpp", generation, fe, epsilon0, preview.rand(), pop);
    }

    int nextReportFE = opt.save > 0 ? std::max(1, opt.maxFE / opt.save) : opt.maxFE + 1;
    int reportIndex = 1;
    while (fe < opt.maxFE)
    {
        double cp = (-std::log(epsilon0) - 6.0) / std::log(1.0 - 0.5);
        double progressLeft = std::max(0.0, 1.0 - (double)fe / (double)opt.maxFE);
        double epsilon = epsilon0 * std::pow(progressLeft, cp);
        int remaining = opt.maxFE - fe;
        OffspringResult offspring = make_offspring(pop, pop.fitness, opt, remaining, fe, opt.maxFE, epsilon, activeTopRace, rng);
        int offspringFE = fe + offspring.evaluations;
        if (activeArchive.enabled) onlineArchive = update_online_archive(onlineArchive, offspring.pop, archiveCap);
        const double archiveProgress = (double)offspringFE / (double)opt.maxFE;
        if (activeArchive.enabled && activeArchive.reinject && offspringFE >= nextArchiveReinjectFE)
        {
            const int remainingAfterOffspring = std::max(0, opt.maxFE - offspringFE);
            Population injected = archive_reinject(onlineArchive, pop, activeArchive, opt.problem, archiveProgress, remainingAfterOffspring, rng);
            if (injected.n > 0)
            {
                offspring.pop = merge_pop(offspring.pop, injected);
                offspringFE += injected.n;
                onlineArchive = update_online_archive(onlineArchive, injected, archiveCap);
            }
            nextArchiveReinjectFE += std::max(1, (int)std::ceil(activeArchive.reinjectInterval * (double)opt.maxFE));
        }
        if (opt.stateTrace && (opt.stateTraceLimit <= 0 || generation + 1 <= opt.stateTraceLimit))
        {
            if ((int)offspring.pop.fitness.size() != offspring.pop.n)
            {
                offspring.pop.fitness.assign(offspring.pop.n, 0.0);
            }
            MatlabRng preview = rng;
            print_state_trace("cpp_offspring", generation + 1, offspringFE, epsilon, preview.rand(), offspring.pop);
        }
        fe = offspringFE;
        Population merged = merge_pop(pop, offspring.pop);
        if (opt.stateTrace && (opt.stateTraceLimit <= 0 || generation + 1 <= opt.stateTraceLimit))
        {
            print_selection_probe("cpp", generation + 1, epsilon, merged, opt.popsize);
        }
        environmental_selection(merged, opt.popsize, epsilon);
        pop = std::move(merged);
        if (activeArchive.enabled) onlineArchive = update_online_archive(onlineArchive, pop, archiveCap);
        ++generation;
        if (opt.stateTrace && (opt.stateTraceLimit <= 0 || generation <= opt.stateTraceLimit))
        {
            MatlabRng preview = rng;
            print_state_trace("cpp", generation, fe, epsilon, preview.rand(), pop);
        }

        while (opt.save > 0 && reportIndex <= opt.save && fe >= reportIndex * nextReportFE)
        {
            const double reportProgress = (double)fe / (double)opt.maxFE;
            Population outputPop = archive_output_population(pop, onlineArchive, opt.popsize, reportProgress, activeArchive);
            print_trace("TRACE", opt, fe, outputPop, ref, now_seconds() - t0);
            ++reportIndex;
        }
    }

    Population finalOutput = archive_output_population(pop, onlineArchive, opt.popsize, 1.0, activeArchive);
    const double finalIgd = opt.igd ? igd_value(finalOutput, ref) : std::numeric_limits<double>::quiet_NaN();
    std::printf("RESULT problem=%d run=%d seed=%d N=%d D=%d maxFE=%d finalFE=%d igd=%.17g score_cv=%.17g mean_cv=%.17g feasible_rate=%.17g best_cv=%.17g archive_size=%d elapsed=%.6f\n",
                opt.problem, opt.run, opt.seed, opt.popsize, opt.dim, opt.maxFE, fe,
                finalIgd, score_cv(finalOutput), mean_cv(finalOutput), feasible_rate(finalOutput), best_cv(finalOutput),
                onlineArchive.n, now_seconds() - t0);
    return 0;
}
