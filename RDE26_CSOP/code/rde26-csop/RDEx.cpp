#include <iostream>
#include <time.h>
#include <fstream>
#include <random>
#include <chrono>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include "cec17_test_COP.cpp"

const int ResTsize1 = 28; // number of functions //28 for CEC 2024(2017)
const int ResTsize2 = 2001; // number of records per function //2000+1 = 2001 for CEC 2024

using namespace std;
/*typedef std::chrono::high_resolution_clock myclock;
myclock::time_point beginning = myclock::now();
myclock::duration d1 = myclock::now() - beginning;
#ifdef __linux__
    unsigned globalseed = d1.count();
#elif _WIN32
    unsigned globalseed = unsigned(time(NULL));
#else

#endif*/
unsigned globalseed = 20260611u;
unsigned seed1 = globalseed+0;
unsigned seed2 = globalseed+100;
unsigned seed3 = globalseed+200;
unsigned seed4 = globalseed+300;
unsigned seed5 = globalseed+400;
unsigned seed6 = globalseed+500;
unsigned seed7 = globalseed+600;
std::mt19937 generator_uni_i(seed1);
std::mt19937 generator_uni_r(seed2);
std::mt19937 generator_norm(seed3);
std::mt19937 generator_cachy(seed4);
std::mt19937 generator_uni_i_3(seed6);
std::uniform_int_distribution<int> uni_int(0,32768);
std::uniform_real_distribution<double> uni_real(0.0,1.0);
std::normal_distribution<double> norm_dist(0.0,1.0);
std::cauchy_distribution<double> cachy_dist(0.0,1.0);

static unsigned SeedBase_global = 20260611u;

static void reset_rng_for_case(int func_num, int run)
{
    unsigned s = SeedBase_global + unsigned(func_num) * 1009u + unsigned(run) * 9176u;
    generator_uni_i.seed(s + 0u);
    generator_uni_r.seed(s + 100u);
    generator_norm.seed(s + 200u);
    generator_cachy.seed(s + 300u);
    generator_uni_i_3.seed(s + 600u);
    srand(s + 700u);
}

const char* algName = "CORDEx_2026_Submission_R39";

// Trust-region parameter block generated from JSON.
const int UDEFusionMode = 3;
const int UDEFusionArchiveRecover = 0;
const double TR_TARGET_TOP_FRAC = 0.080000000000000002;
const double TR_SAMPLE_TOP_FRAC = 0.27000000000000002;
const double TR_PROGRESS_MIN = 0;
const double TR_PROGRESS_MAX = 0.76000000000000001;
const double TR_FEAS_RATE_MIN = 0.48000000000000004;
const double TR_ALT_CR_MIN = 0.26000000000000001;
const double TR_ALT_CR_MAX = 0.85000000000000009;
const double TR_PBEST_PULL_SCALE = 0.62;
const double TR_DIFF_SCALE = 1.1600000000000001;
const int TR_PBEST_RETRY = 7;
const double TR_EQ_EPSILON = 0.0001;
const double TR_EPSILON_CUTOFF_PROGRESS = 0.75;
const double TR_EPSILON_INDEX_PARAM = 0.80000000000000004;
const int TR_MEMORY_SIZE = 11;
const double TR_MEMORY_CR_INIT = 0.94000000000000006;
const double TR_MEMORY_F_INIT = 0.63;
const double TR_MEMORY_LEARNING_RATE = 0.34999999999999998;
const double TR_MEMORY_F_SELECT_PROB = 0.080000000000000002;
const double TR_SUCCESS_RATE_EXPONENT = 0.40000000000000002;
const double TR_F_SIGMA = 0.050000000000000003;
const double TR_F_CAUCHY_SIGMA = 0.13999999999999999;
const double TR_CR_SIGMA = 0.090000000000000011;
const double TR_PSIZE_BASE = 0.5;
const double TR_FRONT_RANK_PRESSURE = 13.98;
const double TR_PRAND_WEIGHT_SCALE = 4.8700000000000001;
const double TR_PSIZE2_BASE = 0.40000000000000002;
const double TR_PSIZE2_DECAY = 0.35000000000000003;
const double TR_EB_RATE_INIT = 0.91000000000000003;
const double TR_EB_DISABLE_BEFORE_PROGRESS = 0.70999999999999996;
const double TR_EB_F_FALLBACK_MEAN = 0.43000000000000005;
const double TR_EB_CR_FALLBACK_MEAN = 0.63;
const double TR_CR_EARLY_PROGRESS_1 = 0.28999999999999998;
const double TR_CR_EARLY_MIN_1 = 0.23000000000000001;
const double TR_CR_EARLY_PROGRESS_2 = 0.63;
const double TR_CR_EARLY_MIN_2 = 0.68999999999999995;
const double TR_PERTURB_PROB = 0.55000000000000004;
const double TR_PERTURB_SCALE = 0.29999999999999999;
const double TR_BOUND_REPAIR_MIX = 0.21000000000000002;
const double TR_POP_SIZE_FACTOR = 9.9800000000000004;
const double UDEFusionTopFrac = TR_TARGET_TOP_FRAC;
const double UDEFusionPBestSampleFrac = TR_SAMPLE_TOP_FRAC;
const double UDEFusionProgressMin = TR_PROGRESS_MIN;
const double UDEFusionProgressMax = TR_PROGRESS_MAX;
const int UDEFusionAltFormula = 2;

inline double TRMemoryUpdate(double observed, double stored)
{
    if(TR_MEMORY_LEARNING_RATE == 0.5)
        return 0.5*(observed + stored);
    return TR_MEMORY_LEARNING_RATE * observed + (1.0 - TR_MEMORY_LEARNING_RATE) * stored;
}

inline double TRBoundRepair(double current, double bound)
{
    if(TR_BOUND_REPAIR_MIX == 0.5)
        return (current + bound)*0.5;
    return (1.0 - TR_BOUND_REPAIR_MIX) * current + TR_BOUND_REPAIR_MIX * bound;
}
double EB_hybrid_rate_init = TR_EB_RATE_INIT;

int IntRandom(int target) {if(target == 0) return 0; return uni_int(generator_uni_i)%target;}
double Random(double minimal, double maximal){return uni_real(generator_uni_r)*(maximal-minimal)+minimal;}
double NormRand(double mu, double sigma){return norm_dist(generator_norm)*sigma + mu;}
double CachyRand(double mu, double sigma){return cachy_dist(generator_cachy)*sigma+mu;}

void cec17_test_COP(double *x, double *f, double *g,double *h, int nx, int mx,int func_num);

const int ng_B[28]={1,1,1,2,2,0,0,0,1,0,1,2,3,1,1,1,1,2,2,2,2,3,1,1,1,1,2,2};
const int nh_B[28]={0,0,1,0,0,6,2,2,1,2,1,0,0,1,1,1,1,1,0,0,0,0,1,1,1,1,1,0};
                        //1   2   3  4  5  6  7   8  9  10  11  12  13  14  15  16  17  18 19  20  21  22  23  24  25  26  27 28
const int border[28] = {100,100,100,10,10,20,50,100,10,100,100,100,100,100,100,100,100,100,50,100,100,100,100,100,100,100,100,50};

double *OShift=NULL,*M=NULL,*M1=NULL,*M2=NULL,*y=NULL,*z=NULL,*z1=NULL,*z2=NULL;
int ini_flag=0,n_flag,func_flag,f5_flag;
int stepsFEval[ResTsize2-1];
double ResultsArray[ResTsize2];
double ResultsArrayG[ResTsize2][3];
double ResultsArrayH[ResTsize2][6];
int LastFEcount;
int NFEval = 0;
int MaxFEval = 0;
int GNVars;
double tempF[1];
double tempG[3];
double tempH[6];
double xopt[100];
double fopt[1];
char buffer[500];
double globalbest;
double globalbestpenalty;
bool globalbestinit;
bool TimeComplexity = false;
double epsilon0001 = TR_EQ_EPSILON;
double PRS_mF[16];
double PRS_sF[16];
double PRS_kF[16];
double Cvalglobal = 4;

int UDESelectBranch(double* branch_prob)
{
    double total = branch_prob[0] + branch_prob[1] + branch_prob[2];
    if(total <= 0.0)
        return 0;
    double r = Random(0, 1) * total;
    if(r < branch_prob[1])
        return 1;
    if(r < branch_prob[1] + branch_prob[2])
        return 2;
    return 0;
}

bool UDEBetterByEpsilon(double fit_new, double pen_new, double fit_old, double pen_old, double epsilon)
{
    double pn = pen_new;
    double po = pen_old;
    bool good_new = false;
    if(pn <= epsilon)
    {
        pn = 0;
        good_new = true;
    }
    if(po <= epsilon)
        po = 0;
    if((good_new && fit_new <= fit_old) || ((pn == po) && (fit_new <= fit_old)))
        return true;
    if(pn < po)
        return true;
    return false;
}

void qSort2int(double* Mass, int* Mass2, int low, int high)
{
    int i=low;
    int j=high;
    double x=Mass[(low+high)>>1];
    do
    {
        while(Mass[i]<x)    ++i;
        while(Mass[j]>x)    --j;
        if(i<=j)
        {
            double temp=Mass[i];
            Mass[i]=Mass[j];
            Mass[j]=temp;
            int temp2=Mass2[i];
            Mass2[i]=Mass2[j];
            Mass2[j]=temp2;
            i++;    j--;
        }
    } while(i<=j);
    if(low<j)   qSort2int(Mass,Mass2,low,j);
    if(i<high)  qSort2int(Mass,Mass2,i,high);
}
void getOptimum(const int func_num)
{
    FILE *fpt=NULL;
	char FileName[30];
    sprintf(FileName, "inputData/shift_data_%d.txt", func_num);
    fpt = fopen(FileName,"r");
    if (fpt==NULL)
        printf("\n Error: Cannot open input file for reading \n");
    for(int k=0;k!=GNVars;k++)
    {
        fscanf(fpt,"%lf",&xopt[k]);
    }
    fclose(fpt);
    cec17_test_COP(xopt, fopt, tempG, tempH, GNVars, 1, func_num);
}
double cec_24_constr(double* HostVector, const int func_num)
{
    tempG[0] = 0;tempG[1] = 0;tempG[2] = 0;
    tempH[0] = 0;tempH[1] = 0;tempH[2] = 0;tempH[3] = 0;tempH[4] = 0;tempH[5] = 0;
    cec17_test_COP(HostVector, tempF, tempG, tempH, GNVars, 1, func_num);
    NFEval++;
    return tempF[0];
}
double cec_24_totalpenalty(const int func_num, double* PopulG, double* PopulH, double* PopulAllConstr)
{
    for(int i=0;i!=3;i++)
        PopulG[i] = tempG[i];
    for(int i=0;i!=6;i++)
    {
        PopulH[i] = tempH[i];
        if(PopulH[i] < epsilon0001 && PopulH[i] > -epsilon0001)
            PopulH[i] = 0;
    }
    for(int i=0;i!=ng_B[func_num-1]+nh_B[func_num-1];i++)
    {
        if(i < ng_B[func_num-1])
            PopulAllConstr[i] = tempG[i];
        else if(PopulH[i-ng_B[func_num-1]] >= epsilon0001)
            PopulAllConstr[i] = tempH[i-ng_B[func_num-1]];
        else if(PopulH[i-ng_B[func_num-1]] <= -epsilon0001)
            PopulAllConstr[i] = -tempH[i-ng_B[func_num-1]];
        else
            PopulAllConstr[i] = 0;
    }
    double total = 0;
    for(int i=0;i!=ng_B[func_num-1];i++)
        if(tempG[i] > 0)
            total += tempG[i];
    for(int i=0;i!=nh_B[func_num-1];i++)
        if(tempH[i] >= epsilon0001)
            total += tempH[i];
        else if(tempH[i] <= -epsilon0001)
            total += -tempH[i];
    return total/(double(ng_B[func_num-1])+double(nh_B[func_num-1]));
}
void SaveBestValues(int func_num, double* BestG, double* BestH)
{
    double temp = globalbest;// - fopt[0];
    if(temp <= 1E-8 && globalbestpenalty <= 1E-8 && ResultsArray[ResTsize2-1] == MaxFEval)
    {
        ResultsArray[ResTsize2-1] = NFEval;
    }
    for(int stepFEcount=LastFEcount;stepFEcount<ResTsize2-1;stepFEcount++)
    {
        if(NFEval == stepsFEval[stepFEcount])
        {
            ResultsArray[stepFEcount] = temp;
            for(int i=0;i!=3;i++)
                ResultsArrayG[stepFEcount][i] = BestG[i];
            for(int i=0;i!=6;i++)
                ResultsArrayH[stepFEcount][i] = BestH[i];
            LastFEcount = stepFEcount;
        }
    }
}
class Optimizer
{
public:
    int MemorySize;
    int MemoryIter;
    int SuccessFilled;
    int MemoryCurrentIndex;
    int MemoryCurrentIndex2;
    int NVars;			    // размерность пространства
    int NIndsCurrent;
    int NIndsFront;
    int NIndsFrontMax;
    int newNIndsFront;
    int PopulSize;
    int func_num;
    int TheChosenOne;
    int Generation;
    int PFIndex;
    int NConstr;

    double bestfit;
    double SuccessRate;
    double F;       /*параметры*/
    double Cr;
    double Right;		    // верхняя граница
    double Left;		    // нижняя граница

    double** Popul;	        // массив для частиц
    double** PopulG;
    double** PopulH;
    double** PopulAllConstr;
    double** PopulFront;
    double** PopulFrontG;
    double** PopulFrontH;
    double** PopulAllConstrFront;
    double** PopulTemp;
    double** PopulTempG;
    double** PopulTempH;
    double** PopulAllConstrTemp;
    double* FitArr;		// значения функции пригодности
    double* FitArrTemp;		// значения функции пригодности
    double* FitArrCopy;
    double* FitArrFront;
    double* Trial;
    double* tempSuccessCr;
    double* tempSuccessF;
    double* MemoryCr;
    double* MemoryF;
    double* FitDelta;
    double* Weights;
    double* PenaltyArr;
    double* PenaltyArrFront;
    double* PenaltyTemp;
    double* BestInd;
    double BestG[3];
    double BestH[6];
    double* EpsLevels;
    double* massvector;
    double* tempvector;
    double* epsvector;

    int* Indices;
    int* IndicesFront;
    int* IndicesConstr;
    int* IndicesConstrFront;

    
    double* ord_best_arch;
    double* ord_medium_arch;
    double* ord_worst_arch;
    double* ord_best_popul;
    double* ord_medium_popul;
    double* ord_worst_popul;
    double* FitMass;
    double EB_hybrid_rate;
    double* FitTemp;

    void Initialize(int _newNInds, int _newNVars, int _newfunc_num);
    void Clean();
    void MainCycle();
    void UpdateMemoryCr();
    double MeanWL(double* Vector, double* TempWeights);
    void RemoveWorst(int NInds, int NewNInds, double epsilon);

    void EB_order(int prand, int Rand1, int Rand2);
    double* EB_hybrid_flag;
    void UpdateEB_hybrid_param(double* EB_hybrid_flag, double* FitArrFront, vector<double> FitTemp);
};

void Optimizer::Initialize(int _newNInds, int _newNVars, int _newfunc_num)
{
    NVars = _newNVars;
    NIndsCurrent = _newNInds;
    NIndsFront = _newNInds;
    NIndsFrontMax = _newNInds;
    PopulSize = _newNInds*2;
    Generation = 0;
    TheChosenOne = 0;
    MemorySize = TR_MEMORY_SIZE;
    MemoryIter = 0;
    SuccessFilled = 0;
    SuccessRate = 0.5;
    func_num = _newfunc_num;
    NConstr = ng_B[func_num-1] + nh_B[func_num-1];
    Left = -border[func_num-1];
    Right = border[func_num-1];
    for(int steps_k=0;steps_k!=ResTsize2-1;steps_k++)
        stepsFEval[steps_k] = 20000.0/double(ResTsize2-1)*GNVars*(steps_k+1);
    EpsLevels = new double[NConstr];

    Popul = new double*[PopulSize];
    for(int i=0;i!=PopulSize;i++)
        Popul[i] = new double[NVars];
    PopulG = new double*[PopulSize];
    for(int i=0;i!=PopulSize;i++)
        PopulG[i] = new double[3];
    PopulH = new double*[PopulSize];
    for(int i=0;i!=PopulSize;i++)
        PopulH[i] = new double[6];
    PopulAllConstr = new double*[PopulSize];
    for(int i=0;i!=PopulSize;i++)
        PopulAllConstr[i] = new double[NConstr];

    PopulFront = new double*[NIndsFront];
    for(int i=0;i!=NIndsFront;i++)
        PopulFront[i] = new double[NVars];
    PopulFrontG = new double*[NIndsFront];
    for(int i=0;i!=NIndsFront;i++)
        PopulFrontG[i] = new double[3];
    PopulFrontH = new double*[NIndsFront];
    for(int i=0;i!=NIndsFront;i++)
        PopulFrontH[i] = new double[6];
    PopulAllConstrFront = new double*[PopulSize];
    for(int i=0;i!=PopulSize;i++)
        PopulAllConstrFront[i] = new double[NConstr];

    PopulTemp = new double*[PopulSize];
    for(int i=0;i!=PopulSize;i++)
        PopulTemp[i] = new double[NVars];
    PopulTempG = new double*[PopulSize];
    for(int i=0;i!=PopulSize;i++)
        PopulTempG[i] = new double[3];
    PopulTempH = new double*[PopulSize];
    for(int i=0;i!=PopulSize;i++)
        PopulTempH[i] = new double[6];
    PopulAllConstrTemp = new double*[PopulSize];
    for(int i=0;i!=PopulSize;i++)
        PopulAllConstrTemp[i] = new double[NConstr];

    FitArr = new double[PopulSize];
    FitArrTemp = new double[PopulSize];
    FitArrCopy = new double[PopulSize];
    FitArrFront = new double[NIndsFront];

    PenaltyArr = new double[PopulSize];
    PenaltyTemp = new double[PopulSize];
    PenaltyArrFront = new double[PopulSize];

    Weights = new double[PopulSize];
    tempSuccessCr = new double[PopulSize];
    tempSuccessF = new double[PopulSize];
    FitDelta = new double[PopulSize];
    MemoryCr = new double[MemorySize];
    MemoryF = new double[MemorySize];
    Trial = new double[NVars];
    BestInd = new double[NVars];
    massvector=new double[NConstr+1];
    tempvector=new double[NConstr+1];
    epsvector =new double[NConstr+1];

    Indices = new int[PopulSize];
    IndicesFront = new int[PopulSize];
    IndicesConstr = new int[PopulSize];
    IndicesConstrFront = new int[PopulSize];

	for (int i = 0; i<PopulSize; i++)
		for (int j = 0; j<NVars; j++)
			Popul[i][j] = Random(Left,Right);
    for(int i=0;i!=PopulSize;i++)
        tempSuccessCr[i] = 0;
    for(int i=0;i!=MemorySize;i++)
        MemoryCr[i] = TR_MEMORY_CR_INIT;
    for(int i=0;i!=MemorySize;i++)
        MemoryF[i] = TR_MEMORY_F_INIT;

    EB_hybrid_flag = new double[NIndsFrontMax];
    FitMass = new double[NIndsFrontMax];
    FitTemp = new double[NIndsFrontMax];
    EB_hybrid_rate = EB_hybrid_rate_init;
}
void Optimizer::Clean()
{
    delete Trial;
    delete BestInd;
    for(int i=0;i!=PopulSize;i++)
    {
        delete Popul[i];
        delete PopulG[i];
        delete PopulH[i];
        delete PopulTemp[i];
        delete PopulTempG[i];
        delete PopulTempH[i];
        delete PopulAllConstr[i];
        delete PopulAllConstrTemp[i];
    }
    for(int i=0;i!=NIndsFrontMax;i++)
    {
        delete PopulFront[i];
        delete PopulFrontG[i];
        delete PopulFrontH[i];
    }
    delete Popul;
    delete PopulG;
    delete PopulH;
    delete PopulTemp;
    delete PopulTempG;
    delete PopulTempH;
    delete PopulAllConstr;
    delete PopulAllConstrTemp;
    delete PopulAllConstrFront;
    delete PopulFront;
    delete PopulFrontG;
    delete PopulFrontH;
    delete FitArr;
    delete FitArrTemp;
    delete FitArrCopy;
    delete FitArrFront;
    delete PenaltyArr;
    delete PenaltyArrFront;
    delete PenaltyTemp;
    delete Indices;
    delete IndicesFront;
    delete IndicesConstr;
    delete IndicesConstrFront;
    delete tempSuccessCr;
    delete tempSuccessF;
    delete FitDelta;
    delete MemoryCr;
    delete MemoryF;
    delete Weights;
    delete EpsLevels;
    delete massvector;
    delete tempvector;
    delete epsvector;

    delete FitMass;
    delete EB_hybrid_flag;
    delete FitTemp;
}
void Optimizer::UpdateMemoryCr()
{
    if(SuccessFilled != 0)
    {
        MemoryCr[MemoryIter] = TRMemoryUpdate(MeanWL(tempSuccessCr,FitDelta), MemoryCr[MemoryIter]);
        MemoryF[MemoryIter] = TRMemoryUpdate(MeanWL(tempSuccessF,FitDelta), MemoryF[MemoryIter]);
        MemoryIter = (MemoryIter+1)%MemorySize;
    }
}
double Optimizer::MeanWL(double* Vector, double* TempWeights)
{
    double SumWeight = 0;
    double SumSquare = 0;
    double Sum = 0;
    for(int i=0;i!=SuccessFilled;i++)
        SumWeight += TempWeights[i];
    for(int i=0;i!=SuccessFilled;i++)
        Weights[i] = TempWeights[i]/SumWeight;
    for(int i=0;i!=SuccessFilled;i++)
        SumSquare += Weights[i]*Vector[i]*Vector[i];
    for(int i=0;i!=SuccessFilled;i++)
        Sum += Weights[i]*Vector[i];
    if(fabs(Sum) > 1e-8)
        return SumSquare/Sum;
    else
        return 1.0;
}

void Optimizer::RemoveWorst(int _NIndsFront, int _newNIndsFront, double epsilon)
{
    int PointsToRemove = _NIndsFront - _newNIndsFront;
    for(int L=0;L!=PointsToRemove;L++)
    {
        double maxfitFront = FitArrFront[0];
        for(int i=0;i!=NIndsFront-L;i++)
            maxfitFront = max(maxfitFront,FitArrFront[i]);
        double maxfit = FitArrFront[0];
        int WorstNum = 0;
        for(int i=0;i!=NIndsFront-L;i++)
        {
            if(PenaltyArrFront[i] > epsilon)
            {
                if(maxfitFront + 1.0 + PenaltyArrFront[i] > maxfit)
                {
                    maxfit = maxfitFront + 1.0 + PenaltyArrFront[i];
                    WorstNum = i;
                }
            }
            else
            {
                if(FitArrFront[i] > maxfit)
                {
                    maxfit = FitArrFront[i];
                    WorstNum = i;
                }
            }
        }
        for(int i=WorstNum;i!=_NIndsFront-1;i++)
        {
            for(int j=0;j!=NVars;j++)
                PopulFront[i][j] = PopulFront[i+1][j];
            for(int j=0;j!=3;j++)
                PopulFrontG[i][j] = PopulFrontG[i+1][j];
            for(int j=0;j!=6;j++)
                PopulFrontH[i][j] = PopulFrontH[i+1][j];
            for(int j=0;j!=NConstr;j++)
                PopulAllConstrFront[i][j] = PopulAllConstrFront[i+1][j];
            FitArrFront[i] = FitArrFront[i+1];
            FitMass[i] = FitMass[i+1];
            PenaltyArrFront[i] = PenaltyArrFront[i+1];
        }
    }
}
void Optimizer::EB_order(int prand, int Rand1, int Rand2) {//double* pos1, double* pos2, double* pos3, double* pos4, double* pos1_Fit, double* pos2_Fit, double* pos3_Fit, double* pos4_Fit
    double* pos1 = Popul[prand];
    double pos1Fit = FitArr[prand];
    double* pos3 = PopulFront[Rand1];
    double pos3Fit = FitArrFront[Rand1];
    double* pos4_arch = Popul[Rand2];
    double pos4Fit_arch = FitArr[Rand2];
    if (pos1Fit <= pos3Fit && pos1Fit <= pos4Fit_arch){
        ord_best_arch = pos1;
        if (pos3Fit <= pos4Fit_arch){
            ord_medium_arch = pos3;
            ord_worst_arch = pos4_arch;
        }
        else{
            ord_medium_arch = pos4_arch;
            ord_worst_arch = pos3;
        }
    }
    else if (pos3Fit <= pos1Fit && pos3Fit <= pos4Fit_arch){
        ord_best_arch = pos3;
        if (pos1Fit <= pos4Fit_arch){
            ord_medium_arch = pos1;
            ord_worst_arch = pos4_arch;
        }
        else{
            ord_medium_arch = pos4_arch;
            ord_worst_arch = pos1;
        }
    }
    else if (pos4Fit_arch <= pos1Fit && pos4Fit_arch <= pos3Fit){
        ord_best_arch = pos4_arch;
        if (pos1Fit <= pos3Fit){
            ord_medium_arch = pos1;
            ord_worst_arch = pos3;
        }
        else{
            ord_medium_arch = pos3;
            ord_worst_arch = pos1;
        }
    }
    double* pos4_popul = Popul[Rand2];
    double pos4Fit_popul = FitArr[Rand2];
    if (pos1Fit <= pos3Fit && pos1Fit <= pos4Fit_popul){
        ord_best_popul = pos1;
        if (pos3Fit <= pos4Fit_popul){
            ord_medium_popul = pos3;
            ord_worst_popul = pos4_popul;
        }
        else{
            ord_medium_popul = pos4_popul;
            ord_worst_popul = pos3;
        }
    }
    else if (pos3Fit <= pos1Fit && pos3Fit <= pos4Fit_popul){
        ord_best_popul = pos3;
        if (pos1Fit <= pos4Fit_popul){
            ord_medium_popul = pos1;
            ord_worst_popul = pos4_popul;
        }
        else{
            ord_medium_popul = pos4_popul;
            ord_worst_popul = pos1;
        }
    }
    else if (pos4Fit_popul <= pos1Fit && pos4Fit_popul <= pos3Fit){
        ord_best_popul = pos4_popul;
        if (pos1Fit <= pos3Fit){
            ord_medium_popul = pos1;
            ord_worst_popul = pos3;
        }
        else{
            ord_medium_popul = pos3;
            ord_worst_popul = pos1;
        }
    }
}
void Optimizer::UpdateEB_hybrid_param(double* EB_hybrid_flag, double* FitArrFront, vector<double> FitTemp) {
    double SumEB_DeltaFit = 0;
    double SumOrigin_DeltaFit = 0;
    for (int ChosenOne = 0; ChosenOne != NIndsFront; ChosenOne++) {
        if (EB_hybrid_flag[ChosenOne] == 1) {
            if (FitTemp[ChosenOne] <= FitArrFront[ChosenOne]) {
                SumEB_DeltaFit += FitArrFront[ChosenOne] - FitTemp[ChosenOne];
            }
        }
        else{
            if (FitTemp[ChosenOne] <= FitArrFront[ChosenOne]) {
                SumOrigin_DeltaFit += FitArrFront[ChosenOne] - FitTemp[ChosenOne];
            }
        }
    }

    if (SumEB_DeltaFit != 0 && SumOrigin_DeltaFit != 0){
        EB_hybrid_rate = SumEB_DeltaFit / (SumEB_DeltaFit + SumOrigin_DeltaFit);
        double EB_limit_min = 0;
        double EB_limit_max = 1;
        if (EB_hybrid_rate > EB_limit_max){
            EB_hybrid_rate = EB_limit_max;
        }
        else if (EB_hybrid_rate < EB_limit_min){
            EB_hybrid_rate = EB_limit_min;
        }
    }
    else{
        EB_hybrid_rate = EB_hybrid_rate_init;
    }
}
void Optimizer::MainCycle()
{
    double epsilon = TR_EQ_EPSILON;
    double epsilonf = TR_EQ_EPSILON;
    double ECutoffParam = TR_EPSILON_CUTOFF_PROGRESS;
    vector<double> FitTemp2;
    vector<double> FitTemp_prand;
    double UDEBranchProb[3] = {0.60, 0.25, 0.15};
    double UDEBranchMF[3] = {0.30, 0.40, 0.50};
    double UDEBranchMCR[3] = {1.00, 0.90, 0.70};
    int UDEStag[600];
    for(int i=0;i!=600;i++) UDEStag[i] = 0;
    double UDEArchiveX[32][100];
    double UDEArchiveG[32][3];
    double UDEArchiveH[32][6];
    double UDEArchiveAll[32][9];
    double UDEArchiveFit[32];
    double UDEArchivePenalty[32];
    int UDEArchiveSize = 0;
    int UDEArchiveCursor = 0;
    for(int IndIter=0;IndIter<NIndsFront;IndIter++)
    {
        FitArr[IndIter] = cec_24_constr(Popul[IndIter],func_num);
        PenaltyArr[IndIter] = cec_24_totalpenalty(func_num,PopulG[IndIter],PopulH[IndIter],PopulAllConstr[IndIter]);
        if(!globalbestinit ||
           ((FitArr[IndIter] <= globalbest) && PenaltyArr[IndIter] <= 0) ||
           ((PenaltyArr[IndIter] <= globalbestpenalty) && PenaltyArr[IndIter] > 0))
        {
            globalbest = FitArr[IndIter];
            globalbestpenalty = PenaltyArr[IndIter];
            globalbestinit = true;
            bestfit = FitArr[IndIter];
            for(int j=0;j!=NVars;j++)
                BestInd[j] = Popul[IndIter][j];
            for(int j=0;j!=3;j++)
                BestG[j] = PopulG[IndIter][j];
            for(int j=0;j!=6;j++)
                BestH[j] = PopulH[IndIter][j];
        }
        SaveBestValues(func_num,BestG,BestH);
    }
    for(int i=0;i!=NIndsFront;i++)
    {
        for(int j=0;j!=NVars;j++)
            PopulFront[i][j] = Popul[i][j];
        for(int j=0;j!=3;j++)
            PopulFrontG[i][j] = PopulG[i][j];
        for(int j=0;j!=6;j++)
            PopulFrontH[i][j] = PopulH[i][j];
        for(int j=0;j!=NConstr;j++)
            PopulAllConstrFront[i][j] = PopulAllConstr[i][j];
        FitArrFront[i] = FitArr[i];
        FitMass[i] = FitArrFront[i];
        PenaltyArrFront[i] = PenaltyArr[i];
    }
    PFIndex = 0;
    double EIndexParam = TR_EPSILON_INDEX_PARAM;
    while(NFEval < MaxFEval)
    {
        double meanF = max(0.0,pow(SuccessRate,TR_SUCCESS_RATE_EXPONENT));
        double sigmaF = TR_F_SIGMA;

        int epsilonindex = (NIndsFront*EIndexParam*
                        (1.0-(double)NFEval/(double)MaxFEval)*
                        (1.0-(double)NFEval/(double)MaxFEval));
        int epsilonfindex = (EIndexParam*NIndsFront*
                        (1.0-(double)NFEval/(double)MaxFEval)*
                        (1.0-(double)NFEval/(double)MaxFEval));

        int psizeval = max(2,int(TR_PSIZE_BASE*NIndsFront));


        double minfit = FitArrFront[0];
        double maxfit = FitArrFront[0];
        for(int i=0;i!=NIndsFront;i++)
        {
            FitArrCopy[i] = FitArrFront[i];
            Indices[i] = i;
            if(FitArrFront[i] >= maxfit)
                maxfit = FitArrFront[i];
            if(FitArrFront[i] <= minfit)
                minfit = FitArrFront[i];
        }
        if(minfit != maxfit)
            qSort2int(FitArrCopy,Indices,0,NIndsFront-1);
        epsilonf = FitArrCopy[epsilonfindex];
        if(NFEval > ECutoffParam*MaxFEval)
            epsilonf = minfit;

        for(int C=0;C!=NConstr;C++)
        {
            for(int i=0;i!=NIndsFront;i++)
            {
                FitArrCopy[i] = PopulAllConstrFront[i][C];
                if(FitArrCopy[i] < epsilon0001)
                    FitArrCopy[i] = 0;
                IndicesConstr[i] = i;
            }
            double penaltymax = FitArrCopy[0];
            double penaltymin = FitArrCopy[0];
            for(int i=0;i!=NIndsFront;i++)
            {
                if(FitArrCopy[i] >= penaltymax)
                    penaltymax = FitArrCopy[i];
                if(FitArrCopy[i] <= penaltymin)
                    penaltymin = FitArrCopy[i];
            }
            if(penaltymin != penaltymax)
                qSort2int(FitArrCopy,IndicesConstr,0,NIndsFront-1);
            if(EpsLevels[C] > FitArrCopy[epsilonindex] || Generation == 0)
                EpsLevels[C] = FitArrCopy[epsilonindex];
            if(NFEval > ECutoffParam*MaxFEval)
                EpsLevels[C] = 0.0;
        }

        int prand = 0;
        int Rand1 = 0;
        int Rand2 = 0;

        double FeasRate = 0;
        for(int i=0;i!=NIndsFront;i++)
        {
            if(PenaltyArrFront[i] < epsilon0001)
                FeasRate += 1;
        }
        FeasRate = FeasRate / double(NIndsFront);

        minfit = PenaltyArrFront[0];
        maxfit = PenaltyArrFront[0];
        for(int i=0;i!=NIndsFront;i++)
        {
            FitArrCopy[i] = PenaltyArrFront[i];
            IndicesConstr[i] = i;
            maxfit = max(maxfit,PenaltyArrFront[i]);
            minfit = min(minfit,PenaltyArrFront[i]);
        }
        if(minfit != maxfit)
            qSort2int(FitArrCopy,IndicesConstr,0,NIndsFront-1);
        epsilon = FitArrCopy[epsilonindex];
        if((double)NFEval/(double)MaxFEval > ECutoffParam)
            epsilon = 0.0;

        double maxfitPop = FitArr[0];
        for(int i=0;i!=NIndsFront;i++)
            maxfitPop = max(maxfitPop,FitArr[i]);
        //get indices for front
        for(int i=0;i!=NIndsFront;i++)
        {
            if(PenaltyArr[i] > epsilon)
                FitArrCopy[i] = maxfitPop + 1.0 + PenaltyArr[i];
            else
                FitArrCopy[i] = FitArr[i];
            if(i == 0)
            {
                minfit = FitArrCopy[0];
                maxfit = FitArrCopy[0];
            }
            Indices[i] = i;
            maxfit = max(maxfit,FitArr[i]);
            minfit = min(minfit,FitArr[i]);
        }
        if(minfit != maxfit)
            qSort2int(FitArrCopy,Indices,0,NIndsFront-1);

        //ranks for selective pressure
        FitTemp2.resize(NIndsFront);
        for(int i=0;i!=NIndsFront;i++)
            FitTemp2[i] = exp(-double(i)/double(NIndsFront)*TR_FRONT_RANK_PRESSURE);
        std::discrete_distribution<int> ComponentSelectorFront (FitTemp2.begin(),FitTemp2.end());

        FitTemp_prand.resize(NIndsFront);
        for (int i = 0; i != NIndsFront; i++)
            FitTemp_prand[i] = TR_PRAND_WEIGHT_SCALE * (NIndsFront - i);
        int psizeval2 = NIndsFront * TR_PSIZE2_BASE * (1 - TR_PSIZE2_DECAY * (double)NFEval / (double)MaxFEval);
        if (psizeval2 <= 1)
            psizeval2 = 2;
        std::discrete_distribution<int> ComponentSelectorFront2 (FitTemp_prand.begin(), FitTemp_prand.begin() + psizeval2);
        std::discrete_distribution<int> ComponentSelectorFront3 (FitTemp_prand.begin(),FitTemp_prand.end());

        int UDEBranchUsed[3] = {0, 0, 0};
        int UDEBranchWin[3] = {0, 0, 0};
        double UDEBranchDelta[3] = {0, 0, 0};
        double UDEBranchFWin[3] = {0, 0, 0};
        double UDEBranchCrWin[3] = {0, 0, 0};
        for(int IndIter=0;IndIter<NIndsFront;IndIter++)
        {
            double maxfitFront = FitArrFront[0];
            for(int i=0;i!=NIndsFront;i++)
                maxfitFront = max(maxfitFront,FitArrFront[i]);
            //get indices for popul
            for(int i=0;i!=NIndsFront;i++)
            {
                if(PenaltyArrFront[i] > epsilon)
                    FitArrCopy[i] = maxfitFront + 1.0 + PenaltyArrFront[i];
                else
                    FitArrCopy[i] = FitArrFront[i];
                if(i == 0)
                {
                    minfit = FitArrCopy[0];
                    maxfit = FitArrCopy[0];
                }
                IndicesFront[i] = i;
                maxfit = max(maxfit,FitArrFront[i]);
                minfit = min(minfit,FitArrFront[i]);
            }
            if(minfit != maxfit)
                qSort2int(FitArrCopy,IndicesFront,0,NIndsFront-1);

            TheChosenOne = IntRandom(NIndsFront);
            MemoryCurrentIndex = IntRandom(MemorySize);
            MemoryCurrentIndex2 = IntRandom(MemorySize+1);
            do
                prand = Indices[IntRandom(psizeval)];
            while(prand == TheChosenOne);
            do
                Rand1 = IndicesFront[ComponentSelectorFront(generator_uni_i_3)];
            while(Rand1 == prand);
            do
                Rand2 = Indices[IntRandom(NIndsFront)];
            while(Rand2 == prand || Rand2 == Rand1);


            double meanF2 = meanF;
            if(Random(0,1) < TR_MEMORY_F_SELECT_PROB*Cvalglobal)
                meanF2 = MemoryF[MemoryCurrentIndex];
            do
                F = NormRand(meanF2,sigmaF);
            while(F < 0.0 || F > 1.0);
            double F2 = F;

            Cr = NormRand(MemoryCr[MemoryCurrentIndex],TR_CR_SIGMA);
            Cr = min(max(Cr,0.0),1.0);
            int WillCrossover = IntRandom(NVars);
            double ActualCr = 0;
            
            int UDELastBranch = 0;
            bool UDETrialEvaluated = false;
            double progress = (double)NFEval / (double)MaxFEval;
            int UDESelectedBranch = 0;
            if(UDEFusionMode == 1 || UDEFusionMode == 2)
            {
                UDESelectedBranch = UDESelectBranch(UDEBranchProb);
            }
            else
            {
                double Rand_EB = Random(0, 1);
                if ((double)NFEval / (double)MaxFEval < TR_EB_DISABLE_BEFORE_PROGRESS){
                    Rand_EB = 2;
                }
                UDESelectedBranch = ((Rand_EB * (1 - NFEval / MaxFEval)) < EB_hybrid_rate) ? 1 : 0;
            }
            UDELastBranch = UDESelectedBranch;
            UDEBranchUsed[UDELastBranch]++;

            if(UDESelectedBranch == 1)
            {
                EB_hybrid_flag[TheChosenOne] = 1;

                do
                    prand = Indices[ComponentSelectorFront2(generator_uni_i_3)];
                while(prand == TheChosenOne);
                do
                    Rand1 = IndicesFront[ComponentSelectorFront3(generator_uni_i_3)];
                while (Rand1 == prand);
                do
                    Rand2 = Indices[ComponentSelectorFront3(generator_uni_i_3)];
                while (Rand2 == prand || Rand2 == Rand1);
                EB_order(prand, Rand1, Rand2);

                if(UDEFusionMode == 2)
                {
                    do
                        F = CachyRand(UDEBranchMF[1], TR_F_CAUCHY_SIGMA);
                    while(F < 0.0 || F > 1.0);
                    Cr = NormRand(UDEBranchMCR[1], TR_CR_SIGMA);
                }
                else
                {
                    do {
                        if (MemoryCurrentIndex2 < MemorySize)
                            F = CachyRand(MemoryF[MemoryCurrentIndex2], TR_F_CAUCHY_SIGMA);
                        else
                            F = CachyRand(TR_EB_F_FALLBACK_MEAN, TR_F_CAUCHY_SIGMA);
                    }  while(F < 0.0 || F > 1.0);

                    if (MemoryCurrentIndex2 < MemorySize) {
                        if (MemoryCr[MemoryCurrentIndex2] < 0)
                            Cr = 0;
                        else
                            Cr = NormRand(MemoryCr[MemoryCurrentIndex2], TR_CR_SIGMA);
                    } else
                        Cr = NormRand(TR_EB_CR_FALLBACK_MEAN, TR_CR_SIGMA);
                }
                if (Cr >= 1)
                    Cr = 1;
                if (Cr <= 0)
                    Cr = 0;

                if ((double)NFEval / (double)MaxFEval < TR_CR_EARLY_PROGRESS_1)
                    Cr = max(Cr, TR_CR_EARLY_MIN_1);
                if ((double)NFEval / (double)MaxFEval < TR_CR_EARLY_PROGRESS_2)
                    Cr = max(Cr, TR_CR_EARLY_MIN_2);

                bool perturbation = rand() / (double)RAND_MAX < TR_PERTURB_PROB;

                for(int j=0;j!=NVars;j++)
                {
                    if(Random(0,1) < Cr || WillCrossover == j)
                    {
                        Trial[j] = PopulFront[TheChosenOne][j] +
                        F * (ord_best_arch[j] - PopulFront[TheChosenOne][j]) +
                        F * (ord_medium_arch[j] - ord_worst_arch[j]);
                        ActualCr++;
                    }
                    else
                    {
                        Trial[j] = perturbation ? CachyRand(PopulFront[TheChosenOne][j], TR_PERTURB_SCALE) : PopulFront[TheChosenOne][j];
                    }
                    if(Trial[j] < Left)
                        Trial[j] = TRBoundRepair(PopulFront[TheChosenOne][j], Left);
                    if(Trial[j] > Right)
                        Trial[j] = TRBoundRepair(PopulFront[TheChosenOne][j], Right);
                    PopulTemp[IndIter][j] = Trial[j];
                }
            }
            else if(UDESelectedBranch == 2)
            {
                EB_hybrid_flag[TheChosenOne] = 0;
                do
                    prand = Indices[IntRandom(NIndsFront)];
                while(prand == TheChosenOne);
                do
                    Rand1 = IndicesFront[IntRandom(NIndsFront)];
                while(Rand1 == TheChosenOne || Rand1 == prand);
                do
                    Rand2 = Indices[IntRandom(NIndsFront)];
                while(Rand2 == prand || Rand2 == Rand1 || Rand2 == TheChosenOne);

                if(UDEFusionMode == 2)
                {
                    do
                        F = CachyRand(UDEBranchMF[2], TR_F_CAUCHY_SIGMA);
                    while(F < 0.0 || F > 1.0);
                    F2 = F;
                    Cr = NormRand(UDEBranchMCR[2], TR_CR_SIGMA);
                    Cr = min(max(Cr,0.0),1.0);
                }
                bool perturbation = rand() / (double)RAND_MAX < TR_PERTURB_PROB;
                for(int j=0;j!=NVars;j++)
                {
                    if(Random(0,1) < Cr || WillCrossover == j)
                    {
                        Trial[j] = PopulFront[TheChosenOne][j] + F*(Popul[prand][j] - PopulFront[TheChosenOne][j]) + F2*(PopulFront[Rand1][j] - Popul[Rand2][j]);
                        ActualCr ++;
                    }
                    else
                    {
                        Trial[j] = perturbation ? CachyRand(PopulFront[TheChosenOne][j], TR_PERTURB_SCALE) : PopulFront[TheChosenOne][j];
                    }
                    if(Trial[j] < Left)
                        Trial[j] = TRBoundRepair(PopulFront[TheChosenOne][j], Left);
                    if(Trial[j] > Right)
                        Trial[j] = TRBoundRepair(PopulFront[TheChosenOne][j], Right);
                    PopulTemp[IndIter][j] = Trial[j];
                }
            }
            else
            {
                EB_hybrid_flag[TheChosenOne] = 0;
                if(UDEFusionMode == 2)
                {
                    do
                        F = NormRand(UDEBranchMF[0], sigmaF);
                    while(F < 0.0 || F > 1.0);
                    F2 = F;
                    Cr = NormRand(UDEBranchMCR[0], TR_CR_SIGMA);
                    Cr = min(max(Cr,0.0),1.0);
                }
                bool perturbation = rand() / (double)RAND_MAX < TR_PERTURB_PROB;
                for(int j=0;j!=NVars;j++)
                {
                    if(Random(0,1) < Cr || WillCrossover == j)
                    {
                        Trial[j] = PopulFront[TheChosenOne][j] + F*(Popul[prand][j] - PopulFront[TheChosenOne][j]) + F2*(PopulFront[Rand1][j] - Popul[Rand2][j]);
                        ActualCr ++;
                    }
                    else
                    {
                        Trial[j] = perturbation ? CachyRand(PopulFront[TheChosenOne][j], TR_PERTURB_SCALE) : PopulFront[TheChosenOne][j];
                    }
                    if(Trial[j] < Left)
                        Trial[j] = TRBoundRepair(PopulFront[TheChosenOne][j], Left);
                    if(Trial[j] > Right)
                        Trial[j] = TRBoundRepair(PopulFront[TheChosenOne][j], Right);
                    PopulTemp[IndIter][j] = Trial[j];
                }
            }

            bool UDEDoRace = false;
            if(UDEFusionMode >= 3 && UDEFusionTopFrac > 0.0 && progress >= UDEFusionProgressMin && progress < UDEFusionProgressMax && FeasRate > TR_FEAS_RATE_MIN && NFEval + 1 < MaxFEval)
            {
                int top_limit = max(1, int(double(NIndsFront) * UDEFusionTopFrac + 0.5));
                for(int rr=0; rr<top_limit; rr++)
                {
                    if(IndicesFront[rr] == TheChosenOne)
                        UDEDoRace = true;
                }
            }
            if(UDEDoRace)
            {
                double bestTrial[100];
                double bestGLocal[3];
                double bestHLocal[6];
                double bestAllLocal[9];
                for(int j=0;j!=NVars;j++)
                    bestTrial[j] = PopulTemp[IndIter][j];
                double bestFitLocal = cec_24_constr(bestTrial, func_num);
                double bestPenaltyLocal = cec_24_totalpenalty(func_num, bestGLocal, bestHLocal, bestAllLocal);
                double bestActualCrLocal = ActualCr;

                double altTrial[100];
                double altGLocal[3];
                double altHLocal[6];
                double altAllLocal[9];
                double altActualCr = 0.0;
                int rA = 0;
                int rB = 0;
                int rC = 0;
                do
                    rA = Indices[IntRandom(NIndsFront)];
                while(rA == TheChosenOne);
                do
                    rB = IndicesFront[IntRandom(NIndsFront)];
                while(rB == rA || rB == TheChosenOne);
                do
                    rC = Indices[IntRandom(NIndsFront)];
                while(rC == rA || rC == rB || rC == TheChosenOne);
                double altF = F;
                if(UDEFusionMode == 2)
                    altF = UDEBranchMF[2];
                double altCr = max(TR_ALT_CR_MIN, min(TR_ALT_CR_MAX, Cr));
                int altPBest = TheChosenOne;
                if(UDEFusionAltFormula == 2)
                {
                    int pbest_limit = max(1, int(double(NIndsFront) * UDEFusionPBestSampleFrac + 0.5));
                    pbest_limit = min(pbest_limit, NIndsFront);
                    int tries = 0;
                    do
                    {
                        altPBest = IndicesFront[IntRandom(pbest_limit)];
                        tries++;
                    }
                    while(altPBest == TheChosenOne && tries < TR_PBEST_RETRY);
                    if(altPBest == TheChosenOne)
                        altPBest = rB;
                }
                if(UDEFusionAltFormula == 3)
                    EB_order(rA, rB, rC);
                int altWillCrossover = IntRandom(NVars);
                for(int j=0;j!=NVars;j++)
                {
                    if(Random(0,1) < altCr || altWillCrossover == j)
                    {
                        if(UDEFusionAltFormula == 1)
                            altTrial[j] = Popul[rA][j] + altF*(PopulFront[rB][j] - Popul[rC][j]);
                        else if(UDEFusionAltFormula == 2)
                            altTrial[j] = PopulFront[TheChosenOne][j] + (TR_PBEST_PULL_SCALE * altF)*(PopulFront[altPBest][j] - PopulFront[TheChosenOne][j]) + (TR_DIFF_SCALE * altF)*(Popul[rA][j] - Popul[rC][j]);
                        else if(UDEFusionAltFormula == 3)
                            altTrial[j] = PopulFront[TheChosenOne][j] + altF*(ord_best_arch[j] - PopulFront[TheChosenOne][j]) + altF*(ord_medium_arch[j] - ord_worst_arch[j]);
                        else
                            altTrial[j] = PopulFront[TheChosenOne][j] + altF*(Popul[rA][j] - PopulFront[TheChosenOne][j]) + altF*(PopulFront[rB][j] - Popul[rC][j]);
                        altActualCr++;
                    }
                    else
                    {
                        altTrial[j] = PopulFront[TheChosenOne][j];
                    }
                    if(altTrial[j] < Left)
                        altTrial[j] = TRBoundRepair(PopulFront[TheChosenOne][j], Left);
                    if(altTrial[j] > Right)
                        altTrial[j] = TRBoundRepair(PopulFront[TheChosenOne][j], Right);
                }
                double altFitLocal = cec_24_constr(altTrial, func_num);
                double altPenaltyLocal = cec_24_totalpenalty(func_num, altGLocal, altHLocal, altAllLocal);
                if(UDEBetterByEpsilon(altFitLocal, altPenaltyLocal, bestFitLocal, bestPenaltyLocal, epsilon))
                {
                    for(int j=0;j!=NVars;j++)
                        bestTrial[j] = altTrial[j];
                    for(int j=0;j!=3;j++)
                        bestGLocal[j] = altGLocal[j];
                    for(int j=0;j!=6;j++)
                        bestHLocal[j] = altHLocal[j];
                    for(int j=0;j!=NConstr;j++)
                        bestAllLocal[j] = altAllLocal[j];
                    bestFitLocal = altFitLocal;
                    bestPenaltyLocal = altPenaltyLocal;
                    bestActualCrLocal = altActualCr;
                    UDELastBranch = 2;
                }

                if(UDEFusionMode == 4 && NFEval + 1 < MaxFEval)
                {
                    double randTrial[100];
                    double randGLocal[3];
                    double randHLocal[6];
                    double randAllLocal[9];
                    double randActualCr = 0.0;
                    int rBase = 0;
                    do
                        rBase = Indices[IntRandom(NIndsFront)];
                    while(rBase == TheChosenOne);
                    do
                        rA = IndicesFront[IntRandom(NIndsFront)];
                    while(rA == rBase || rA == TheChosenOne);
                    do
                        rB = Indices[IntRandom(NIndsFront)];
                    while(rB == rBase || rB == rA || rB == TheChosenOne);
                    int randWillCrossover = IntRandom(NVars);
                    for(int j=0;j!=NVars;j++)
                    {
                        if(Random(0,1) < altCr || randWillCrossover == j)
                        {
                            randTrial[j] = Popul[rBase][j] + altF*(PopulFront[rA][j] - Popul[rB][j]);
                            randActualCr++;
                        }
                        else
                        {
                            randTrial[j] = PopulFront[TheChosenOne][j];
                        }
                        if(randTrial[j] < Left)
                            randTrial[j] = TRBoundRepair(PopulFront[TheChosenOne][j], Left);
                        if(randTrial[j] > Right)
                            randTrial[j] = TRBoundRepair(PopulFront[TheChosenOne][j], Right);
                    }
                    double randFitLocal = cec_24_constr(randTrial, func_num);
                    double randPenaltyLocal = cec_24_totalpenalty(func_num, randGLocal, randHLocal, randAllLocal);
                    if(UDEBetterByEpsilon(randFitLocal, randPenaltyLocal, bestFitLocal, bestPenaltyLocal, epsilon))
                    {
                        for(int j=0;j!=NVars;j++)
                            bestTrial[j] = randTrial[j];
                        for(int j=0;j!=3;j++)
                            bestGLocal[j] = randGLocal[j];
                        for(int j=0;j!=6;j++)
                            bestHLocal[j] = randHLocal[j];
                        for(int j=0;j!=NConstr;j++)
                            bestAllLocal[j] = randAllLocal[j];
                        bestFitLocal = randFitLocal;
                        bestPenaltyLocal = randPenaltyLocal;
                        bestActualCrLocal = randActualCr;
                        UDELastBranch = 2;
                    }
                }

                for(int j=0;j!=NVars;j++)
                {
                    Trial[j] = bestTrial[j];
                    PopulTemp[IndIter][j] = bestTrial[j];
                }
                for(int j=0;j!=3;j++)
                    PopulTempG[IndIter][j] = bestGLocal[j];
                for(int j=0;j!=6;j++)
                    PopulTempH[IndIter][j] = bestHLocal[j];
                for(int j=0;j!=NConstr;j++)
                    PopulAllConstrTemp[IndIter][j] = bestAllLocal[j];
                FitArrTemp[IndIter] = bestFitLocal;
                PenaltyTemp[IndIter] = bestPenaltyLocal;
                ActualCr = bestActualCrLocal;
                UDETrialEvaluated = true;
            }
            ActualCr = ActualCr / double(NVars);
            if(!UDETrialEvaluated)
            {
                FitArrTemp[IndIter] = cec_24_constr(Trial,func_num);
                PenaltyTemp[IndIter] = cec_24_totalpenalty(func_num,PopulTempG[IndIter],PopulTempH[IndIter],PopulAllConstrTemp[IndIter]);
            }
            FitTemp[TheChosenOne] = FitArrTemp[IndIter];
            if(!globalbestinit ||
               ((FitArrTemp[IndIter] <= globalbest) && PenaltyTemp[IndIter] <= 0) ||
                ((PenaltyTemp[IndIter] <= globalbestpenalty) && PenaltyTemp[IndIter] > 0))
            {
                globalbest = FitArrTemp[IndIter];
                globalbestpenalty = PenaltyTemp[IndIter];
                globalbestinit = true;
                bestfit = FitArrTemp[IndIter];
                for(int j=0;j!=NVars;j++)
                    BestInd[j] = PopulTemp[IndIter][j];
                for(int j=0;j!=3;j++)
                    BestG[j] = PopulTempG[IndIter][j];
                for(int j=0;j!=6;j++)
                    BestH[j] = PopulTempH[IndIter][j];
            }

            double improvement = 0;
            bool change = false;
            if(Random(0,1) < 1.0)
            {
                double temppenalty = PenaltyTemp[IndIter];
                double frontpenalty = PenaltyArrFront[TheChosenOne];
                bool goodtemp = false;
                if(temppenalty <= epsilon) //new is epsilon-feasible
                {
                    temppenalty = 0;
                    goodtemp = true;
                }
                if(frontpenalty <= epsilon) // current is epsilon-feasible
                    frontpenalty = 0;
                if(
                   ((goodtemp)&&(FitArrTemp[IndIter]<=FitArrFront[TheChosenOne])) || //new is epsilon-feasible, and with better fitness
                   ((temppenalty==frontpenalty) && (FitArrTemp[IndIter]<=FitArrFront[TheChosenOne])) //new has same penalty and better fitness
                   )
                {
                    change = true;
                    improvement = FitArrFront[TheChosenOne] - FitArrTemp[IndIter];
                }
                else if(temppenalty < frontpenalty) //both are epsilon-infeasible, but new has smaller violation
                {
                    change = true;
                    improvement = frontpenalty - temppenalty;
                }
            }
            else
            {
                ////////////////////////////////////////////////epsilonf
                double temppenalty = PenaltyTemp[IndIter];
                double frontpenalty = PenaltyArrFront[TheChosenOne];
                bool goodtemp = false;
                if(FitArrTemp[IndIter] <= epsilonf)
                {
                    temppenalty = 0;
                    goodtemp = true;
                }
                if(FitArrFront[TheChosenOne] <= epsilonf)
                    frontpenalty = 0;
                if((goodtemp && FitArrTemp[IndIter] <= FitArrFront[TheChosenOne] ) ||
                   ((temppenalty==frontpenalty) && (FitArrTemp[IndIter]<=FitArrFront[TheChosenOne])))
                {
                    change = true;
                    improvement = FitArrFront[TheChosenOne] - FitArrTemp[IndIter];
                }
                else if(temppenalty < frontpenalty)
                {
                    change = true;
                    improvement = frontpenalty - temppenalty;
                }
                ////////////////////////////////////////////////epsilonf
            }
            if(change)
            {
                UDEBranchWin[UDELastBranch]++;
                UDEBranchDelta[UDELastBranch] += fabs(improvement);
                UDEBranchFWin[UDELastBranch] += F;
                UDEBranchCrWin[UDELastBranch] += ActualCr;
                UDEStag[TheChosenOne] = 0;
                for(int j=0;j!=NVars;j++)
                {
                    Popul[NIndsCurrent+SuccessFilled][j] = PopulTemp[IndIter][j];
                    PopulFront[PFIndex][j] = PopulTemp[IndIter][j];
                }
                for(int j=0;j!=3;j++)
                {
                    PopulG[NIndsCurrent+SuccessFilled][j] = PopulTempG[IndIter][j];
                    PopulFrontG[PFIndex][j] = PopulTempG[IndIter][j];
                }
                for(int j=0;j!=6;j++)
                {
                    PopulH[NIndsCurrent+SuccessFilled][j] = PopulTempH[IndIter][j];
                    PopulFrontH[PFIndex][j] = PopulTempH[IndIter][j];
                }
                for(int j=0;j!=NConstr;j++)
                {
                    PopulAllConstr[NIndsCurrent+SuccessFilled][j] = PopulAllConstrTemp[IndIter][j];
                    PopulAllConstrFront[PFIndex][j] = PopulAllConstrTemp[IndIter][j];
                }
                FitArr[NIndsCurrent+SuccessFilled] = FitArrTemp[IndIter];
                FitArrFront[PFIndex] = FitArrTemp[IndIter];
                PenaltyArr[NIndsCurrent+SuccessFilled] = PenaltyTemp[IndIter];
                PenaltyArrFront[PFIndex] = PenaltyTemp[IndIter];
                tempSuccessCr[SuccessFilled] = ActualCr;//Cr
                tempSuccessF[SuccessFilled] = F;
                FitDelta[SuccessFilled] = improvement;
                SuccessFilled++;
                PFIndex = (PFIndex + 1)%NIndsFront;
            }
            else
            {
                if(TheChosenOne >= 0 && TheChosenOne < 600)
                    UDEStag[TheChosenOne]++;
                if(UDEFusionArchiveRecover && UDEArchiveSize < 32 && PenaltyTemp[IndIter] <= max(PenaltyArrFront[TheChosenOne], epsilon + 1e-12))
                {
                    int slot = UDEArchiveCursor % 32;
                    for(int j=0;j!=NVars;j++) UDEArchiveX[slot][j] = PopulTemp[IndIter][j];
                    for(int j=0;j!=3;j++) UDEArchiveG[slot][j] = PopulTempG[IndIter][j];
                    for(int j=0;j!=6;j++) UDEArchiveH[slot][j] = PopulTempH[IndIter][j];
                    for(int j=0;j!=NConstr;j++) UDEArchiveAll[slot][j] = PopulAllConstrTemp[IndIter][j];
                    UDEArchiveFit[slot] = FitArrTemp[IndIter];
                    UDEArchivePenalty[slot] = PenaltyTemp[IndIter];
                    UDEArchiveCursor++;
                    if(UDEArchiveSize < 32) UDEArchiveSize++;
                }
            }
            SaveBestValues(func_num,BestG,BestH);
        }
        if(UDEFusionMode == 1 || UDEFusionMode == 2)
        {
            double score_sum = 0.0;
            double branch_score[3];
            for(int b=0;b!=3;b++)
            {
                branch_score[b] = 0.02;
                if(UDEBranchUsed[b] > 0)
                    branch_score[b] += double(UDEBranchWin[b]) / double(UDEBranchUsed[b]);
                score_sum += branch_score[b];
                if(UDEFusionMode == 2 && UDEBranchWin[b] > 0)
                {
                    UDEBranchMF[b] = 0.7 * UDEBranchMF[b] + 0.3 * (UDEBranchFWin[b] / double(UDEBranchWin[b]));
                    UDEBranchMCR[b] = 0.7 * UDEBranchMCR[b] + 0.3 * (UDEBranchCrWin[b] / double(UDEBranchWin[b]));
                    UDEBranchMF[b] = min(max(UDEBranchMF[b], 0.05), 0.95);
                    UDEBranchMCR[b] = min(max(UDEBranchMCR[b], 0.05), 1.0);
                }
            }
            if(score_sum > 0.0)
            {
                for(int b=0;b!=3;b++)
                    UDEBranchProb[b] = 0.15 * UDEBranchProb[b] + 0.85 * (branch_score[b] / score_sum);
                for(int b=0;b!=3;b++)
                    UDEBranchProb[b] = max(0.05, min(0.90, UDEBranchProb[b]));
            }
        }
        if(UDEFusionArchiveRecover && UDEArchiveSize > 0 && Generation > 20 && FeasRate > 0.45)
        {
            int replace_id = -1;
            for(int i=0;i!=NIndsFront;i++)
            {
                if(UDEStag[i] > 35)
                {
                    replace_id = i;
                    break;
                }
            }
            if(replace_id >= 0)
            {
                int slot = IntRandom(UDEArchiveSize);
                for(int j=0;j!=NVars;j++)
                {
                    PopulFront[replace_id][j] = UDEArchiveX[slot][j];
                    Popul[replace_id][j] = UDEArchiveX[slot][j];
                }
                for(int j=0;j!=3;j++)
                {
                    PopulFrontG[replace_id][j] = UDEArchiveG[slot][j];
                    PopulG[replace_id][j] = UDEArchiveG[slot][j];
                }
                for(int j=0;j!=6;j++)
                {
                    PopulFrontH[replace_id][j] = UDEArchiveH[slot][j];
                    PopulH[replace_id][j] = UDEArchiveH[slot][j];
                }
                for(int j=0;j!=NConstr;j++)
                {
                    PopulAllConstrFront[replace_id][j] = UDEArchiveAll[slot][j];
                    PopulAllConstr[replace_id][j] = UDEArchiveAll[slot][j];
                }
                FitArrFront[replace_id] = UDEArchiveFit[slot];
                FitArr[replace_id] = UDEArchiveFit[slot];
                PenaltyArrFront[replace_id] = UDEArchivePenalty[slot];
                PenaltyArr[replace_id] = UDEArchivePenalty[slot];
                UDEStag[replace_id] = 0;
            }
        }
        for (int ChosenOne = 0; ChosenOne != NIndsFront; ChosenOne++) {
            FitTemp_prand[ChosenOne] = FitTemp[ChosenOne];            
        }
        UpdateEB_hybrid_param(EB_hybrid_flag, FitMass, FitTemp_prand);
        for (int ChosenOne = 0; ChosenOne != NIndsFront; ChosenOne++) {
            FitMass[ChosenOne] = FitArrFront[ChosenOne];            
        }
        SuccessRate = double(SuccessFilled)/double(NIndsFront);
        newNIndsFront = int(double(4-NIndsFrontMax)/double(MaxFEval)*NFEval + NIndsFrontMax);
        RemoveWorst(NIndsFront,newNIndsFront,epsilon);
        NIndsFront = newNIndsFront;
        UpdateMemoryCr();
        NIndsCurrent = NIndsFront + SuccessFilled;
        SuccessFilled = 0;
        Generation++;
        if(NIndsCurrent > NIndsFront)
        {
            double maxfitPop = FitArr[0];
            for(int i=0;i!=NIndsCurrent;i++)
                maxfitPop = max(maxfitPop,FitArr[i]);

            for(int i=0;i!=NIndsCurrent;i++)
            {
                if(PenaltyArr[i] > epsilon)
                    FitArrCopy[i] = maxfitPop + 1.0 + PenaltyArr[i];
                else
                    FitArrCopy[i] = FitArr[i];
                if(i == 0)
                {
                    minfit = FitArrCopy[0];
                    maxfit = FitArrCopy[0];
                }
                Indices[i] = i;
                maxfit = max(maxfit,FitArr[i]);
                minfit = min(minfit,FitArr[i]);
            }
            if(minfit != maxfit)
                qSort2int(FitArrCopy,Indices,0,NIndsCurrent-1);

            NIndsCurrent = NIndsFront;
            for(int i=0;i!=NIndsCurrent;i++)
            {
                for(int j=0;j!=NVars;j++)
                    PopulTemp[i][j] = Popul[Indices[i]][j];
                for(int j=0;j!=3;j++)
                    PopulTempG[i][j] = PopulG[Indices[i]][j];
                for(int j=0;j!=6;j++)
                    PopulTempH[i][j] = PopulH[Indices[i]][j];
                for(int j=0;j!=NConstr;j++)
                    PopulAllConstrTemp[i][j] = PopulAllConstr[Indices[i]][j];
                FitArrTemp[i] = FitArr[Indices[i]];
                PenaltyTemp[i] = PenaltyArr[Indices[i]];
            }
            for(int i=0;i!=NIndsCurrent;i++)
            {
                for(int j=0;j!=NVars;j++)
                    Popul[i][j] = PopulTemp[i][j];
                for(int j=0;j!=3;j++)
                    PopulG[i][j] = PopulTempG[i][j];
                for(int j=0;j!=6;j++)
                    PopulH[i][j] = PopulTempH[i][j];
                for(int j=0;j!=NConstr;j++)
                    PopulAllConstr[i][j] = PopulAllConstrTemp[i][j];
                FitArr[i] = FitArrTemp[i];
                PenaltyArr[i] = PenaltyTemp[i];
            }
        }
    }
}

int main(int argc, char** argv)
{
    unsigned t0g=clock(),t1g;
    int TotalNRuns = 25;
    int FuncStart = 1;
    int FuncEnd = 28;
    if(argc >= 2)
    {
        FuncStart = atoi(argv[1]);
        if(FuncStart < 1) FuncStart = 1;
        if(FuncStart > 28) FuncStart = 28;
    }
    if(argc >= 3)
    {
        FuncEnd = atoi(argv[2]);
        if(FuncEnd < 1) FuncEnd = 1;
        if(FuncEnd > 28) FuncEnd = 28;
    }
    if(FuncStart > FuncEnd)
    {
        int tmp = FuncStart;
        FuncStart = FuncEnd;
        FuncEnd = tmp;
    }
    if(argc >= 4)
    {
        TotalNRuns = atoi(argv[3]);
        if(TotalNRuns < 1) TotalNRuns = 1;
    }
    if(argc >= 5)
    {
        SeedBase_global = unsigned(atoi(argv[4]));
        if(SeedBase_global == 0u) SeedBase_global = 20260611u;
    }

    if(TimeComplexity)
	{
		ofstream fout_t("time_complexity.txt");
		cout<<"Running time complexity code"<<endl;
		double T1, T2;
		unsigned t1=clock(),t0;
		GNVars = 30;
		MaxFEval = 10000;
		double* xtmp = new double[GNVars];
		for(int j=0;j!=GNVars;j++)
			xtmp[j] = 0;

		t0=clock();
		for(int func_num=1;func_num!=29;func_num++)
		{
			for(int j=0;j!=MaxFEval;j++)
			{
				cec_24_constr(xtmp, func_num);
			}
			//cout<<func_num<<endl;
		}
		t1=clock()-t0;
		T1 = double(t1)/28.0;
		cout<<"T1 = "<<T1<<endl;
		fout_t<<"T1 = "<<T1<<endl;

		t0=clock();
		for(int func_num=1;func_num!=29;func_num++)
		{
            globalbestinit = false;
            LastFEcount = 0;
            NFEval = 0;
            Cvalglobal = 0;
            int NewPopSize = 600;
            Optimizer OptZ;
            OptZ.Initialize(NewPopSize, GNVars, func_num);
            OptZ.MainCycle();
            OptZ.Clean();
			//cout<<func_num<<endl;
		}
		t1=clock()-t0;
		T2 = double(t1)/28.0;
		cout<<"T2 = " << T2 << endl;
		fout_t<<"T2 = " << T2 << endl;

		delete xtmp;
	}

    for(int GNVarsIter = 1;GNVarsIter!=2;GNVarsIter++)
    {
        if(GNVarsIter == 0)
            GNVars = 10;
        if(GNVarsIter == 1)
            GNVars = 30;
        if(GNVarsIter == 2)
            GNVars = 50;
        if(GNVarsIter == 3)
            GNVars = 100;
        MaxFEval = GNVars*20000;
        for(int func_num = FuncStart; func_num <= FuncEnd; func_num++)
        {
            sprintf(buffer, "%s_f_F%d_D%d_.txt",algName,func_num,GNVars);
            ofstream fout(buffer);
            sprintf(buffer, "%s_H_F%d_D%d.txt",algName,func_num,GNVars);
            ofstream foutH(buffer);
            sprintf(buffer, "%s_G_F%d_D%d.txt",algName,func_num,GNVars);
            ofstream foutG(buffer);
            sprintf(buffer, "%s_F%d_D%d.txt",algName,func_num,GNVars);
            ofstream foutC3(buffer);
            getOptimum(func_num);
            for (int run = 0;run!=TotalNRuns;run++)
            {
                cout<<"func\t"<<func_num<<"\trun\t"<<run<<endl;
                reset_rng_for_case(func_num, run);
                ResultsArray[ResTsize2-1] = MaxFEval;
                globalbestinit = false;
                LastFEcount = 0;
                NFEval = 0;
                Cvalglobal = 0;
                int NewPopSize = max(4, int(TR_POP_SIZE_FACTOR * GNVars));
                Optimizer OptZ;
                OptZ.Initialize(NewPopSize, GNVars, func_num);
                OptZ.MainCycle();
                OptZ.Clean();

                for(int k=0;k!=ResTsize2;k++)
                    foutC3<<ResultsArray[k]<<"\t";
                foutC3<<endl;

                for(int k=0;k!=ResTsize2;k++)
                {
                    fout<<ResultsArray[k]<<"\t";

                    for(int L=0;L!=6;L++)
                        foutH<<ResultsArrayH[k][L]<<"\t";

                    double sumConstr = 0;
                    for(int L=0;L!=ng_B[func_num-1];L++)
                    {
                        foutG<<ResultsArrayG[k][L]<<"\t";
                        sumConstr += ResultsArrayG[k][L];
                    }
                    if(k == ResTsize2-1)
                        sumConstr = ResultsArray[k];
                    foutC3<<sumConstr<<"\t";
                }

				fout<<endl;
				foutH<<endl;
				foutG<<endl;
				foutC3<<endl;                
            }
			fout.close();
			foutH.close();
			foutG.close();
			foutC3.close();
        }
    }
    t1g=clock()-t0g;
    double T0g = t1g/double(CLOCKS_PER_SEC);
    cout << "Time spent: " << T0g << endl;
	return 0;
}
