//
//  RobotVisionTracker.h
//  RobotVisionTracker
//
//  Created by Iaroslav Omelianenko on 3/30/15.
//  Copyright (c) 2015 Nologin. All rights reserved.
//

#ifndef RobotVisionTracker_RobotVisionTracker_h
#define RobotVisionTracker_RobotVisionTracker_h
#define LOCAL true

#define USE_RF_REGRESSION
//#define USE_REGERESSION

#ifdef LOCAL
#define SAVE_MODELS
#define SAVE_DATA
#include "stdc++.h"
#else
#include <bits/stdc++.h>
#endif

#include <iostream>
#include <sys/time.h>

using namespace std;

#define FOR(i,a,b)  for(int i=(a);i<(b);++i)
#define LL          long long
#define ULL         unsigned long long
#define LD          long double
#define MP          make_pair
#define VC          vector
#define PII         pair <int, int>
#define VI          VC < int >
#define VVI         VC < VI >
#define VVVI        VC < VVI >
#define VPII        VC < PII >
#define VD          VC < double >
#define VVD         VC < VD >
#define VF          VC < float >
#define VVF         VC < VF >
#define VS          VC < string >
#define VVS         VC < VS >
#define VE          VC < Entry >
#define VVE         VC < VE >
#define VB          VC < bool >

template<class T> void print(VC < T > v) {cerr << "[";if (v.size()) cerr << v[0];FOR(i, 1, v.size()) cerr << ", " << v[i];cerr << "]" << endl;}
template<class T> void printWithIndex(VC < T > v) {cerr << "[";if (v.size()) cerr << "0:" <<  v[0];FOR(i, 1, v.size()) cerr << ", " << i << ":" <<  v[i];cerr << "]" << endl;}

#ifdef LOCAL
static bool LOG_DEBUG = true;
#else
static bool LOG_DEBUG = false;
#endif
/*! the message buffer length */
const int kPrintBuffer = 1 << 12;
inline void Printf(const char *fmt, ...) {
    if (LOG_DEBUG) {
        std::string msg(kPrintBuffer, '\0');
        va_list args;
        va_start(args, fmt);
        vsnprintf(&msg[0], kPrintBuffer, fmt, args);
        va_end(args);
        fprintf(stderr, "%s", msg.c_str());
    }
}

inline void Assert(bool exp, const char *fmt, ...) {
    if (!exp) {
        std::string msg(kPrintBuffer, '\0');
        va_list args;
        va_start(args, fmt);
        vsnprintf(&msg[0], kPrintBuffer, fmt, args);
        va_end(args);
        fprintf(stderr, "AssertError:%s\n", msg.c_str());
        exit(-1);
    }
}

inline double getTime() {
    timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + 1e-6 * tv.tv_usec;
}

//
// ----------------------------
//
class HoG{
public:
    HoG();
    int wx;
    int wy;
    int nbin;
    void HOGdescriptor(vector<vector<double>>& Im, vector<double>& descriptor);
    
private:
    double PI = 3.14159;
    void imfilterGx(vector<vector<double>>& Im, vector<vector<double>>& grad_xr);
    void imfilterGy(vector<vector<double>>& Im, vector<vector<double>>& grad_yu);
    void GetAnglesAndMagnit(vector<vector<double>>& grad_yu, vector<vector<double>>& grad_xr, vector<vector<double>>& angles, vector<vector<double>>& magnit);
    void GetVector2Range(vector<vector<double>>& inVec, int L1, int L2, int C1, int C2, vector<vector<double>>& outVec);
    void StraitenVector(vector<vector<double>>& inVec, vector<double>& outLine);
    float L2NormVec1(vector<double>& inVec);
};

// set default values of wx, wy and nbin
HoG::HoG() {
    wx = 5;
    wy = 5;
    nbin = 10;
}


// gradients on x direction
void HoG::imfilterGx(vector<vector<double>>& Im, vector<vector<double>>& grad_xr) {
    // hx = [-1, 0, 1];
    grad_xr.clear();
    grad_xr = Im;
    float TempLeft, TempRight;
    for (int ii = 0; ii < Im.size(); ii++) {
        for (int jj = 0; jj < Im[0].size(); jj++) {
            TempLeft = (jj - 1<0 ? 0 : Im[ii][jj - 1]);
            TempRight = (jj + 1 >= Im[0].size() ? 0 : Im[ii][jj + 1]);
            grad_xr[ii][jj] = TempRight - TempLeft;
        }
    }
    return;
}

// gradients on y direction
void HoG::imfilterGy(vector<vector<double>>& Im, vector<vector<double>>& grad_yu) {
    // hy = [1 0 -1]^T
    grad_yu.clear();
    grad_yu = Im;
    float TempUp, TempDown;
    for (int ii = 0; ii < Im.size(); ii++) {
        for (int jj = 0; jj < Im[0].size(); jj++) {
            TempUp = (ii - 1<0 ? 0 : Im[ii - 1][jj]);
            TempDown = (ii + 1 >= Im.size() ? 0 : Im[ii + 1][jj]);
            grad_yu[ii][jj] = TempUp - TempDown;
        }
    }
    return;
}


// compute angle and magnitude
void HoG::GetAnglesAndMagnit(vector<vector<double>>& grad_yu, vector<vector<double>>& grad_xr,  \
                             vector<vector<double>>& angles, vector<vector<double>>& magnit) {
    angles.clear();
    angles = grad_xr;
    for (int ii = 0; ii < grad_xr.size(); ii++) {
        for (int jj = 0; jj < grad_xr[0].size(); jj++) {
            angles[ii][jj] = atan2(grad_yu[ii][jj], grad_xr[ii][jj]);
        }
    }
    
    magnit.clear();
    magnit = grad_xr;
    for (int ii = 0; ii < grad_xr.size(); ii++) {
        for (int jj = 0; jj < grad_xr[0].size(); jj++) {
            magnit[ii][jj] = sqrt(pow(grad_yu[ii][jj], 2) + pow(grad_xr[ii][jj], 2));
        }
    }
    
    return;
}

void HoG::GetVector2Range(vector<vector<double>>& inVec, int L1, int L2, int C1, int C2, vector<vector<double>>& outVec) {
    outVec.clear();
    int Lnum = L2 - L1 + 1;
    int Cnum = C2 - C1 + 1;
    vector<vector<double>> tmpVec(Lnum, vector<double>(Cnum));
    for (int ii = L1 - 1; ii < L2; ii++) {
        for (int jj = C1 - 1; jj < C2; jj++) {
            tmpVec[ii - L1 + 1][jj - C1 + 1] = inVec[ii][jj];
        }
    }
    outVec = tmpVec;
    
    return;
}

void HoG::StraitenVector(vector<vector<double>>& inVec, vector<double>& outLine) {
    outLine.clear();
    for (int jj = 0; jj < inVec[0].size(); jj++) {
        for (int ii = 0; ii < inVec.size(); ii++) {
            outLine.push_back(inVec[ii][jj]);
        }
    }
    return;
}

float HoG::L2NormVec1(vector<double>& inVec) {
    float value = 0;
    for (int ii = 0; ii < inVec.size(); ii++) {
        value += pow(inVec[ii], 2);
    }
    return sqrt(value);
}


// project angle and magnitude into bins
void HoG::HOGdescriptor(vector<vector<double>>& Im, vector<double>& descriptor) {
    int nwin_x = wx;
    int nwin_y = wy;
    int B = nbin;
    int L = (int)Im.size();
    int C = (int)Im[0].size();
    vector<double> H(nwin_x * nwin_y * B, 0);
    Assert(C > 1, "Error: Input Im has only one column");
    
    int step_x = floor(C / (nwin_x + 1));
    int step_y = floor(L / (nwin_y + 1));
    int cont = 0;
    //    cout << nwin_x << " " << nwin_y << " " << B << " " << L << " " << C << " " << step_x << " " << step_y << endl;
    vector<vector<double>> grad_xr;
    imfilterGx(Im, grad_xr);
    vector<vector<double>> grad_yu;
    imfilterGy(Im, grad_yu);
    vector<vector<double>> angles;
    vector<vector<double>> magnit;
    GetAnglesAndMagnit(grad_yu, grad_xr, angles, magnit);
    
    for (int n = 0; n < nwin_y; n++) {
        for (int m = 0; m < nwin_x; m++) {
            cont++;
            vector<vector<double>> angles2;
            GetVector2Range(angles, n*step_y + 1, (n + 2)*step_y, m*step_x + 1, (m + 2)*step_x, angles2);
            vector<vector<double>> magnit2;
            GetVector2Range(magnit, n*step_y + 1, (n + 2)*step_y, m*step_x + 1, (m + 2)*step_x, magnit2);
            vector<double> v_angles;
            StraitenVector(angles2, v_angles);
            vector<double> v_magnit;
            StraitenVector(magnit2, v_magnit);
            int K = (int)v_angles.size();
            int bin = -1;
            vector<double> H2(B, 0);
            for (float ang_lim = -PI + 2 * PI / B; ang_lim <= PI + 0.01; ang_lim += 2 * PI / B) {
                //cout << ang_lim << "     " << 2*PI/B << endl;
                bin++;
                for (int k = 0; k < K; k++) {
                    if (v_angles[k] < ang_lim) {
                        v_angles[k] = 100;
                        H2[bin] += v_magnit[k];
                    }
                }
            }
            double nH2 = L2NormVec1(H2);
            for (int ss = 0; ss < H2.size(); ss++) {
                H2[ss] = H2[ss] / (nH2 + 0.01);
            }
            
            for (int tt = (cont - 1)*B; tt < cont*B; tt++) {
                H[tt] = H2[tt - (cont - 1)*B];
            }
        }
    }
    
    descriptor.clear();
    descriptor = H;
    
    return;
    
}


//
// ----------------------------
//
struct MT_RNG {
    typedef unsigned int uint32;
    
#define hiBit(u)       ((u) & 0x80000000U)   // mask all but highest   bit of u
#define loBit(u)       ((u) & 0x00000001U)   // mask all but lowest    bit of u
#define loBits(u)      ((u) & 0x7FFFFFFFU)   // mask     the highest   bit of u
#define mixBits(u, v)  (hiBit(u)|loBits(v))  // move hi bit of u to hi bit of v
    
    const uint32 K = 0x9908B0DFU; // a magic constant
    const uint32 N = 624; // length of state vector
    const uint32 M = 397; // a period parameter
    
    uint32   state[624 + 1];     // state vector + 1 extra to not violate ANSI C
    uint32   *next;          // next random value is computed from here
    int      left = -1;      // can *next++ this many times before reloading
    
    
    void seedMT(uint32 seed) {
        uint32 x = (seed | 1U) & 0xFFFFFFFFU, *s = state;
        int    j;
        
        for (left = 0, *s++ = x, j = N; --j;
             *s++ = (x*=69069U) & 0xFFFFFFFFU);
    }
    
    
    uint32 reloadMT(void) {
        uint32 *p0 = state, *p2 = state + 2, *pM = state + M, s0, s1;
        int    j;
        
        if (left < -1)
            seedMT(4357U);
        
        left = N - 1, next = state + 1;
        
        for (s0 = state[0], s1 = state[1], j = N - M + 1; --j; s0 = s1, s1 = *p2++)
            *p0++ = *pM++ ^ (mixBits(s0, s1) >> 1) ^ (loBit(s1) ? K : 0U);
        
        for (pM = state, j = M; --j; s0 = s1, s1 = *p2++)
            *p0++ = *pM++ ^ (mixBits(s0, s1) >> 1) ^ (loBit(s1) ? K : 0U);
        
        s1=state[0], *p0 = *pM ^ (mixBits(s0, s1) >> 1) ^ (loBit(s1) ? K : 0U);
        s1 ^= (s1 >> 11);
        s1 ^= (s1 <<  7) & 0x9D2C5680U;
        s1 ^= (s1 << 15) & 0xEFC60000U;
        return(s1 ^ (s1 >> 18));
    }
    
    
    uint32 randomMT() {
        uint32 y;
        
        if(--left < 0)
            return(reloadMT());
        
        y  = *next++;
        y ^= (y >> 11);
        y ^= (y <<  7) & 0x9D2C5680U;
        y ^= (y << 15) & 0xEFC60000U;
        y ^= (y >> 18);
        return(y);
    }
    inline int nextInt() {
        return randomMT();
    }
    
    inline int nextInt(int x) {
        return randomMT() % x;
    }
    
    inline int nextInt(int a, int b) {
        return a + (randomMT() % (b - a));
    }
    
#define MAX_UINT_COKUS 4294967295  //basically 2^32-1
    double unif_rand(){
        return (((double)randomMT())/((double)MAX_UINT_COKUS));
    }
};

#define qsort_Index
void R_qsort_I(double *v, int *I, int i, int j) {
    
    int il[31], iu[31];
    double vt, vtt;
    double R = 0.375;
    int ii, ij, k, l, m;
#ifdef qsort_Index
    int it, tt;
#endif
    
    --v;
#ifdef qsort_Index
    --I;
#endif
    
    ii = i;/* save */
    m = 1;
    
L10:
    if (i < j) {
        if (R < 0.5898437) R += 0.0390625; else R -= 0.21875;
    L20:
        k = i;
        /* ij = (j + i) >> 1; midpoint */
        ij = i + (int)((j - i)*R);
#ifdef qsort_Index
        it = I[ij];
#endif
        vt = v[ij];
        if (v[i] > vt) {
#ifdef qsort_Index
            I[ij] = I[i]; I[i] = it; it = I[ij];
#endif
            v[ij] = v[i]; v[i] = vt; vt = v[ij];
        }
        /* L30:*/
        l = j;
        if (v[j] < vt) {
#ifdef qsort_Index
            I[ij] = I[j]; I[j] = it; it = I[ij];
#endif
            v[ij] = v[j]; v[j] = vt; vt = v[ij];
            if (v[i] > vt) {
#ifdef qsort_Index
                I[ij] = I[i]; I[i] = it; it = I[ij];
#endif
                v[ij] = v[i]; v[i] = vt; vt = v[ij];
            }
        }
        
        for(;;) { /*L50:*/
            //do l--;  while (v[l] > vt);
            l--;for(;v[l]>vt;l--);
            
            
#ifdef qsort_Index
            tt = I[l];
#endif
            vtt = v[l];
            /*L60:*/
            //do k++;  while (v[k] < vt);
            k=k+1;for(;v[k]<vt;k++);
            
            if (k > l) break;
            
            /* else (k <= l) : */
#ifdef qsort_Index
            I[l] = I[k]; I[k] =  tt;
#endif
            v[l] = v[k]; v[k] = vtt;
        }
        
        m++;
        if (l - i <= j - k) {
            /*L70: */
            il[m] = k;
            iu[m] = j;
            j = l;
        }
        else {
            il[m] = i;
            iu[m] = l;
            i = k;
        }
    }else { /* i >= j : */
        
    L80:
        if (m == 1)     return;
        
        /* else */
        i = il[m];
        j = iu[m];
        m--;
    }
    
    if (j - i > 10)     goto L20;
    
    if (i == ii)        goto L10;
    
    --i;
L100:
    do {
        ++i;
        if (i == j) {
            goto L80;
        }
#ifdef qsort_Index
        it = I[i + 1];
#endif
        vt = v[i + 1];
    } while (v[i] <= vt);
    
    k = i;
    
    do { /*L110:*/
#ifdef qsort_Index
        I[k + 1] = I[k];
#endif
        v[k + 1] = v[k];
        --k;
    } while (vt < v[k]);
    
#ifdef qsort_Index
    I[k + 1] = it;
#endif
    v[k + 1] = vt;
    goto L100;
}



struct RF_config {
    // number of trees in run.  200-500 gives pretty good results
    int nTree = 500;
    // number of variables to pick to split on at each node.  mdim/3 seems to give genrally good performance, but it can be altered up or down
    int mtry;
    
    // 0 or 1 (default is 1) sampling with or without replacement
    bool replace = true;
    // Minimum size of terminal nodes. Setting this number larger causes smaller trees to be grown (and thus take less time). Note that
    // the default values are different for classification (1) and regression (5).
    int nodesize = 5;
    // Should importance of predictors be assessed?
    bool importance = false;
    // Should casewise importance measure be computed? (Setting this to TRUE will override importance.)
    bool localImp = false;
    
    // Should proximity measure among the rows be calculated?
    bool proximity = false;
    // Should proximity be calculated only on 'out-of-bag' data?
    bool oob_prox = false;
    // Should an n by ntree matrix be returned that keeps track of which samples are 'in-bag' in which trees (but not how many times, if sampling with replacement)
    bool keep_inbag = false;
    // If set to TRUE, give a more verbose output as randomForest is run. If set to some integer, then running output is printed for every do_trace trees.
    int do_trace = 1;
    
    // which happens only for regression. perform bias correction for regression? Note: Experimental.risk.
    bool corr_bias = false;
    // Number of times the OOB data are permuted per tree for assessing variable
    // importance. Number larger than 1 gives slightly more stable estimate, but not
    // very effective. Currently only implemented for regression.
    int nPerm = 1;
    // a 1xD true/false vector to say which features are categorical (true), which are numeric (false)
    // maximum of 32 categories per feature is permitted
    VB categorical_feature;
    
    // whether to run run test data prediction during training against current tree
    bool testdat = false;//true;
    // the number of test trees for test data predicitions
    int nts = 10;
    // controls whether to save test set MSE labels
    bool labelts = true;
};

#define swapInt(a, b) ((a ^= b), (b ^= a), (a ^= b))

#if !defined(ARRAY_SIZE)
#define ARRAY_SIZE(x) (sizeof((x)) / sizeof((x)[0]))
#endif

/**
 * Do random forest regression.
 */
class RF_Regression {
    
    typedef enum {
        NODE_TERMINAL = -1,
        NODE_TOSPLIT  = -2,
        NODE_INTERIOR = -3
    } NodeStatus;
    
    typedef char small_int;
    
    
    // the random number generator
    MT_RNG rnd;
    
    //Global to  handle mem in findBestSplit
    int in_findBestSplit = 0; // 0 -initialize and normal.  1-normal  , -99 release
    int in_regTree = 0; //// 0 -initialize and normal.  1-normal  , -99 release
    
    //
    // the model definitions
    //
    /*  a matrix with nclass + 2 (for classification) or two (for regression) columns.
     For classification, the first nclass columns are the class-specific measures
     computed as mean decrease in accuracy. The nclass + 1st column is the
     mean decrease in accuracy over all classes. The last column is the mean decrease
     in Gini index. For Regression, the first column is the mean decrease in
     accuracy and the second the mean decrease in MSE. If importance=FALSE,
     the last measure is still returned as a vector. */
    double *impout = NULL;
    /*  The 'standard errors' of the permutation-based importance measure. For classification,
     a p by nclass + 1 matrix corresponding to the first nclass + 1
     columns of the importance matrix. For regression, a length p vector. */
    double *impSD = NULL;
    /*  a p by n matrix containing the casewise importance measures, the [i,j] element
     of which is the importance of i-th variable on the j-th case. NULL if
     localImp=FALSE. */
    double *impmat = NULL;
    // number of trees grown.
    int ntree;
    // number of predictors sampled for spliting at each node.
    int mtry;
    // the number of nodes to be created
    int nrnodes;
    // vector of mean square errors: sum of squared residuals divided by n.
    double *mse = NULL;
    // number of times cases are 'out-of-bag' (and thus used in computing OOB error estimate)
    int *nout = NULL;
    /*  if proximity=TRUE when randomForest is called, a matrix of proximity
     measures among the input (based on the frequency that pairs of data points are
     in the same terminal nodes). */
    double *prox = NULL;
    
    int *ndtree = NULL;
    small_int *nodestatus = NULL;
    int *lDaughter = NULL;
    int *rDaughter = NULL;
    double *avnode = NULL;
    int *mbest = NULL;
    double *upper = NULL;
    int *inbag = NULL;
    double *coef = NULL;
    double *y_pred_trn = NULL;
    
    // the number of categories per feature if any
    int *ncat = NULL;
    // the maximal number of categories in any feature
    int maxcat;
    // the original uniques per feature
    int **orig_uniques_in_feature = NULL;
    
public:
    
    void train(const VVD &input_X, const VD &input_Y, const RF_config &config) {
        int n_size = (int)input_X.size(); // rows
        int p_size = (int)input_X[0].size(); // cols
        
        int sampsize = n_size;
        int nodesize = config.nodesize;
        int nsum = sampsize;
        nrnodes = 2 * (int)((float)floor((float)(sampsize / ( 1 > (nodesize - 4) ? 1 : (nodesize - 4))))) + 1;
        ntree = config.nTree;
        
        Printf("sampsize: %d, nodesize: %d, nsum %d, nrnodes %d\n", sampsize, nodesize, nsum, nrnodes);
        Printf("doprox: %i, oobProx %i, biascorr %i\n", config.proximity, config.oob_prox, config.corr_bias);
        
        Assert(sampsize == input_Y.size(), "Number of samples must be equal to number of observations");
        Assert(config.mtry > 0, "Please specify number of variables to pick to split on at each node.");
        
        mtry = config.mtry;
        
        // prepare categorical inputs
        ncat = (int*) calloc(p_size, sizeof(int));
        if (config.categorical_feature.size() > 0) {
            Assert(config.categorical_feature.size() == p_size, "If provided, than list of categorical features marks must have size equal to the features dimension");
            orig_uniques_in_feature = (int **)malloc(p_size * sizeof(int *));
            for (int i = 0; i < p_size; i++) {
                if (config.categorical_feature[i]) {
                    // map categorical features
                    ncat[i] = findSortedUniqueFeaturesAndMap(input_X, i, orig_uniques_in_feature[i]);
                } else {
                    // just numerical value
                    ncat[i] = 1;
                }
            }
        } else {
            // all features numerical - set all values just to ones
            for (int i = 0; i < p_size; i++) ncat[i] = 1;
        }
        // find max of categroies
        maxcat = 1;
        for (int i = 0; i < p_size; i++) {
            maxcat = max(maxcat, ncat[i]);
        }
        
        //double y_pred_trn[n_size];
        y_pred_trn = (double*) calloc(n_size, sizeof(double));
        
        
        int imp[] = {config.importance, config.localImp, config.nPerm};
        if (imp[0] == 1) {
            impout = (double*) calloc(p_size * 2, sizeof(double));
        } else {
            impout = (double*) calloc(p_size, sizeof(double));
        }
        if (imp[1] == 1) {
            impmat = (double*) calloc(p_size * n_size, sizeof(double));
        } else {
            impmat = (double*) calloc(1, sizeof(double));
            impmat[0] = 0;
        }
        if (imp[0] == 1) {
            impSD = (double*)calloc(p_size, sizeof(double));
        } else {
            impSD = (double*)calloc(1, sizeof(double));
            impSD[0]=0;
        }
        
        // Should an n by ntree matrix be returned that keeps track of which samples are 'in-bag' in which trees (but not how many times, if sampling with replacement)
        int keepf[2];
        keepf[0] = 1;
        keepf[1] = config.keep_inbag;
        int nt;
        if (keepf[0] == 1){
            nt = ntree;
        } else {
            nt = 1;
        }
        
        // create ouput proximity matrix
        if (!config.proximity) {
            prox = (double*)calloc(1, sizeof(double));
            prox[0] = 0;
        } else {
            prox = (double*)calloc(n_size * n_size, sizeof(double));
        }
        
        //int ndtree[ntree];
        ndtree = (int*)calloc(ntree, sizeof(int));
        
        //int nodestatus[nrnodes * nt];
        nodestatus = (small_int*)calloc(nrnodes*nt, sizeof(small_int));
        
        //int lDaughter[nrnodes * nt];
        lDaughter = (int*)calloc(nrnodes*nt, sizeof(int));
        
        //int rDaughter[nrnodes * nt];
        rDaughter = (int*)calloc(nrnodes*nt, sizeof(int));
        
        //double avnode[nrnodes * nt];
        avnode = (double*) calloc(nrnodes*nt, sizeof(double));
        
        //int mbest[nrnodes * nt];
        mbest=(int*)calloc(nrnodes*nt, sizeof(int));
        
        //double upper[nrnodes * nt];
        upper = (double*) calloc(nrnodes*nt, sizeof(double));
        // vector of mean square errors: sum of squared residuals divided by n.
        mse = (double*)calloc(ntree, sizeof(double));
        
        // copy data
        //        double X[n_size * p_size], Y[n_size];
        double Y[n_size];
        
        // allocate on heap
        double *X = (double *) calloc(n_size * p_size, sizeof(double));
        
        int dimx[2];
        dimx[0] = n_size;
        dimx[1] = p_size;
        
        for (int i = 0; i < n_size; i++) {
            for (int j = 0; j < p_size; j++){
                if (ncat[j] == 1) {
                    // just ordinary numeric feature
                    double fval = input_X[i][j];
                    X[i * p_size + j] = fval;
                } else {
                    // store mapped value
                    int val = input_X[i][j];
                    for (int k = 0; k < ARRAY_SIZE(orig_uniques_in_feature[j]); k++) {
                        if (val == orig_uniques_in_feature[j][k]) {
                            val = k;
                            break;
                        }
                    }
                    X[i * p_size + j] = val;
                }
            }
            Y[i] = input_Y[i];
        }
        
        int replace = config.replace;
        int testdat = config.testdat;
        int nts = config.nts;
        
        double *xts = X;
        double *yts = Y;
        int labelts = config.labelts;
        
        //double yTestPred[nts];
        double *yTestPred; yTestPred = (double*)calloc(nts, sizeof(double));
        double proxts[] = {1};
        
        double *msets;
        if (labelts == 1) {
            msets = (double*)calloc(ntree, sizeof(double));
        } else {
            msets = (double*)calloc(ntree, sizeof(double));
            msets[0] = 1;
        }
        
        coef = (double*)calloc(2, sizeof(double));
        
        //int nout[n_size];
        nout = (int*)calloc(n_size, sizeof(int));
        
        if (keepf[1] == 1) {
            inbag = (int*)calloc(n_size * ntree, sizeof(int));
        } else {
            inbag = (int*)calloc(1, sizeof(int));
            inbag[0] = 1;
        }
        
        int jprint = config.do_trace;
        bool print_verbose_tree_progression = false;
        
        //train the RF
        regRF(X, Y, dimx, &sampsize,
              &nodesize, &nrnodes, &ntree, &mtry,
              imp, ncat, maxcat, &jprint,
              config.proximity, config.oob_prox, config.corr_bias, y_pred_trn,
              impout, impmat, impSD, prox,
              ndtree, nodestatus, lDaughter, rDaughter,
              avnode, mbest, upper, mse,
              keepf, &replace, testdat, xts,
              &nts, yts, labelts, yTestPred,
              proxts, msets, coef, nout,
              inbag, print_verbose_tree_progression) ;
        
        // let the train variables go free
        free(yTestPred);
        free(msets);
        free(X);
    }
    
    VD predict(const VVD &test_X, const RF_config &config) {
        int n_size = (int)test_X.size(); // rows
        int p_size = (int)test_X[0].size(); // cols
        double* ypred = (double*)calloc(n_size, sizeof(double));
        int mdim = p_size;
        
        double* xsplit = upper;
        double* avnodes = avnode;
        int* treeSize = ndtree;
        int keepPred = 0;
        double allPred = 0;
        int nodes = 0;
        int *nodex; nodex = (int*)calloc(n_size, sizeof(int));
        
        double* proxMat;
        if (!config.proximity) {
            proxMat = (double*)calloc(1, sizeof(double));
            proxMat[0] = 0;
        } else {
            proxMat = (double*)calloc(n_size * n_size, sizeof(double));
        }
        
        double X_test[n_size * p_size];
        for (int i = 0; i < n_size; i++) {
            for (int j = 0; j < p_size; j++){
                if (ncat[j] == 1) {
                    // just ordinary numeric feature
                    X_test[i * p_size + j] = test_X[i][j];
                } else {
                    // store mapped value
                    int val = test_X[i][j];
                    for (int k = 0; k < ARRAY_SIZE(orig_uniques_in_feature[j]); k++) {
                        if (val == orig_uniques_in_feature[j][k]) {
                            val = k;
                            break;
                        }
                    }
                    X_test[i * p_size + j] = val;
                }
            }
        }
        
        regForest(X_test, ypred, &mdim, &n_size,
                  &ntree, lDaughter, rDaughter,
                  nodestatus, &nrnodes, xsplit,
                  avnodes, mbest, treeSize, ncat,
                  maxcat, &keepPred, &allPred, config.proximity,
                  proxMat, &nodes, nodex);
        
        VD res(n_size, 0);
        for (int i = 0;i < n_size;i++) {
            res[i] = ypred[i];
        }
        
        free(ypred);
        free(nodex);
        free(proxMat);
        
        return res;
    }
    
    /**
     * Invoked to clear stored model state
     */
    void release() {
        // let the model variables go free
        free(nout);
        free(inbag);
        free(y_pred_trn);
        free(impout);
        free(impmat);
        free(impSD);
        free(mse);
        free(ndtree);
        free(nodestatus);
        free(lDaughter);
        free(rDaughter);
        free(upper);
        free(avnode);
        free(mbest);
        free(ncat);
        
        if (orig_uniques_in_feature) {
            int N = ARRAY_SIZE(orig_uniques_in_feature);
            for(int i = 0; i < N; i++) {
                free(orig_uniques_in_feature[i]);
            }
            free(orig_uniques_in_feature);
        }
    }
    
private:
    
    inline int findSortedUniqueFeaturesAndMap(const VVD input_x, const int fIndex, int *features) const {
        size_t rows = input_x.size();
        VD fTmp(rows, 0);
        for (int i = 0; i < rows; i++) {
            fTmp[i] = input_x[i][fIndex];
        }
        sort(fTmp.begin(), fTmp.end());
        VD unique;
        int previous = numeric_limits<int>::min();
        for (int i = 0; i < rows; i++) {
            if (fTmp[i] != previous) {
                previous = fTmp[i];
                unique.push_back(fTmp[i]);
            }
        }
        int catNum = (int)unique.size();
        features = (int *)malloc(catNum * sizeof(int));
        
        for (int i = 0; i < catNum; i++) {
            features[i] = (int)unique[i];
        }
        return catNum;
    }
    
    void regRF(double *x, double *y, int *xdim, int *sampsize,
               int *nthsize, int *nrnodes, int *nTree, int *mtry, int *imp,
               int *cat, int maxcat, int *jprint, int doProx, int oobprox,
               int biasCorr, double *yptr, double *errimp, double *impmat,
               double *impSD, double *prox, int *treeSize, small_int *nodestatus,
               int *lDaughter, int *rDaughter, double *avnode, int *mbest,
               double *upper, double *mse, const int *keepf, int *replace,
               int testdat, double *xts, int *nts, double *yts, int labelts,
               double *yTestPred, double *proxts, double *msets, double *coef,
               int *nout, int *inbag, int print_verbose_tree_progression) {
        
        
        double errts = 0.0, averrb, meanY, meanYts, varY, varYts, r, xrand,
        errb = 0.0, resid=0.0, ooberr, ooberrperm, delta, *resOOB;
        
        double *yb, *xtmp, *xb, *ytr, *ytree, *tgini;
        
        int k, m, mr, n, nOOB, j, jout, idx, ntest, last, ktmp, nPerm, nsample, mdim, keepF, keepInbag;
        int *oobpair, varImp, localImp, *varUsed;
        
        int *in, *nind, *nodex, *nodexts;
        
        //Abhi:temp variable
        double tmp_d;
        int tmp_i;
        small_int tmp_c;
        
        //Do initialization for COKUS's Random generator
        rnd.seedMT(2*rand()+1);  //works well with odd number so why don't use that
        
        nsample = xdim[0];
        mdim = xdim[1];
        ntest = *nts;
        varImp = imp[0];
        localImp = imp[1];
        nPerm = imp[2]; //printf("nPerm %d\n",nPerm);
        keepF = keepf[0];
        keepInbag = keepf[1];
        
        if (*jprint == 0) *jprint = *nTree + 1;
        
        yb         = (double *) calloc(*sampsize, sizeof(double));
        xb         = (double *) calloc(mdim * *sampsize, sizeof(double));
        ytr        = (double *) calloc(nsample, sizeof(double));
        xtmp       = (double *) calloc(nsample, sizeof(double));
        resOOB     = (double *) calloc(nsample, sizeof(double));
        in        = (int *) calloc(nsample, sizeof(int));
        nodex      = (int *) calloc(nsample, sizeof(int));
        varUsed    = (int *) calloc(mdim, sizeof(int));
        nind = *replace ? NULL : (int *) calloc(nsample, sizeof(int));
        
        oobpair = (doProx && oobprox) ?
        (int *) calloc(nsample * nsample, sizeof(int)) : NULL;
        
        /* If variable importance is requested, tgini points to the second
         "column" of errimp, otherwise it's just the same as errimp. */
        tgini = varImp ? errimp + mdim : errimp;
        
        averrb = 0.0;
        meanY = 0.0;
        varY = 0.0;
        
        zeroDouble(yptr, nsample);
        zeroInt(nout, nsample);
        for (n = 0; n < nsample; ++n) {
            varY += n * (y[n] - meanY) * (y[n] - meanY) / (n + 1);
            meanY = (n * meanY + y[n]) / (n + 1);
        }
        varY /= nsample;
        
        varYts = 0.0;
        meanYts = 0.0;
        if (testdat) {
            for (n = 0; n <= ntest; ++n) {
                varYts += n * (yts[n] - meanYts) * (yts[n] - meanYts) / (n + 1);
                meanYts = (n * meanYts + yts[n]) / (n + 1);
            }
            varYts /= ntest;
        }
        
        if (doProx) {
            zeroDouble(prox, nsample * nsample);
            if (testdat) zeroDouble(proxts, ntest * (nsample + ntest));
        }
        
        if (varImp) {
            zeroDouble(errimp, mdim * 2);
            if (localImp) zeroDouble(impmat, nsample * mdim);
        } else {
            zeroDouble(errimp, mdim);
        }
        if (labelts) zeroDouble(yTestPred, ntest);
        
        /* print header for running output */
        if (*jprint <= *nTree) {
            Printf("     |      Out-of-bag   ");
            if (testdat) Printf("|       Test set    ");
            Printf("|\n");
            Printf("Tree |      MSE  %%Var(y) ");
            if (testdat) Printf("|      MSE  %%Var(y) ");
            Printf("|\n");
        }
        /*************************************
         * Start the loop over trees.
         *************************************/
        
        time_t curr_time;
        if (testdat) {
            ytree = (double *) calloc(ntest, sizeof(double));
            nodexts = (int *) calloc(ntest, sizeof(int));
        }
        
        for (j = 0; j < *nTree; ++j) {
            
            idx = keepF ? j * *nrnodes : 0;
            zeroInt(in, nsample);
            zeroInt(varUsed, mdim);
            /* Draw a random sample for growing a tree. */
            
            if (*replace) { /* sampling with replacement */
                for (n = 0; n < *sampsize; ++n) {
                    xrand = rnd.unif_rand();
                    k = xrand * nsample;
                    in[k] = 1;
                    yb[n] = y[k];
                    for(m = 0; m < mdim; ++m) {
                        xb[m + n * mdim] = x[m + k * mdim];
                    }
                }
            } else { /* sampling w/o replacement */
                for (n = 0; n < nsample; ++n) nind[n] = n;
                last = nsample - 1;
                for (n = 0; n < *sampsize; ++n) {
                    ktmp = (int) (rnd.unif_rand() * (last+1));
                    k = nind[ktmp];
                    swapInt(nind[ktmp], nind[last]);
                    last--;
                    in[k] = 1;
                    yb[n] = y[k];
                    for(m = 0; m < mdim; ++m) {
                        xb[m + n * mdim] = x[m + k * mdim];
                    }
                }
            }
            
            if (keepInbag) {
                for (n = 0; n < nsample; ++n) inbag[n + j * nsample] = in[n];
            }
            
            /* grow the regression tree */
            regTree(xb, yb, mdim, *sampsize, lDaughter + idx, rDaughter + idx,
                    upper + idx, avnode + idx, nodestatus + idx, *nrnodes,
                    treeSize + j, *nthsize, *mtry, mbest + idx, cat, tgini,
                    varUsed);
            
            /* predict the OOB data with the current tree */
            /* ytr is the prediction on OOB data by the current tree */
            predictRegTree(x, nsample, mdim, lDaughter + idx,
                           rDaughter + idx, nodestatus + idx, ytr, upper + idx,
                           avnode + idx, mbest + idx, treeSize[j], cat, maxcat,
                           nodex);
            /* yptr is the aggregated prediction by all trees grown so far */
            errb = 0.0;
            ooberr = 0.0;
            jout = 0; /* jout is the number of cases that has been OOB so far */
            nOOB = 0; /* nOOB is the number of OOB samples for this tree */
            for (n = 0; n < nsample; ++n) {
                if (in[n] == 0) {
                    nout[n]++;
                    nOOB++;
                    yptr[n] = ((nout[n]-1) * yptr[n] + ytr[n]) / nout[n];
                    resOOB[n] = ytr[n] - y[n];
                    ooberr += resOOB[n] * resOOB[n];
                }
                if (nout[n]) {
                    jout++;
                    errb += (y[n] - yptr[n]) * (y[n] - yptr[n]);
                }
            }
            errb /= jout;
            /* Do simple linear regression of y on yhat for bias correction. */
            if (biasCorr) simpleLinReg(nsample, yptr, y, coef, &errb, nout);
            
            /* predict testset data with the current tree */
            if (testdat) {
                predictRegTree(xts, ntest, mdim, lDaughter + idx,
                               rDaughter + idx, nodestatus + idx, ytree,
                               upper + idx, avnode + idx,
                               mbest + idx, treeSize[j], cat, maxcat, nodexts);
                /* ytree is the prediction for test data by the current tree */
                /* yTestPred is the average prediction by all trees grown so far */
                errts = 0.0;
                for (n = 0; n < ntest; ++n) {
                    yTestPred[n] = (j * yTestPred[n] + ytree[n]) / (j + 1);
                }
                /* compute testset MSE */
                if (labelts) {
                    for (n = 0; n < ntest; ++n) {
                        resid = biasCorr ? yts[n] - (coef[0] + coef[1] * yTestPred[n]) : yts[n] - yTestPred[n];
                        errts += resid * resid;
                    }
                    errts /= ntest;
                }
            }
            
            /* Print running output. */
            if ((j + 1) % *jprint == 0) {
                Printf("%4d |", j + 1);
                Printf(" %8.4g %8.2f ", errb, 100 * errb / varY);
                if(labelts == 1)
                    Printf("| %8.4g %8.2f ", errts, 100.0 * errts / varYts);
                Printf("|\n");
            }
            
            mse[j] = errb;
            if (labelts) msets[j] = errts;
            
            /*  DO PROXIMITIES */
            if (doProx) {
                computeProximity(prox, oobprox, nodex, in, oobpair, nsample);
                /* proximity for test data */
                if (testdat) {
                    /* In the next call, in and oobpair are not used. */
                    computeProximity(proxts, 0, nodexts, in, oobpair, ntest);
                    for (n = 0; n < ntest; ++n) {
                        for (k = 0; k < nsample; ++k) {
                            if (nodexts[n] == nodex[k]) {
                                proxts[n + ntest * (k+ntest)] += 1.0;
                            }
                        }
                    }
                }
            }
            
            /* Variable importance */
            if (varImp) {
                for (mr = 0; mr < mdim; ++mr) {
                    if (varUsed[mr]) { /* Go ahead if the variable is used */
                        /* make a copy of the m-th variable into xtmp */
                        for (n = 0; n < nsample; ++n)
                            xtmp[n] = x[mr + n * mdim];
                        ooberrperm = 0.0;
                        for (k = 0; k < nPerm; ++k) {
                            permuteOOB(mr, x, in, nsample, mdim);
                            predictRegTree(x, nsample, mdim, lDaughter + idx,
                                           rDaughter + idx, nodestatus + idx, ytr,
                                           upper + idx, avnode + idx, mbest + idx,
                                           treeSize[j], cat, maxcat, nodex);
                            for (n = 0; n < nsample; ++n) {
                                if (in[n] == 0) {
                                    r = ytr[n] - y[n];
                                    ooberrperm += r * r;
                                    if (localImp) {
                                        impmat[mr + n * mdim] += (r * r - resOOB[n] * resOOB[n]) / nPerm;
                                    }
                                }
                            }
                        }
                        delta = (ooberrperm / nPerm - ooberr) / nOOB;
                        errimp[mr] += delta;
                        impSD[mr] += delta * delta;
                        /* copy original data back */
                        for (n = 0; n < nsample; ++n)
                            x[mr + n * mdim] = xtmp[n];
                    }
                    
                }
                
            }
        }
        /* end of tree iterations=======================================*/
        
        if (biasCorr) {  /* bias correction for predicted values */
            for (n = 0; n < nsample; ++n) {
                if (nout[n]) yptr[n] = coef[0] + coef[1] * yptr[n];
            }
            if (testdat) {
                for (n = 0; n < ntest; ++n) {
                    yTestPred[n] = coef[0] + coef[1] * yTestPred[n];
                }
            }
        }
        
        if (doProx) {
            for (n = 0; n < nsample; ++n) {
                for (k = n + 1; k < nsample; ++k) {
                    prox[nsample*k + n] /= oobprox ?
                    (oobpair[nsample*k + n] > 0 ? oobpair[nsample*k + n] : 1) :
                    *nTree;
                    prox[nsample * n + k] = prox[nsample * k + n];
                }
                prox[nsample * n + n] = 1.0;
            }
            if (testdat) {
                for (n = 0; n < ntest; ++n)
                    for (k = 0; k < ntest + nsample; ++k)
                        proxts[ntest*k + n] /= *nTree;
            }
        }
        
        if (varImp) {
            for (m = 0; m < mdim; ++m) {
                errimp[m] = errimp[m] / *nTree;
                impSD[m] = sqrt( ((impSD[m] / *nTree) -
                                  (errimp[m] * errimp[m])) / *nTree );
                if (localImp) {
                    for (n = 0; n < nsample; ++n) {
                        impmat[m + n * mdim] /= nout[n];
                    }
                }
            }
        }
        for (m = 0; m < mdim; ++m) tgini[m] /= *nTree;
        
        in_findBestSplit=-99;
        findBestSplit(&tmp_d, &tmp_i, &tmp_d, tmp_i, tmp_i,
                      tmp_i, tmp_i, &tmp_i, &tmp_d,
                      &tmp_d, &tmp_i, &tmp_i, tmp_i,
                      tmp_d, tmp_i, &tmp_i);
        
        //do the same mxFreeing of space by calling with -99
        in_regTree=-99;
        regTree(&tmp_d, &tmp_d, tmp_i, tmp_i, &tmp_i,
                &tmp_i,
                &tmp_d, &tmp_d, &tmp_c, tmp_i,
                &tmp_i, tmp_i, tmp_i, &tmp_i, &tmp_i,
                &tmp_d, &tmp_i);
        
        
        free(yb);
        free(xb);
        free(ytr);
        free(xtmp);
        free(resOOB);
        free(in);
        free(nodex);
        free(varUsed);
        if (!(*replace)  )
            free(nind);
        
        if (testdat) {
            free(ytree);
            free(nodexts);
        }
        
        if (doProx && oobprox)
            free(oobpair) ;
    }
    
    /*----------------------------------------------------------------------*/
    void regForest(double *x, double *ypred, int *mdim, int *n,
                   int *ntree, int *lDaughter, int *rDaughter,
                   small_int *nodestatus, int *nrnodes, double *xsplit,
                   double *avnodes, int *mbest, int *treeSize, int *cat,
                   int maxcat, int *keepPred, double *allpred, int doProx,
                   double *proxMat, int *nodes, int *nodex) {
        int i, j, idx1, idx2, *junk;
        double *ytree;
        
        junk = NULL;
        ytree = (double *) calloc(*n, sizeof(double));
        if (*nodes) {
            zeroInt(nodex, *n * *ntree);
        } else {
            zeroInt(nodex, *n);
        }
        if (doProx) zeroDouble(proxMat, *n * *n);
        if (*keepPred) zeroDouble(allpred, *n * *ntree);
        idx1 = 0;
        idx2 = 0;
        for (i = 0; i < *ntree; ++i) {
            zeroDouble(ytree, *n);
            predictRegTree(x, *n, *mdim, lDaughter + idx1, rDaughter + idx1,
                           nodestatus + idx1, ytree, xsplit + idx1,
                           avnodes + idx1, mbest + idx1, treeSize[i], cat, maxcat,
                           nodex + idx2);
            
            for (j = 0; j < *n; ++j) ypred[j] += ytree[j];
            if (*keepPred) {
                for (j = 0; j < *n; ++j) allpred[j + i * *n] = ytree[j];
            }
            /* if desired, do proximities for this round */
            if (doProx) computeProximity(proxMat, 0, nodex + idx2, junk,
                                         junk, *n);
            idx1 += *nrnodes; /* increment the offset */
            if (*nodes) idx2 += *n;
        }
        for (i = 0; i < *n; ++i) ypred[i] /= *ntree;
        if (doProx) {
            for (i = 0; i < *n; ++i) {
                for (j = i + 1; j < *n; ++j) {
                    proxMat[i + j * *n] /= *ntree;
                    proxMat[j + i * *n] = proxMat[i + j * *n];
                }
                proxMat[i + i * *n] = 1.0;
            }
        }
        free(ytree);
    }
    
    void simpleLinReg(int nsample, double *x, double *y, double *coef,
                      double *mse, int *hasPred) {
        /* Compute simple linear regression of y on x, returning the coefficients,
         the average squared residual, and the predicted values (overwriting y). */
        int i, nout = 0;
        double sxx=0.0, sxy=0.0, xbar=0.0, ybar=0.0;
        double dx = 0.0, dy = 0.0, py=0.0;
        
        for (i = 0; i < nsample; ++i) {
            if (hasPred[i]) {
                nout++;
                xbar += x[i];
                ybar += y[i];
            }
        }
        xbar /= nout;
        ybar /= nout;
        
        for (i = 0; i < nsample; ++i) {
            if (hasPred[i]) {
                dx = x[i] - xbar;
                dy = y[i] - ybar;
                sxx += dx * dx;
                sxy += dx * dy;
            }
        }
        coef[1] = sxy / sxx;
        coef[0] = ybar - coef[1] * xbar;
        
        *mse = 0.0;
        for (i = 0; i < nsample; ++i) {
            if (hasPred[i]) {
                py = coef[0] + coef[1] * x[i];
                dy = y[i] - py;
                *mse += dy * dy;
                /* y[i] = py; */
            }
        }
        *mse /= nout;
        return;
    }
    
    
    void regTree(double *x, double *y, int mdim, int nsample, int *lDaughter,
                 int *rDaughter,
                 double *upper, double *avnode, small_int *nodestatus, int nrnodes,
                 int *treeSize, int nthsize, int mtry, int *mbest, int *cat,
                 double *tgini, int *varUsed) {
        int i, j, k, m, ncur;
        static int *jdex, *nodestart, *nodepop;
        int ndstart, ndend, ndendl, nodecnt, jstat, msplit;
        double d, ss, av, decsplit, ubest, sumnode;
        
        if (in_regTree==-99){
            free(nodestart);
            free(jdex);
            free(nodepop);
            //      Printf("giving up mem in in_regTree\n");
            return;
        }
        
        if (in_regTree==0){
            in_regTree=1;
            nodestart = (int *) calloc(nrnodes, sizeof(int));
            nodepop   = (int *) calloc(nrnodes, sizeof(int));
            jdex = (int *) calloc(nsample, sizeof(int));
        }
        
        /* initialize some arrays for the tree */
        zeroSMALLInt(nodestatus, nrnodes);
        zeroInt(nodestart, nrnodes);
        zeroInt(nodepop, nrnodes);
        zeroDouble(avnode, nrnodes);
        
        for (i = 1; i <= nsample; ++i) jdex[i-1] = i;
        
        ncur = 0;
        nodestart[0] = 0;
        nodepop[0] = nsample;
        nodestatus[0] = NODE_TOSPLIT;
        
        /* compute mean and sum of squares for Y */
        av = 0.0;
        ss = 0.0;
        for (i = 0; i < nsample; ++i) {
            d = y[jdex[i] - 1];
            ss += i * (av - d) * (av - d) / (i + 1);
            av = (i * av + d) / (i + 1);
        }
        avnode[0] = av;
        
        /* start main loop */
        for (k = 0; k < nrnodes - 2; ++k) {
            if (k > ncur || ncur >= nrnodes - 2) break;
            /* skip if the node is not to be split */
            if (nodestatus[k] != NODE_TOSPLIT) continue;
            
            /* initialize for next call to findbestsplit */
            ndstart = nodestart[k];
            ndend = ndstart + nodepop[k] - 1;
            nodecnt = nodepop[k];
            sumnode = nodecnt * avnode[k];
            jstat = 0;
            decsplit = 0.0;
            
            findBestSplit(x, jdex, y, mdim, nsample, ndstart, ndend, &msplit,
                          &decsplit, &ubest, &ndendl, &jstat, mtry, sumnode,
                          nodecnt, cat);
            if (jstat == 1) {
                /* Node is terminal: Mark it as such and move on to the next. */
                nodestatus[k] = NODE_TERMINAL;
                continue;
            }
            /* Found the best split. */
            mbest[k] = msplit;
            varUsed[msplit - 1] = 1;
            upper[k] = ubest;
            tgini[msplit - 1] += decsplit;
            nodestatus[k] = NODE_INTERIOR;
            
            /* leftnode no.= ncur+1, rightnode no. = ncur+2. */
            nodepop[ncur + 1] = ndendl - ndstart + 1;
            nodepop[ncur + 2] = ndend - ndendl;
            nodestart[ncur + 1] = ndstart;
            nodestart[ncur + 2] = ndendl + 1;
            
            /* compute mean and sum of squares for the left daughter node */
            av = 0.0;
            ss = 0.0;
            for (j = ndstart; j <= ndendl; ++j) {
                d = y[jdex[j]-1];
                m = j - ndstart;
                ss += m * (av - d) * (av - d) / (m + 1);
                av = (m * av + d) / (m+1);
            }
            avnode[ncur+1] = av;
            nodestatus[ncur+1] = NODE_TOSPLIT;
            if (nodepop[ncur + 1] <= nthsize) {
                nodestatus[ncur + 1] = NODE_TERMINAL;
            }
            
            /* compute mean and sum of squares for the right daughter node */
            av = 0.0;
            ss = 0.0;
            for (j = ndendl + 1; j <= ndend; ++j) {
                d = y[jdex[j]-1];
                m = j - (ndendl + 1);
                ss += m * (av - d) * (av - d) / (m + 1);
                av = (m * av + d) / (m + 1);
            }
            avnode[ncur + 2] = av;
            nodestatus[ncur + 2] = NODE_TOSPLIT;
            if (nodepop[ncur + 2] <= nthsize) {
                nodestatus[ncur + 2] = NODE_TERMINAL;
            }
            
            /* map the daughter nodes */
            lDaughter[k] = ncur + 1 + 1;
            rDaughter[k] = ncur + 2 + 1;
            /* Augment the tree by two nodes. */
            ncur += 2;
        }
        *treeSize = nrnodes;
        for (k = nrnodes - 1; k >= 0; --k) {
            if (nodestatus[k] == 0) (*treeSize)--;
            if (nodestatus[k] == NODE_TOSPLIT) {
                nodestatus[k] = NODE_TERMINAL;
            }
        }
        
    }
    
    /*--------------------------------------------------------------*/
    
    void findBestSplit(double *x, int *jdex, double *y, int mdim, int nsample,
                       int ndstart, int ndend, int *msplit, double *decsplit,
                       double *ubest, int *ndendl, int *jstat, int mtry,
                       double sumnode, int nodecnt, int *cat) {
        int last, ncat[32], icat[32], lc, nl, nr, npopl, npopr;
        int i, j, kv, l;
        static int *mind, *ncase;
        static double *xt, *ut, *v, *yl;
        double sumcat[32], avcat[32], tavcat[32], ubestt;
        double crit, critmax, critvar, suml, sumr, d, critParent;
        
        
        if (in_findBestSplit==-99){
            free(ncase);
            free(mind); //had to remove this so that it wont crash for when mdim=0, strangely happened for replace=0
            free(v);
            free(yl);
            free(xt);
            free(ut);
            // Printf("giving up mem in findBestSplit\n");
            return;
        }
        
        if (in_findBestSplit==0){
            in_findBestSplit=1;
            ut = (double *) calloc(nsample, sizeof(double));
            xt = (double *) calloc(nsample, sizeof(double));
            v  = (double *) calloc(nsample, sizeof(double));
            yl = (double *) calloc(nsample, sizeof(double));
            mind  = (int *) calloc(mdim+1, sizeof(int));   //seems that the sometimes i am asking for kv[10] and that causes problesmms
            //so allocate 1 more. helps with not crashing in windows
            ncase = (int *) calloc(nsample, sizeof(int));
        }
        zeroDouble(ut, nsample);
        zeroDouble(xt, nsample);
        zeroDouble(v, nsample);
        zeroDouble(yl, nsample);
        zeroInt(mind, mdim);
        zeroInt(ncase, nsample);
        
        zeroDouble(avcat, 32);
        zeroDouble(tavcat, 32);
        
        /* START BIG LOOP */
        *msplit = -1;
        *decsplit = 0.0;
        critmax = 0.0;
        ubestt = 0.0;
        for (i=0; i < mdim; ++i) mind[i] = i;
        
        last = mdim - 1;
        for (i = 0; i < mtry; ++i) {
            critvar = 0.0;
            j = (int) (rnd.unif_rand() * (last+1));
            //Printf("j=%d, last=%d mind[j]=%d\n", j, last, mind[j]);fflush(stdout);
            kv = mind[j];
            //if(kv>100){
            //      1;
            //      getchar();
            //}
            swapInt(mind[j], mind[last]);
            /* mind[j] = mind[last];
             * mind[last] = kv; */
            last--;
            
            lc = cat[kv];
            if (lc == 1) {
                /* numeric variable */
                for (j = ndstart; j <= ndend; ++j) {
                    xt[j] = x[kv + (jdex[j] - 1) * mdim];
                    yl[j] = y[jdex[j] - 1];
                }
            } else {
                /* categorical variable */
                zeroInt(ncat, 32);
                zeroDouble(sumcat, 32);
                for (j = ndstart; j <= ndend; ++j) {
                    l = (int) x[kv + (jdex[j] - 1) * mdim];
                    sumcat[l - 1] += y[jdex[j] - 1];
                    ncat[l - 1] ++;
                }
                /* Compute means of Y by category. */
                for (j = 0; j < lc; ++j) {
                    avcat[j] = ncat[j] ? sumcat[j] / ncat[j] : 0.0;
                }
                /* Make the category mean the `pseudo' X data. */
                for (j = 0; j < nsample; ++j) {
                    xt[j] = avcat[(int) x[kv + (jdex[j] - 1) * mdim] - 1];
                    yl[j] = y[jdex[j] - 1];
                }
            }
            /* copy the x data in this node. */
            for (j = ndstart; j <= ndend; ++j) v[j] = xt[j];
            for (j = 1; j <= nsample; ++j) ncase[j - 1] = j;
            R_qsort_I(v, ncase, ndstart + 1, ndend + 1);
            if (v[ndstart] >= v[ndend]) continue;
            /* ncase(n)=case number of v nth from bottom */
            /* Start from the right and search to the left. */
            critParent = sumnode * sumnode / nodecnt;
            suml = 0.0;
            sumr = sumnode;
            npopl = 0;
            npopr = nodecnt;
            crit = 0.0;
            /* Search through the "gaps" in the x-variable. */
            for (j = ndstart; j <= ndend - 1; ++j) {
                d = yl[ncase[j] - 1];
                suml += d;
                sumr -= d;
                npopl++;
                npopr--;
                if (v[j] < v[j+1]) {
                    crit = (suml * suml / npopl) + (sumr * sumr / npopr) -
                    critParent;
                    if (crit > critvar) {
                        ubestt = (v[j] + v[j+1]) / 2.0;
                        critvar = crit;
                    }
                }
            }
            if (critvar > critmax) {
                *ubest = ubestt;
                *msplit = kv + 1;
                critmax = critvar;
                for (j = ndstart; j <= ndend; ++j) {
                    ut[j] = xt[j];
                }
                if (cat[kv] > 1) {
                    for (j = 0; j < cat[kv]; ++j) tavcat[j] = avcat[j];
                }
            }
        }
        *decsplit = critmax;
        
        /* If best split can not be found, set to terminal node and return. */
        if (*msplit != -1) {
            nl = ndstart;
            for (j = ndstart; j <= ndend; ++j) {
                if (ut[j] <= *ubest) {
                    nl++;
                    ncase[nl-1] = jdex[j];
                }
            }
            *ndendl = imax2(nl - 1, ndstart);
            nr = *ndendl + 1;
            for (j = ndstart; j <= ndend; ++j) {
                if (ut[j] > *ubest) {
                    if (nr >= nsample) break;
                    nr++;
                    ncase[nr - 1] = jdex[j];
                }
            }
            if (*ndendl >= ndend) *ndendl = ndend - 1;
            for (j = ndstart; j <= ndend; ++j) jdex[j] = ncase[j];
            
            lc = cat[*msplit - 1];
            if (lc > 1) {
                for (j = 0; j < lc; ++j) {
                    icat[j] = (tavcat[j] < *ubest) ? 1 : 0;
                }
                *ubest = pack(lc, icat);
            }
        } else *jstat = 1;
        
    }
    /*====================================================================*/
    void predictRegTree(double *x, int nsample, int mdim,
                        int *lDaughter, int *rDaughter, small_int *nodestatus,
                        double *ypred, double *split, double *nodepred,
                        int *splitVar, int treeSize, int *cat, int maxcat,
                        int *nodex) {
        int i, j, k, m, *cbestsplit = NULL;
        unsigned int npack;
        
        /* decode the categorical splits */
        if (maxcat > 1) {
            cbestsplit = (int *) calloc(maxcat * treeSize, sizeof(int));
            zeroInt(cbestsplit, maxcat * treeSize);
            for (i = 0; i < treeSize; ++i) {
                if (nodestatus[i] != NODE_TERMINAL && cat[splitVar[i] - 1] > 1) {
                    npack = (unsigned int) split[i];
                    /* unpack `npack' into bits */
                    for (j = 0; npack; npack >>= 1, ++j) {
                        cbestsplit[j + i*maxcat] = npack & 1;
                    }
                }
            }
        }
        
        for (i = 0; i < nsample; ++i) {
            k = 0;
            while (nodestatus[k] != NODE_TERMINAL) { /* go down the tree */
                m = splitVar[k] - 1;
                if (cat[m] == 1) {
                    k = (x[m + i*mdim] <= split[k]) ?
                    lDaughter[k] - 1 : rDaughter[k] - 1;
                } else if (cbestsplit){
                    /* Split by a categorical predictor */
                    k = cbestsplit[(int) x[m + i * mdim] - 1 + k * maxcat] ?
                    lDaughter[k] - 1 : rDaughter[k] - 1;
                }
            }
            /* terminal node: assign prediction and move on to next */
            ypred[i] = nodepred[k];
            nodex[i] = k + 1;
        }
        if (maxcat > 1) free(cbestsplit);
    }
    
    void zeroSMALLInt(void *x, int length) {
        memset(x, 0, length * sizeof(small_int));
    }
    void zeroInt(int *x, int length) {
        memset(x, 0, length * sizeof(int));
    }
    
    void zeroDouble(double *x, int length) {
        memset(x, 0, length * sizeof(double));
    }
    
    int imax2(int x, int y) {
        return (x < y) ? y : x;
    }
    
    
    int pack(int nBits, int *bits) {
        int i = nBits, pack = 0;
        while (--i >= 0) pack += bits[i] << i;
        return(pack);
    }
    
    void unpack(unsigned int pack, int *bits) {
        /* pack is a 4-byte integer.  The sub. returns icat, an integer array of
         zeroes and ones corresponding to the coefficients in the binary expansion
         of pack. */
        int i;
        for (i = 0; pack != 0; pack >>= 1, ++i) bits[i] = pack & 1;
    }
    
    /* Compute proximity. */
    void computeProximity(double *prox, int oobprox, int *node, int *inbag,
                          int *oobpair, int n) {
        /* Accumulate the number of times a pair of points fall in the same node.
         prox:    n x n proximity matrix
         oobprox: should the accumulation only count OOB cases? (0=no, 1=yes)
         node:    vector of terminal node labels
         inbag:   indicator of whether a case is in-bag
         oobpair: matrix to accumulate the number of times a pair is OOB together
         n:       total number of cases
         */
        int i, j;
        for (i = 0; i < n; ++i) {
            for (j = i+1; j < n; ++j) {
                if (oobprox) {
                    if ((inbag[i] > 0) ^ (inbag[j] > 0)) {
                        oobpair[j*n + i] ++;
                        oobpair[i*n + j] ++;
                        if (node[i] == node[j]) {
                            prox[j*n + i] += 1.0;
                            prox[i*n + j] += 1.0;
                        }
                    }
                } else {
                    if (node[i] == node[j]) {
                        prox[j*n + i] += 1.0;
                        prox[i*n + j] += 1.0;
                    }
                }
            }
        }
    }
    
    void permuteOOB(int m, double *x, int *in, int nsample, int mdim) {
        /* Permute the OOB part of a variable in x.
         * Argument:
         *   m: the variable to be permuted
         *   x: the data matrix (variables in rows)
         *   in: vector indicating which case is OOB
         *   nsample: number of cases in the data
         *   mdim: number of variables in the data
         */
        double *tp, tmp;
        int i, last, k, nOOB = 0;
        
        tp = (double *) calloc(nsample , sizeof(double));
        
        for (i = 0; i < nsample; ++i) {
            /* make a copy of the OOB part of the data into tp (for permuting) */
            if (in[i] == 0) {
                tp[nOOB] = x[m + i*mdim];
                nOOB++;
            }
        }
        /* Permute tp */
        last = nOOB;
        for (i = 0; i < nOOB; ++i) {
            k = (int) (last * rnd.unif_rand());
            tmp = tp[last - 1];
            tp[last - 1] = tp[k];
            tp[k] = tmp;
            last--;
        }
        
        /* Copy the permuted OOB data back into x. */
        nOOB = 0;
        for (i = 0; i < nsample; ++i) {
            if (in[i] == 0) {
                x[m + i*mdim] = tp[nOOB];
                nOOB++;
            }
        }
        free(tp);
    }
    
    void print_regRF_params( int *xdim, int *sampsize,
                            int *nthsize, int *nrnodes, int *nTree, int *mtry, int *imp,
                            int *cat, int maxcat, int *jprint, int doProx, int oobprox,
                            int biasCorr, double *yptr, double *errimp, double *impmat,
                            double *impSD, double *prox, int *treeSize, small_int *nodestatus,
                            int *lDaughter, int *rDaughter, double *avnode, int *mbest,
                            double *upper, double *mse, int *keepf, int *replace,
                            int testdat, double *xts, int *nts, double *yts, int labelts,
                            double *yTestPred, double *proxts, double *msets, double *coef,
                            int *nout, int *inbag)  {
        Printf("n_size %d p_size %d\n", xdim[0], xdim[1]);
        Printf("sampsize %d, nodesize %d nrnodes %d\n", *sampsize, *nthsize, *nrnodes);
        Printf("ntree %d, mtry %d, impor %d, localimp %d, nPerm %d\n", *nTree, *mtry, imp[0], imp[1], imp[2]);
        Printf("maxcat %d, jprint %d, doProx %d, oobProx %d, biasCorr %d\n", maxcat, *jprint, doProx, oobprox, biasCorr);
        Printf("prox %f, keep.forest %d, keep.inbag %d\n", *prox, keepf[0], keepf[1]);
        Printf("replace %d, labelts %d, proxts %f\n", *replace, labelts, *proxts);
    }
};
//=============================================================================================
//
//-----------------------------Start SVM Header------------------------------------------------
//
struct svm_node
{
    int index;
    double value;
};

struct svm_problem
{
    int l;
    double *y;
    struct svm_node **x;
};

enum { C_SVC, NU_SVC, ONE_CLASS, EPSILON_SVR, NU_SVR };	/* svm_type */
enum { LINEAR, POLY, RBF, SIGMOID, PRECOMPUTED }; /* kernel_type */

struct svm_parameter
{
    int svm_type;
    int kernel_type;
    int degree;	/* for poly */
    double gamma;	/* for poly/rbf/sigmoid */
    double coef0;	/* for poly/sigmoid */
    
    /* these are for training only */
    double cache_size; /* in MB */
    double eps;	/* stopping criteria */
    double C;	/* for C_SVC, EPSILON_SVR and NU_SVR */
    int nr_weight;		/* for C_SVC */
    int *weight_label;	/* for C_SVC */
    double* weight;		/* for C_SVC */
    double nu;	/* for NU_SVC, ONE_CLASS, and NU_SVR */
    double p;	/* for EPSILON_SVR */
    int shrinking;	/* use the shrinking heuristics */
    int probability; /* do probability estimates */
};

//
// svm_model
//
struct svm_model
{
    struct svm_parameter param;	/* parameter */
    int nr_class;		/* number of classes, = 2 in regression/one class svm */
    int l;			/* total #SV */
    struct svm_node **SV;		/* SVs (SV[l]) */
    double **sv_coef;	/* coefficients for SVs in decision functions (sv_coef[k-1][l]) */
    double *rho;		/* constants in decision functions (rho[k*(k-1)/2]) */
    double *probA;		/* pariwise probability information */
    double *probB;
    int *sv_indices;        /* sv_indices[0,...,nSV-1] are values in [1,...,num_traning_data] to indicate SVs in the training set */
    
    /* for classification only */
    
    int *label;		/* label of each class (label[k]) */
    int *nSV;		/* number of SVs for each class (nSV[k]) */
				/* nSV[0] + nSV[1] + ... + nSV[k-1] = l */
    /* XXX */
    int free_sv;		/* 1 if svm_model is created by svm_load_model*/
				/* 0 if svm_model is created by svm_train */
};
struct svm_model *svm_train(const struct svm_problem *prob, const struct svm_parameter *param);
void svm_cross_validation(const struct svm_problem *prob, const struct svm_parameter *param, int nr_fold, double *target);

int svm_save_model(const char *model_file_name, const struct svm_model *model);
struct svm_model *svm_load_model(const char *model_file_name);

int svm_get_svm_type(const struct svm_model *model);
int svm_get_nr_class(const struct svm_model *model);
void svm_get_labels(const struct svm_model *model, int *label);
void svm_get_sv_indices(const struct svm_model *model, int *sv_indices);
int svm_get_nr_sv(const struct svm_model *model);
double svm_get_svr_probability(const struct svm_model *model);

double svm_predict_values(const struct svm_model *model, const struct svm_node *x, double* dec_values);
double svm_predict(const struct svm_model *model, const struct svm_node *x);
double svm_predict_probability(const struct svm_model *model, const struct svm_node *x, double* prob_estimates);

void svm_free_model_content(struct svm_model *model_ptr);
void svm_free_and_destroy_model(struct svm_model **model_ptr_ptr);
void svm_destroy_param(struct svm_parameter *param);

const char *svm_check_parameter(const struct svm_problem *prob, const struct svm_parameter *param);
int svm_check_probability_model(const struct svm_model *model);

void svm_set_print_string_function(void (*print_func)(const char *));
//
//-----------------------------Start SVM Implementation----------------------------------------
//

typedef float Qfloat;
typedef signed char schar;

template <class S, class T> static inline void clone(T*& dst, S* src, int n)
{
    dst = new T[n];
    memcpy((void *)dst,(void *)src,sizeof(T)*n);
}
static inline double powi(double base, int times)
{
    double tmp = base, ret = 1.0;
    
    for(int t=times; t>0; t/=2)
    {
        if(t%2==1) ret*=tmp;
        tmp = tmp * tmp;
    }
    return ret;
}
#define INF HUGE_VAL
#define TAU 1e-12
#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

static void print_string_stdout(const char *s)
{
    fputs(s,stdout);
    fflush(stdout);
}
static void (*svm_print_string) (const char *) = &print_string_stdout;
#if 1
static void info(const char *fmt,...)
{
    char buf[BUFSIZ];
    va_list ap;
    va_start(ap,fmt);
    vsprintf(buf,fmt,ap);
    va_end(ap);
    (*svm_print_string)(buf);
}
#else
static void info(const char *fmt,...) {}
#endif

//
// Kernel Cache
//
// l is the number of total data items
// size is the cache size limit in bytes
//
class Cache
{
public:
    Cache(int l,long int size);
    ~Cache();
    
    // request data [0,len)
    // return some position p where [p,len) need to be filled
    // (p >= len if nothing needs to be filled)
    int get_data(const int index, Qfloat **data, int len);
    void swap_index(int i, int j);
private:
    int l;
    long int size;
    struct head_t
    {
        head_t *prev, *next;	// a circular list
        Qfloat *data;
        int len;		// data[0,len) is cached in this entry
    };
    
    head_t *head;
    head_t lru_head;
    void lru_delete(head_t *h);
    void lru_insert(head_t *h);
};

Cache::Cache(int l_,long int size_):l(l_),size(size_)
{
    head = (head_t *)calloc(l,sizeof(head_t));	// initialized to 0
    size /= sizeof(Qfloat);
    size -= l * sizeof(head_t) / sizeof(Qfloat);
    size = max(size, 2 * (long int) l);	// cache must be large enough for two columns
    lru_head.next = lru_head.prev = &lru_head;
}

Cache::~Cache()
{
    for(head_t *h = lru_head.next; h != &lru_head; h=h->next)
        free(h->data);
    free(head);
}

void Cache::lru_delete(head_t *h)
{
    // delete from current location
    h->prev->next = h->next;
    h->next->prev = h->prev;
}

void Cache::lru_insert(head_t *h)
{
    // insert to last position
    h->next = &lru_head;
    h->prev = lru_head.prev;
    h->prev->next = h;
    h->next->prev = h;
}

int Cache::get_data(const int index, Qfloat **data, int len)
{
    head_t *h = &head[index];
    if(h->len) lru_delete(h);
    int more = len - h->len;
    
    if(more > 0)
    {
        // free old space
        while(size < more)
        {
            head_t *old = lru_head.next;
            lru_delete(old);
            free(old->data);
            size += old->len;
            old->data = 0;
            old->len = 0;
        }
        
        // allocate new space
        h->data = (Qfloat *)realloc(h->data,sizeof(Qfloat)*len);
        size -= more;
        swap(h->len,len);
    }
    
    lru_insert(h);
    *data = h->data;
    return len;
}

void Cache::swap_index(int i, int j)
{
    if(i==j) return;
    
    if(head[i].len) lru_delete(&head[i]);
    if(head[j].len) lru_delete(&head[j]);
    swap(head[i].data,head[j].data);
    swap(head[i].len,head[j].len);
    if(head[i].len) lru_insert(&head[i]);
    if(head[j].len) lru_insert(&head[j]);
    
    if(i>j) swap(i,j);
    for(head_t *h = lru_head.next; h!=&lru_head; h=h->next)
    {
        if(h->len > i)
        {
            if(h->len > j)
                swap(h->data[i],h->data[j]);
            else
            {
                // give up
                lru_delete(h);
                free(h->data);
                size += h->len;
                h->data = 0;
                h->len = 0;
            }
        }
    }
}

//
// Kernel evaluation
//
// the static method k_function is for doing single kernel evaluation
// the constructor of Kernel prepares to calculate the l*l kernel matrix
// the member function get_Q is for getting one column from the Q Matrix
//
class QMatrix {
public:
    virtual Qfloat *get_Q(int column, int len) const = 0;
    virtual double *get_QD() const = 0;
    virtual void swap_index(int i, int j) const = 0;
    virtual ~QMatrix() {}
};

class Kernel: public QMatrix {
public:
    Kernel(int l, svm_node * const * x, const svm_parameter& param);
    virtual ~Kernel();
    
    static double k_function(const svm_node *x, const svm_node *y,
                             const svm_parameter& param);
    virtual Qfloat *get_Q(int column, int len) const = 0;
    virtual double *get_QD() const = 0;
    virtual void swap_index(int i, int j) const	// no so const...
    {
        swap(x[i],x[j]);
        if(x_square) swap(x_square[i],x_square[j]);
    }
protected:
    
    double (Kernel::*kernel_function)(int i, int j) const;
    
private:
    const svm_node **x;
    double *x_square;
    
    // svm_parameter
    const int kernel_type;
    const int degree;
    const double gamma;
    const double coef0;
    
    static double dot(const svm_node *px, const svm_node *py);
    double kernel_linear(int i, int j) const
    {
        return dot(x[i],x[j]);
    }
    double kernel_poly(int i, int j) const
    {
        return powi(gamma*dot(x[i],x[j])+coef0,degree);
    }
    double kernel_rbf(int i, int j) const
    {
        return exp(-gamma*(x_square[i]+x_square[j]-2*dot(x[i],x[j])));
    }
    double kernel_sigmoid(int i, int j) const
    {
        return tanh(gamma*dot(x[i],x[j])+coef0);
    }
    double kernel_precomputed(int i, int j) const
    {
        return x[i][(int)(x[j][0].value)].value;
    }
};

Kernel::Kernel(int l, svm_node * const * x_, const svm_parameter& param)
:kernel_type(param.kernel_type), degree(param.degree),
gamma(param.gamma), coef0(param.coef0)
{
    switch(kernel_type)
    {
        case LINEAR:
            kernel_function = &Kernel::kernel_linear;
            break;
        case POLY:
            kernel_function = &Kernel::kernel_poly;
            break;
        case RBF:
            kernel_function = &Kernel::kernel_rbf;
            break;
        case SIGMOID:
            kernel_function = &Kernel::kernel_sigmoid;
            break;
        case PRECOMPUTED:
            kernel_function = &Kernel::kernel_precomputed;
            break;
    }
    
    clone(x,x_,l);
    
    if(kernel_type == RBF)
    {
        x_square = new double[l];
        for(int i=0;i<l;i++)
            x_square[i] = dot(x[i],x[i]);
    }
    else
        x_square = 0;
}

Kernel::~Kernel()
{
    delete[] x;
    delete[] x_square;
}

double Kernel::dot(const svm_node *px, const svm_node *py)
{
    double sum = 0;
    while(px->index != -1 && py->index != -1)
    {
        if(px->index == py->index)
        {
            sum += px->value * py->value;
            ++px;
            ++py;
        }
        else
        {
            if(px->index > py->index)
                ++py;
            else
                ++px;
        }
    }
    return sum;
}

double Kernel::k_function(const svm_node *x, const svm_node *y,
                          const svm_parameter& param)
{
    switch(param.kernel_type)
    {
        case LINEAR:
            return dot(x,y);
        case POLY:
            return powi(param.gamma*dot(x,y)+param.coef0,param.degree);
        case RBF:
        {
            double sum = 0;
            while(x->index != -1 && y->index !=-1)
            {
                if(x->index == y->index)
                {
                    double d = x->value - y->value;
                    sum += d*d;
                    ++x;
                    ++y;
                }
                else
                {
                    if(x->index > y->index)
                    {
                        sum += y->value * y->value;
                        ++y;
                    }
                    else
                    {
                        sum += x->value * x->value;
                        ++x;
                    }
                }
            }
            
            while(x->index != -1)
            {
                sum += x->value * x->value;
                ++x;
            }
            
            while(y->index != -1)
            {
                sum += y->value * y->value;
                ++y;
            }
            
            return exp(-param.gamma*sum);
        }
        case SIGMOID:
            return tanh(param.gamma*dot(x,y)+param.coef0);
        case PRECOMPUTED:  //x: test (validation), y: SV
            return x[(int)(y->value)].value;
        default:
            return 0;  // Unreachable
    }
}

// An SMO algorithm in Fan et al., JMLR 6(2005), p. 1889--1918
// Solves:
//
//	min 0.5(\alpha^T Q \alpha) + p^T \alpha
//
//		y^T \alpha = \delta
//		y_i = +1 or -1
//		0 <= alpha_i <= Cp for y_i = 1
//		0 <= alpha_i <= Cn for y_i = -1
//
// Given:
//
//	Q, p, y, Cp, Cn, and an initial feasible point \alpha
//	l is the size of vectors and matrices
//	eps is the stopping tolerance
//
// solution will be put in \alpha, objective value will be put in obj
//
class Solver {
public:
    Solver() {};
    virtual ~Solver() {};
    
    struct SolutionInfo {
        double obj;
        double rho;
        double upper_bound_p;
        double upper_bound_n;
        double r;	// for Solver_NU
    };
    
    void Solve(int l, const QMatrix& Q, const double *p_, const schar *y_,
               double *alpha_, double Cp, double Cn, double eps,
               SolutionInfo* si, int shrinking);
protected:
    int active_size;
    schar *y;
    double *G;		// gradient of objective function
    enum { LOWER_BOUND, UPPER_BOUND, FREE };
    char *alpha_status;	// LOWER_BOUND, UPPER_BOUND, FREE
    double *alpha;
    const QMatrix *Q;
    const double *QD;
    double eps;
    double Cp,Cn;
    double *p;
    int *active_set;
    double *G_bar;		// gradient, if we treat free variables as 0
    int l;
    bool unshrink;	// XXX
    
    double get_C(int i)
    {
        return (y[i] > 0)? Cp : Cn;
    }
    void update_alpha_status(int i)
    {
        if(alpha[i] >= get_C(i))
            alpha_status[i] = UPPER_BOUND;
        else if(alpha[i] <= 0)
            alpha_status[i] = LOWER_BOUND;
        else alpha_status[i] = FREE;
    }
    bool is_upper_bound(int i) { return alpha_status[i] == UPPER_BOUND; }
    bool is_lower_bound(int i) { return alpha_status[i] == LOWER_BOUND; }
    bool is_free(int i) { return alpha_status[i] == FREE; }
    void swap_index(int i, int j);
    void reconstruct_gradient();
    virtual int select_working_set(int &i, int &j);
    virtual double calculate_rho();
    virtual void do_shrinking();
private:
    bool be_shrunk(int i, double Gmax1, double Gmax2);
};

void Solver::swap_index(int i, int j)
{
    Q->swap_index(i,j);
    swap(y[i],y[j]);
    swap(G[i],G[j]);
    swap(alpha_status[i],alpha_status[j]);
    swap(alpha[i],alpha[j]);
    swap(p[i],p[j]);
    swap(active_set[i],active_set[j]);
    swap(G_bar[i],G_bar[j]);
}

void Solver::reconstruct_gradient()
{
    // reconstruct inactive elements of G from G_bar and free variables
    
    if(active_size == l) return;
    
    int i,j;
    int nr_free = 0;
    
    for(j=active_size;j<l;j++)
        G[j] = G_bar[j] + p[j];
    
    for(j=0;j<active_size;j++)
        if(is_free(j))
            nr_free++;
    
    if(2*nr_free < active_size)
        info("\nWARNING: using -h 0 may be faster\n");
    
    if (nr_free*l > 2*active_size*(l-active_size))
    {
        for(i=active_size;i<l;i++)
        {
            const Qfloat *Q_i = Q->get_Q(i,active_size);
            for(j=0;j<active_size;j++)
                if(is_free(j))
                    G[i] += alpha[j] * Q_i[j];
        }
    }
    else
    {
        for(i=0;i<active_size;i++)
            if(is_free(i))
            {
                const Qfloat *Q_i = Q->get_Q(i,l);
                double alpha_i = alpha[i];
                for(j=active_size;j<l;j++)
                    G[j] += alpha_i * Q_i[j];
            }
    }
}

void Solver::Solve(int l, const QMatrix& Q, const double *p_, const schar *y_,
                   double *alpha_, double Cp, double Cn, double eps,
                   SolutionInfo* si, int shrinking)
{
    this->l = l;
    this->Q = &Q;
    QD=Q.get_QD();
    clone(p, p_,l);
    clone(y, y_,l);
    clone(alpha,alpha_,l);
    this->Cp = Cp;
    this->Cn = Cn;
    this->eps = eps;
    unshrink = false;
    
    // initialize alpha_status
    {
        alpha_status = new char[l];
        for(int i=0;i<l;i++)
            update_alpha_status(i);
    }
    
    // initialize active set (for shrinking)
    {
        active_set = new int[l];
        for(int i=0;i<l;i++)
            active_set[i] = i;
        active_size = l;
    }
    
    // initialize gradient
    {
        G = new double[l];
        G_bar = new double[l];
        int i;
        for(i=0;i<l;i++)
        {
            G[i] = p[i];
            G_bar[i] = 0;
        }
        for(i=0;i<l;i++)
            if(!is_lower_bound(i))
            {
                const Qfloat *Q_i = Q.get_Q(i,l);
                double alpha_i = alpha[i];
                int j;
                for(j=0;j<l;j++)
                    G[j] += alpha_i*Q_i[j];
                if(is_upper_bound(i))
                    for(j=0;j<l;j++)
                        G_bar[j] += get_C(i) * Q_i[j];
            }
    }
    
    // optimization step
    
    int iter = 0;
    int max_iter = max(10000000, l>INT_MAX/100 ? INT_MAX : 100*l);
    int counter = min(l,1000)+1;
    
    while(iter < max_iter)
    {
        // show progress and do shrinking
        
        if(--counter == 0)
        {
            counter = min(l,1000);
            if(shrinking) do_shrinking();
            info(".");
        }
        
        int i,j;
        if(select_working_set(i,j)!=0)
        {
            // reconstruct the whole gradient
            reconstruct_gradient();
            // reset active set size and check
            active_size = l;
            info("*");
            if(select_working_set(i,j)!=0)
                break;
            else
                counter = 1;	// do shrinking next iteration
        }
        
        ++iter;
        
        // update alpha[i] and alpha[j], handle bounds carefully
        
        const Qfloat *Q_i = Q.get_Q(i,active_size);
        const Qfloat *Q_j = Q.get_Q(j,active_size);
        
        double C_i = get_C(i);
        double C_j = get_C(j);
        
        double old_alpha_i = alpha[i];
        double old_alpha_j = alpha[j];
        
        if(y[i]!=y[j])
        {
            double quad_coef = QD[i]+QD[j]+2*Q_i[j];
            if (quad_coef <= 0)
                quad_coef = TAU;
            double delta = (-G[i]-G[j])/quad_coef;
            double diff = alpha[i] - alpha[j];
            alpha[i] += delta;
            alpha[j] += delta;
            
            if(diff > 0)
            {
                if(alpha[j] < 0)
                {
                    alpha[j] = 0;
                    alpha[i] = diff;
                }
            }
            else
            {
                if(alpha[i] < 0)
                {
                    alpha[i] = 0;
                    alpha[j] = -diff;
                }
            }
            if(diff > C_i - C_j)
            {
                if(alpha[i] > C_i)
                {
                    alpha[i] = C_i;
                    alpha[j] = C_i - diff;
                }
            }
            else
            {
                if(alpha[j] > C_j)
                {
                    alpha[j] = C_j;
                    alpha[i] = C_j + diff;
                }
            }
        }
        else
        {
            double quad_coef = QD[i]+QD[j]-2*Q_i[j];
            if (quad_coef <= 0)
                quad_coef = TAU;
            double delta = (G[i]-G[j])/quad_coef;
            double sum = alpha[i] + alpha[j];
            alpha[i] -= delta;
            alpha[j] += delta;
            
            if(sum > C_i)
            {
                if(alpha[i] > C_i)
                {
                    alpha[i] = C_i;
                    alpha[j] = sum - C_i;
                }
            }
            else
            {
                if(alpha[j] < 0)
                {
                    alpha[j] = 0;
                    alpha[i] = sum;
                }
            }
            if(sum > C_j)
            {
                if(alpha[j] > C_j)
                {
                    alpha[j] = C_j;
                    alpha[i] = sum - C_j;
                }
            }
            else
            {
                if(alpha[i] < 0)
                {
                    alpha[i] = 0;
                    alpha[j] = sum;
                }
            }
        }
        
        // update G
        
        double delta_alpha_i = alpha[i] - old_alpha_i;
        double delta_alpha_j = alpha[j] - old_alpha_j;
        
        for(int k=0;k<active_size;k++)
        {
            G[k] += Q_i[k]*delta_alpha_i + Q_j[k]*delta_alpha_j;
        }
        
        // update alpha_status and G_bar
        
        {
            bool ui = is_upper_bound(i);
            bool uj = is_upper_bound(j);
            update_alpha_status(i);
            update_alpha_status(j);
            int k;
            if(ui != is_upper_bound(i))
            {
                Q_i = Q.get_Q(i,l);
                if(ui)
                    for(k=0;k<l;k++)
                        G_bar[k] -= C_i * Q_i[k];
                else
                    for(k=0;k<l;k++)
                        G_bar[k] += C_i * Q_i[k];
            }
            
            if(uj != is_upper_bound(j))
            {
                Q_j = Q.get_Q(j,l);
                if(uj)
                    for(k=0;k<l;k++)
                        G_bar[k] -= C_j * Q_j[k];
                else
                    for(k=0;k<l;k++)
                        G_bar[k] += C_j * Q_j[k];
            }
        }
    }
    
    if(iter >= max_iter)
    {
        if(active_size < l)
        {
            // reconstruct the whole gradient to calculate objective value
            reconstruct_gradient();
            active_size = l;
            info("*");
        }
        fprintf(stderr,"\nWARNING: reaching max number of iterations\n");
    }
    
    // calculate rho
    
    si->rho = calculate_rho();
    
    // calculate objective value
    {
        double v = 0;
        int i;
        for(i=0;i<l;i++)
            v += alpha[i] * (G[i] + p[i]);
        
        si->obj = v/2;
    }
    
    // put back the solution
    {
        for(int i=0;i<l;i++)
            alpha_[active_set[i]] = alpha[i];
    }
    
    // juggle everything back
    /*{
     for(int i=0;i<l;i++)
     while(active_set[i] != i)
     swap_index(i,active_set[i]);
     // or Q.swap_index(i,active_set[i]);
     }*/
    
    si->upper_bound_p = Cp;
    si->upper_bound_n = Cn;
    
    info("\noptimization finished, #iter = %d\n",iter);
    
    delete[] p;
    delete[] y;
    delete[] alpha;
    delete[] alpha_status;
    delete[] active_set;
    delete[] G;
    delete[] G_bar;
}

// return 1 if already optimal, return 0 otherwise
int Solver::select_working_set(int &out_i, int &out_j)
{
    // return i,j such that
    // i: maximizes -y_i * grad(f)_i, i in I_up(\alpha)
    // j: minimizes the decrease of obj value
    //    (if quadratic coefficeint <= 0, replace it with tau)
    //    -y_j*grad(f)_j < -y_i*grad(f)_i, j in I_low(\alpha)
    
    double Gmax = -INF;
    double Gmax2 = -INF;
    int Gmax_idx = -1;
    int Gmin_idx = -1;
    double obj_diff_min = INF;
    
    for(int t=0;t<active_size;t++)
        if(y[t]==+1)
        {
            if(!is_upper_bound(t))
                if(-G[t] >= Gmax)
                {
                    Gmax = -G[t];
                    Gmax_idx = t;
                }
        }
        else
        {
            if(!is_lower_bound(t))
                if(G[t] >= Gmax)
                {
                    Gmax = G[t];
                    Gmax_idx = t;
                }
        }
    
    int i = Gmax_idx;
    const Qfloat *Q_i = NULL;
    if(i != -1) // NULL Q_i not accessed: Gmax=-INF if i=-1
        Q_i = Q->get_Q(i,active_size);
    
    for(int j=0;j<active_size;j++)
    {
        if(y[j]==+1)
        {
            if (!is_lower_bound(j))
            {
                double grad_diff=Gmax+G[j];
                if (G[j] >= Gmax2)
                    Gmax2 = G[j];
                if (grad_diff > 0)
                {
                    double obj_diff;
                    double quad_coef = QD[i]+QD[j]-2.0*y[i]*Q_i[j];
                    if (quad_coef > 0)
                        obj_diff = -(grad_diff*grad_diff)/quad_coef;
                    else
                        obj_diff = -(grad_diff*grad_diff)/TAU;
                    
                    if (obj_diff <= obj_diff_min)
                    {
                        Gmin_idx=j;
                        obj_diff_min = obj_diff;
                    }
                }
            }
        }
        else
        {
            if (!is_upper_bound(j))
            {
                double grad_diff= Gmax-G[j];
                if (-G[j] >= Gmax2)
                    Gmax2 = -G[j];
                if (grad_diff > 0)
                {
                    double obj_diff;
                    double quad_coef = QD[i]+QD[j]+2.0*y[i]*Q_i[j];
                    if (quad_coef > 0)
                        obj_diff = -(grad_diff*grad_diff)/quad_coef;
                    else
                        obj_diff = -(grad_diff*grad_diff)/TAU;
                    
                    if (obj_diff <= obj_diff_min)
                    {
                        Gmin_idx=j;
                        obj_diff_min = obj_diff;
                    }
                }
            }
        }
    }
    
    if(Gmax+Gmax2 < eps)
        return 1;
    
    out_i = Gmax_idx;
    out_j = Gmin_idx;
    return 0;
}

bool Solver::be_shrunk(int i, double Gmax1, double Gmax2)
{
    if(is_upper_bound(i))
    {
        if(y[i]==+1)
            return(-G[i] > Gmax1);
        else
            return(-G[i] > Gmax2);
    }
    else if(is_lower_bound(i))
    {
        if(y[i]==+1)
            return(G[i] > Gmax2);
        else
            return(G[i] > Gmax1);
    }
    else
        return(false);
}

void Solver::do_shrinking()
{
    int i;
    double Gmax1 = -INF;		// max { -y_i * grad(f)_i | i in I_up(\alpha) }
    double Gmax2 = -INF;		// max { y_i * grad(f)_i | i in I_low(\alpha) }
    
    // find maximal violating pair first
    for(i=0;i<active_size;i++)
    {
        if(y[i]==+1)
        {
            if(!is_upper_bound(i))
            {
                if(-G[i] >= Gmax1)
                    Gmax1 = -G[i];
            }
            if(!is_lower_bound(i))
            {
                if(G[i] >= Gmax2)
                    Gmax2 = G[i];
            }
        }
        else
        {
            if(!is_upper_bound(i))
            {
                if(-G[i] >= Gmax2)
                    Gmax2 = -G[i];
            }
            if(!is_lower_bound(i))
            {
                if(G[i] >= Gmax1)
                    Gmax1 = G[i];
            }
        }
    }
    
    if(unshrink == false && Gmax1 + Gmax2 <= eps*10)
    {
        unshrink = true;
        reconstruct_gradient();
        active_size = l;
        info("*");
    }
    
    for(i=0;i<active_size;i++)
        if (be_shrunk(i, Gmax1, Gmax2))
        {
            active_size--;
            while (active_size > i)
            {
                if (!be_shrunk(active_size, Gmax1, Gmax2))
                {
                    swap_index(i,active_size);
                    break;
                }
                active_size--;
            }
        }
}

double Solver::calculate_rho()
{
    double r;
    int nr_free = 0;
    double ub = INF, lb = -INF, sum_free = 0;
    for(int i=0;i<active_size;i++)
    {
        double yG = y[i]*G[i];
        
        if(is_upper_bound(i))
        {
            if(y[i]==-1)
                ub = min(ub,yG);
            else
                lb = max(lb,yG);
        }
        else if(is_lower_bound(i))
        {
            if(y[i]==+1)
                ub = min(ub,yG);
            else
                lb = max(lb,yG);
        }
        else
        {
            ++nr_free;
            sum_free += yG;
        }
    }
    
    if(nr_free>0)
        r = sum_free/nr_free;
    else
        r = (ub+lb)/2;
    
    return r;
}

//
// Solver for nu-svm classification and regression
//
// additional constraint: e^T \alpha = constant
//
class Solver_NU: public Solver
{
public:
    Solver_NU() {}
    void Solve(int l, const QMatrix& Q, const double *p, const schar *y,
               double *alpha, double Cp, double Cn, double eps,
               SolutionInfo* si, int shrinking)
    {
        this->si = si;
        Solver::Solve(l,Q,p,y,alpha,Cp,Cn,eps,si,shrinking);
    }
private:
    SolutionInfo *si;
    int select_working_set(int &i, int &j);
    double calculate_rho();
    bool be_shrunk(int i, double Gmax1, double Gmax2, double Gmax3, double Gmax4);
    void do_shrinking();
};

// return 1 if already optimal, return 0 otherwise
int Solver_NU::select_working_set(int &out_i, int &out_j)
{
    // return i,j such that y_i = y_j and
    // i: maximizes -y_i * grad(f)_i, i in I_up(\alpha)
    // j: minimizes the decrease of obj value
    //    (if quadratic coefficeint <= 0, replace it with tau)
    //    -y_j*grad(f)_j < -y_i*grad(f)_i, j in I_low(\alpha)
    
    double Gmaxp = -INF;
    double Gmaxp2 = -INF;
    int Gmaxp_idx = -1;
    
    double Gmaxn = -INF;
    double Gmaxn2 = -INF;
    int Gmaxn_idx = -1;
    
    int Gmin_idx = -1;
    double obj_diff_min = INF;
    
    for(int t=0;t<active_size;t++)
        if(y[t]==+1)
        {
            if(!is_upper_bound(t))
                if(-G[t] >= Gmaxp)
                {
                    Gmaxp = -G[t];
                    Gmaxp_idx = t;
                }
        }
        else
        {
            if(!is_lower_bound(t))
                if(G[t] >= Gmaxn)
                {
                    Gmaxn = G[t];
                    Gmaxn_idx = t;
                }
        }
    
    int ip = Gmaxp_idx;
    int in = Gmaxn_idx;
    const Qfloat *Q_ip = NULL;
    const Qfloat *Q_in = NULL;
    if(ip != -1) // NULL Q_ip not accessed: Gmaxp=-INF if ip=-1
        Q_ip = Q->get_Q(ip,active_size);
    if(in != -1)
        Q_in = Q->get_Q(in,active_size);
    
    for(int j=0;j<active_size;j++)
    {
        if(y[j]==+1)
        {
            if (!is_lower_bound(j))
            {
                double grad_diff=Gmaxp+G[j];
                if (G[j] >= Gmaxp2)
                    Gmaxp2 = G[j];
                if (grad_diff > 0)
                {
                    double obj_diff;
                    double quad_coef = QD[ip]+QD[j]-2*Q_ip[j];
                    if (quad_coef > 0)
                        obj_diff = -(grad_diff*grad_diff)/quad_coef;
                    else
                        obj_diff = -(grad_diff*grad_diff)/TAU;
                    
                    if (obj_diff <= obj_diff_min)
                    {
                        Gmin_idx=j;
                        obj_diff_min = obj_diff;
                    }
                }
            }
        }
        else
        {
            if (!is_upper_bound(j))
            {
                double grad_diff=Gmaxn-G[j];
                if (-G[j] >= Gmaxn2)
                    Gmaxn2 = -G[j];
                if (grad_diff > 0)
                {
                    double obj_diff;
                    double quad_coef = QD[in]+QD[j]-2*Q_in[j];
                    if (quad_coef > 0)
                        obj_diff = -(grad_diff*grad_diff)/quad_coef;
                    else
                        obj_diff = -(grad_diff*grad_diff)/TAU;
                    
                    if (obj_diff <= obj_diff_min)
                    {
                        Gmin_idx=j;
                        obj_diff_min = obj_diff;
                    }
                }
            }
        }
    }
    
    if(max(Gmaxp+Gmaxp2,Gmaxn+Gmaxn2) < eps)
        return 1;
    
    if (y[Gmin_idx] == +1)
        out_i = Gmaxp_idx;
    else
        out_i = Gmaxn_idx;
    out_j = Gmin_idx;
    
    return 0;
}

bool Solver_NU::be_shrunk(int i, double Gmax1, double Gmax2, double Gmax3, double Gmax4)
{
    if(is_upper_bound(i))
    {
        if(y[i]==+1)
            return(-G[i] > Gmax1);
        else
            return(-G[i] > Gmax4);
    }
    else if(is_lower_bound(i))
    {
        if(y[i]==+1)
            return(G[i] > Gmax2);
        else
            return(G[i] > Gmax3);
    }
    else
        return(false);
}

void Solver_NU::do_shrinking()
{
    double Gmax1 = -INF;	// max { -y_i * grad(f)_i | y_i = +1, i in I_up(\alpha) }
    double Gmax2 = -INF;	// max { y_i * grad(f)_i | y_i = +1, i in I_low(\alpha) }
    double Gmax3 = -INF;	// max { -y_i * grad(f)_i | y_i = -1, i in I_up(\alpha) }
    double Gmax4 = -INF;	// max { y_i * grad(f)_i | y_i = -1, i in I_low(\alpha) }
    
    // find maximal violating pair first
    int i;
    for(i=0;i<active_size;i++)
    {
        if(!is_upper_bound(i))
        {
            if(y[i]==+1)
            {
                if(-G[i] > Gmax1) Gmax1 = -G[i];
            }
            else	if(-G[i] > Gmax4) Gmax4 = -G[i];
        }
        if(!is_lower_bound(i))
        {
            if(y[i]==+1)
            {
                if(G[i] > Gmax2) Gmax2 = G[i];
            }
            else	if(G[i] > Gmax3) Gmax3 = G[i];
        }
    }
    
    if(unshrink == false && max(Gmax1+Gmax2,Gmax3+Gmax4) <= eps*10)
    {
        unshrink = true;
        reconstruct_gradient();
        active_size = l;
    }
    
    for(i=0;i<active_size;i++)
        if (be_shrunk(i, Gmax1, Gmax2, Gmax3, Gmax4))
        {
            active_size--;
            while (active_size > i)
            {
                if (!be_shrunk(active_size, Gmax1, Gmax2, Gmax3, Gmax4))
                {
                    swap_index(i,active_size);
                    break;
                }
                active_size--;
            }
        }
}

double Solver_NU::calculate_rho()
{
    int nr_free1 = 0,nr_free2 = 0;
    double ub1 = INF, ub2 = INF;
    double lb1 = -INF, lb2 = -INF;
    double sum_free1 = 0, sum_free2 = 0;
    
    for(int i=0;i<active_size;i++)
    {
        if(y[i]==+1)
        {
            if(is_upper_bound(i))
                lb1 = max(lb1,G[i]);
            else if(is_lower_bound(i))
                ub1 = min(ub1,G[i]);
            else
            {
                ++nr_free1;
                sum_free1 += G[i];
            }
        }
        else
        {
            if(is_upper_bound(i))
                lb2 = max(lb2,G[i]);
            else if(is_lower_bound(i))
                ub2 = min(ub2,G[i]);
            else
            {
                ++nr_free2;
                sum_free2 += G[i];
            }
        }
    }
    
    double r1,r2;
    if(nr_free1 > 0)
        r1 = sum_free1/nr_free1;
    else
        r1 = (ub1+lb1)/2;
    
    if(nr_free2 > 0)
        r2 = sum_free2/nr_free2;
    else
        r2 = (ub2+lb2)/2;
    
    si->r = (r1+r2)/2;
    return (r1-r2)/2;
}

//
// Q matrices for various formulations
//
class SVC_Q: public Kernel
{
public:
    SVC_Q(const svm_problem& prob, const svm_parameter& param, const schar *y_)
    :Kernel(prob.l, prob.x, param)
    {
        clone(y,y_,prob.l);
        cache = new Cache(prob.l,(long int)(param.cache_size*(1<<20)));
        QD = new double[prob.l];
        for(int i=0;i<prob.l;i++)
            QD[i] = (this->*kernel_function)(i,i);
    }
    
    Qfloat *get_Q(int i, int len) const
    {
        Qfloat *data;
        int start, j;
        if((start = cache->get_data(i,&data,len)) < len)
        {
            for(j=start;j<len;j++)
                data[j] = (Qfloat)(y[i]*y[j]*(this->*kernel_function)(i,j));
        }
        return data;
    }
    
    double *get_QD() const
    {
        return QD;
    }
    
    void swap_index(int i, int j) const
    {
        cache->swap_index(i,j);
        Kernel::swap_index(i,j);
        swap(y[i],y[j]);
        swap(QD[i],QD[j]);
    }
    
    ~SVC_Q()
    {
        delete[] y;
        delete cache;
        delete[] QD;
    }
private:
    schar *y;
    Cache *cache;
    double *QD;
};

class ONE_CLASS_Q: public Kernel
{
public:
    ONE_CLASS_Q(const svm_problem& prob, const svm_parameter& param)
    :Kernel(prob.l, prob.x, param)
    {
        cache = new Cache(prob.l,(long int)(param.cache_size*(1<<20)));
        QD = new double[prob.l];
        for(int i=0;i<prob.l;i++)
            QD[i] = (this->*kernel_function)(i,i);
    }
    
    Qfloat *get_Q(int i, int len) const
    {
        Qfloat *data;
        int start, j;
        if((start = cache->get_data(i,&data,len)) < len)
        {
            for(j=start;j<len;j++)
                data[j] = (Qfloat)(this->*kernel_function)(i,j);
        }
        return data;
    }
    
    double *get_QD() const
    {
        return QD;
    }
    
    void swap_index(int i, int j) const
    {
        cache->swap_index(i,j);
        Kernel::swap_index(i,j);
        swap(QD[i],QD[j]);
    }
    
    ~ONE_CLASS_Q()
    {
        delete cache;
        delete[] QD;
    }
private:
    Cache *cache;
    double *QD;
};

class SVR_Q: public Kernel
{
public:
    SVR_Q(const svm_problem& prob, const svm_parameter& param)
    :Kernel(prob.l, prob.x, param)
    {
        l = prob.l;
        cache = new Cache(l,(long int)(param.cache_size*(1<<20)));
        QD = new double[2*l];
        sign = new schar[2*l];
        index = new int[2*l];
        for(int k=0;k<l;k++)
        {
            sign[k] = 1;
            sign[k+l] = -1;
            index[k] = k;
            index[k+l] = k;
            QD[k] = (this->*kernel_function)(k,k);
            QD[k+l] = QD[k];
        }
        buffer[0] = new Qfloat[2*l];
        buffer[1] = new Qfloat[2*l];
        next_buffer = 0;
    }
    
    void swap_index(int i, int j) const
    {
        swap(sign[i],sign[j]);
        swap(index[i],index[j]);
        swap(QD[i],QD[j]);
    }
    
    Qfloat *get_Q(int i, int len) const
    {
        Qfloat *data;
        int j, real_i = index[i];
        if(cache->get_data(real_i,&data,l) < l)
        {
            for(j=0;j<l;j++)
                data[j] = (Qfloat)(this->*kernel_function)(real_i,j);
        }
        
        // reorder and copy
        Qfloat *buf = buffer[next_buffer];
        next_buffer = 1 - next_buffer;
        schar si = sign[i];
        for(j=0;j<len;j++)
            buf[j] = (Qfloat) si * (Qfloat) sign[j] * data[index[j]];
        return buf;
    }
    
    double *get_QD() const
    {
        return QD;
    }
    
    ~SVR_Q()
    {
        delete cache;
        delete[] sign;
        delete[] index;
        delete[] buffer[0];
        delete[] buffer[1];
        delete[] QD;
    }
private:
    int l;
    Cache *cache;
    schar *sign;
    int *index;
    mutable int next_buffer;
    Qfloat *buffer[2];
    double *QD;
};

//
// construct and solve various formulations
//
static void solve_c_svc(
                        const svm_problem *prob, const svm_parameter* param,
                        double *alpha, Solver::SolutionInfo* si, double Cp, double Cn)
{
    int l = prob->l;
    double *minus_ones = new double[l];
    schar *y = new schar[l];
    
    int i;
    
    for(i=0;i<l;i++)
    {
        alpha[i] = 0;
        minus_ones[i] = -1;
        if(prob->y[i] > 0) y[i] = +1; else y[i] = -1;
    }
    
    Solver s;
    s.Solve(l, SVC_Q(*prob,*param,y), minus_ones, y,
            alpha, Cp, Cn, param->eps, si, param->shrinking);
    
    double sum_alpha=0;
    for(i=0;i<l;i++)
        sum_alpha += alpha[i];
    
    if (Cp==Cn)
        info("nu = %f\n", sum_alpha/(Cp*prob->l));
    
    for(i=0;i<l;i++)
        alpha[i] *= y[i];
    
    delete[] minus_ones;
    delete[] y;
}

static void solve_nu_svc(
                         const svm_problem *prob, const svm_parameter *param,
                         double *alpha, Solver::SolutionInfo* si)
{
    int i;
    int l = prob->l;
    double nu = param->nu;
    
    schar *y = new schar[l];
    
    for(i=0;i<l;i++)
        if(prob->y[i]>0)
            y[i] = +1;
        else
            y[i] = -1;
    
    double sum_pos = nu*l/2;
    double sum_neg = nu*l/2;
    
    for(i=0;i<l;i++)
        if(y[i] == +1)
        {
            alpha[i] = min(1.0,sum_pos);
            sum_pos -= alpha[i];
        }
        else
        {
            alpha[i] = min(1.0,sum_neg);
            sum_neg -= alpha[i];
        }
    
    double *zeros = new double[l];
    
    for(i=0;i<l;i++)
        zeros[i] = 0;
    
    Solver_NU s;
    s.Solve(l, SVC_Q(*prob,*param,y), zeros, y,
            alpha, 1.0, 1.0, param->eps, si,  param->shrinking);
    double r = si->r;
    
    info("C = %f\n",1/r);
    
    for(i=0;i<l;i++)
        alpha[i] *= y[i]/r;
    
    si->rho /= r;
    si->obj /= (r*r);
    si->upper_bound_p = 1/r;
    si->upper_bound_n = 1/r;
    
    delete[] y;
    delete[] zeros;
}

static void solve_one_class(
                            const svm_problem *prob, const svm_parameter *param,
                            double *alpha, Solver::SolutionInfo* si)
{
    int l = prob->l;
    double *zeros = new double[l];
    schar *ones = new schar[l];
    int i;
    
    int n = (int)(param->nu*prob->l);	// # of alpha's at upper bound
    
    for(i=0;i<n;i++)
        alpha[i] = 1;
    if(n<prob->l)
        alpha[n] = param->nu * prob->l - n;
    for(i=n+1;i<l;i++)
        alpha[i] = 0;
    
    for(i=0;i<l;i++)
    {
        zeros[i] = 0;
        ones[i] = 1;
    }
    
    Solver s;
    s.Solve(l, ONE_CLASS_Q(*prob,*param), zeros, ones,
            alpha, 1.0, 1.0, param->eps, si, param->shrinking);
    
    delete[] zeros;
    delete[] ones;
}

static void solve_epsilon_svr(
                              const svm_problem *prob, const svm_parameter *param,
                              double *alpha, Solver::SolutionInfo* si)
{
    int l = prob->l;
    double *alpha2 = new double[2*l];
    double *linear_term = new double[2*l];
    schar *y = new schar[2*l];
    int i;
    
    for(i=0;i<l;i++)
    {
        alpha2[i] = 0;
        linear_term[i] = param->p - prob->y[i];
        y[i] = 1;
        
        alpha2[i+l] = 0;
        linear_term[i+l] = param->p + prob->y[i];
        y[i+l] = -1;
    }
    
    Solver s;
    s.Solve(2*l, SVR_Q(*prob,*param), linear_term, y,
            alpha2, param->C, param->C, param->eps, si, param->shrinking);
    
    double sum_alpha = 0;
    for(i=0;i<l;i++)
    {
        alpha[i] = alpha2[i] - alpha2[i+l];
        sum_alpha += fabs(alpha[i]);
    }
    info("nu = %f\n",sum_alpha/(param->C*l));
    
    delete[] alpha2;
    delete[] linear_term;
    delete[] y;
}

static void solve_nu_svr(
                         const svm_problem *prob, const svm_parameter *param,
                         double *alpha, Solver::SolutionInfo* si)
{
    int l = prob->l;
    double C = param->C;
    double *alpha2 = new double[2*l];
    double *linear_term = new double[2*l];
    schar *y = new schar[2*l];
    int i;
    
    double sum = C * param->nu * l / 2;
    for(i=0;i<l;i++)
    {
        alpha2[i] = alpha2[i+l] = min(sum,C);
        sum -= alpha2[i];
        
        linear_term[i] = - prob->y[i];
        y[i] = 1;
        
        linear_term[i+l] = prob->y[i];
        y[i+l] = -1;
    }
    
    Solver_NU s;
    s.Solve(2*l, SVR_Q(*prob,*param), linear_term, y,
            alpha2, C, C, param->eps, si, param->shrinking);
    
    info("epsilon = %f\n",-si->r);
    
    for(i=0;i<l;i++)
        alpha[i] = alpha2[i] - alpha2[i+l];
    
    delete[] alpha2;
    delete[] linear_term;
    delete[] y;
}

//
// decision_function
//
struct decision_function
{
    double *alpha;
    double rho;
};

static decision_function svm_train_one(
                                       const svm_problem *prob, const svm_parameter *param,
                                       double Cp, double Cn)
{
    double *alpha = Malloc(double,prob->l);
    Solver::SolutionInfo si;
    switch(param->svm_type)
    {
        case C_SVC:
            solve_c_svc(prob,param,alpha,&si,Cp,Cn);
            break;
        case NU_SVC:
            solve_nu_svc(prob,param,alpha,&si);
            break;
        case ONE_CLASS:
            solve_one_class(prob,param,alpha,&si);
            break;
        case EPSILON_SVR:
            solve_epsilon_svr(prob,param,alpha,&si);
            break;
        case NU_SVR:
            solve_nu_svr(prob,param,alpha,&si);
            break;
    }
    
    info("obj = %f, rho = %f\n",si.obj,si.rho);
    
    // output SVs
    
    int nSV = 0;
    int nBSV = 0;
    for(int i=0;i<prob->l;i++)
    {
        if(fabs(alpha[i]) > 0)
        {
            ++nSV;
            if(prob->y[i] > 0)
            {
                if(fabs(alpha[i]) >= si.upper_bound_p)
                    ++nBSV;
            }
            else
            {
                if(fabs(alpha[i]) >= si.upper_bound_n)
                    ++nBSV;
            }
        }
    }
    
    info("nSV = %d, nBSV = %d\n",nSV,nBSV);
    
    decision_function f;
    f.alpha = alpha;
    f.rho = si.rho;
    return f;
}

// Platt's binary SVM Probablistic Output: an improvement from Lin et al.
static void sigmoid_train(
                          int l, const double *dec_values, const double *labels,
                          double& A, double& B)
{
    double prior1=0, prior0 = 0;
    int i;
    
    for (i=0;i<l;i++)
        if (labels[i] > 0) prior1+=1;
        else prior0+=1;
    
    int max_iter=100;	// Maximal number of iterations
    double min_step=1e-10;	// Minimal step taken in line search
    double sigma=1e-12;	// For numerically strict PD of Hessian
    double eps=1e-5;
    double hiTarget=(prior1+1.0)/(prior1+2.0);
    double loTarget=1/(prior0+2.0);
    double *t=Malloc(double,l);
    double fApB,p,q,h11,h22,h21,g1,g2,det,dA,dB,gd,stepsize;
    double newA,newB,newf,d1,d2;
    int iter;
    
    // Initial Point and Initial Fun Value
    A=0.0; B=log((prior0+1.0)/(prior1+1.0));
    double fval = 0.0;
    
    for (i=0;i<l;i++)
    {
        if (labels[i]>0) t[i]=hiTarget;
        else t[i]=loTarget;
        fApB = dec_values[i]*A+B;
        if (fApB>=0)
            fval += t[i]*fApB + log(1+exp(-fApB));
        else
            fval += (t[i] - 1)*fApB +log(1+exp(fApB));
    }
    for (iter=0;iter<max_iter;iter++)
    {
        // Update Gradient and Hessian (use H' = H + sigma I)
        h11=sigma; // numerically ensures strict PD
        h22=sigma;
        h21=0.0;g1=0.0;g2=0.0;
        for (i=0;i<l;i++)
        {
            fApB = dec_values[i]*A+B;
            if (fApB >= 0)
            {
                p=exp(-fApB)/(1.0+exp(-fApB));
                q=1.0/(1.0+exp(-fApB));
            }
            else
            {
                p=1.0/(1.0+exp(fApB));
                q=exp(fApB)/(1.0+exp(fApB));
            }
            d2=p*q;
            h11+=dec_values[i]*dec_values[i]*d2;
            h22+=d2;
            h21+=dec_values[i]*d2;
            d1=t[i]-p;
            g1+=dec_values[i]*d1;
            g2+=d1;
        }
        
        // Stopping Criteria
        if (fabs(g1)<eps && fabs(g2)<eps)
            break;
        
        // Finding Newton direction: -inv(H') * g
        det=h11*h22-h21*h21;
        dA=-(h22*g1 - h21 * g2) / det;
        dB=-(-h21*g1+ h11 * g2) / det;
        gd=g1*dA+g2*dB;
        
        
        stepsize = 1;		// Line Search
        while (stepsize >= min_step)
        {
            newA = A + stepsize * dA;
            newB = B + stepsize * dB;
            
            // New function value
            newf = 0.0;
            for (i=0;i<l;i++)
            {
                fApB = dec_values[i]*newA+newB;
                if (fApB >= 0)
                    newf += t[i]*fApB + log(1+exp(-fApB));
                else
                    newf += (t[i] - 1)*fApB +log(1+exp(fApB));
            }
            // Check sufficient decrease
            if (newf<fval+0.0001*stepsize*gd)
            {
                A=newA;B=newB;fval=newf;
                break;
            }
            else
                stepsize = stepsize / 2.0;
        }
        
        if (stepsize < min_step)
        {
            info("Line search fails in two-class probability estimates\n");
            break;
        }
    }
    
    if (iter>=max_iter)
        info("Reaching maximal iterations in two-class probability estimates\n");
    free(t);
}

static double sigmoid_predict(double decision_value, double A, double B)
{
    double fApB = decision_value*A+B;
    // 1-p used later; avoid catastrophic cancellation
    if (fApB >= 0)
        return exp(-fApB)/(1.0+exp(-fApB));
    else
        return 1.0/(1+exp(fApB)) ;
}

// Method 2 from the multiclass_prob paper by Wu, Lin, and Weng
static void multiclass_probability(int k, double **r, double *p)
{
    int t,j;
    int iter = 0, max_iter=max(100,k);
    double **Q=Malloc(double *,k);
    double *Qp=Malloc(double,k);
    double pQp, eps=0.005/k;
    
    for (t=0;t<k;t++)
    {
        p[t]=1.0/k;  // Valid if k = 1
        Q[t]=Malloc(double,k);
        Q[t][t]=0;
        for (j=0;j<t;j++)
        {
            Q[t][t]+=r[j][t]*r[j][t];
            Q[t][j]=Q[j][t];
        }
        for (j=t+1;j<k;j++)
        {
            Q[t][t]+=r[j][t]*r[j][t];
            Q[t][j]=-r[j][t]*r[t][j];
        }
    }
    for (iter=0;iter<max_iter;iter++)
    {
        // stopping condition, recalculate QP,pQP for numerical accuracy
        pQp=0;
        for (t=0;t<k;t++)
        {
            Qp[t]=0;
            for (j=0;j<k;j++)
                Qp[t]+=Q[t][j]*p[j];
            pQp+=p[t]*Qp[t];
        }
        double max_error=0;
        for (t=0;t<k;t++)
        {
            double error=fabs(Qp[t]-pQp);
            if (error>max_error)
                max_error=error;
        }
        if (max_error<eps) break;
        
        for (t=0;t<k;t++)
        {
            double diff=(-Qp[t]+pQp)/Q[t][t];
            p[t]+=diff;
            pQp=(pQp+diff*(diff*Q[t][t]+2*Qp[t]))/(1+diff)/(1+diff);
            for (j=0;j<k;j++)
            {
                Qp[j]=(Qp[j]+diff*Q[t][j])/(1+diff);
                p[j]/=(1+diff);
            }
        }
    }
    if (iter>=max_iter)
        info("Exceeds max_iter in multiclass_prob\n");
    for(t=0;t<k;t++) free(Q[t]);
    free(Q);
    free(Qp);
}

// Cross-validation decision values for probability estimates
static void svm_binary_svc_probability(
                                       const svm_problem *prob, const svm_parameter *param,
                                       double Cp, double Cn, double& probA, double& probB)
{
    int i;
    int nr_fold = 5;
    int *perm = Malloc(int,prob->l);
    double *dec_values = Malloc(double,prob->l);
    
    // random shuffle
    for(i=0;i<prob->l;i++) perm[i]=i;
    for(i=0;i<prob->l;i++)
    {
        int j = i+rand()%(prob->l-i);
        swap(perm[i],perm[j]);
    }
    for(i=0;i<nr_fold;i++)
    {
        int begin = i*prob->l/nr_fold;
        int end = (i+1)*prob->l/nr_fold;
        int j,k;
        struct svm_problem subprob;
        
        subprob.l = prob->l-(end-begin);
        subprob.x = Malloc(struct svm_node*,subprob.l);
        subprob.y = Malloc(double,subprob.l);
        
        k=0;
        for(j=0;j<begin;j++)
        {
            subprob.x[k] = prob->x[perm[j]];
            subprob.y[k] = prob->y[perm[j]];
            ++k;
        }
        for(j=end;j<prob->l;j++)
        {
            subprob.x[k] = prob->x[perm[j]];
            subprob.y[k] = prob->y[perm[j]];
            ++k;
        }
        int p_count=0,n_count=0;
        for(j=0;j<k;j++)
            if(subprob.y[j]>0)
                p_count++;
            else
                n_count++;
        
        if(p_count==0 && n_count==0)
            for(j=begin;j<end;j++)
                dec_values[perm[j]] = 0;
        else if(p_count > 0 && n_count == 0)
            for(j=begin;j<end;j++)
                dec_values[perm[j]] = 1;
        else if(p_count == 0 && n_count > 0)
            for(j=begin;j<end;j++)
                dec_values[perm[j]] = -1;
        else
        {
            svm_parameter subparam = *param;
            subparam.probability=0;
            subparam.C=1.0;
            subparam.nr_weight=2;
            subparam.weight_label = Malloc(int,2);
            subparam.weight = Malloc(double,2);
            subparam.weight_label[0]=+1;
            subparam.weight_label[1]=-1;
            subparam.weight[0]=Cp;
            subparam.weight[1]=Cn;
            struct svm_model *submodel = svm_train(&subprob,&subparam);
            for(j=begin;j<end;j++)
            {
                svm_predict_values(submodel,prob->x[perm[j]],&(dec_values[perm[j]]));
                // ensure +1 -1 order; reason not using CV subroutine
                dec_values[perm[j]] *= submodel->label[0];
            }
            svm_free_and_destroy_model(&submodel);
            svm_destroy_param(&subparam);
        }
        free(subprob.x);
        free(subprob.y);
    }
    sigmoid_train(prob->l,dec_values,prob->y,probA,probB);
    free(dec_values);
    free(perm);
}

// Return parameter of a Laplace distribution
static double svm_svr_probability(const svm_problem *prob, const svm_parameter *param)
{
    int i;
    int nr_fold = 5;
    double *ymv = Malloc(double,prob->l);
    double mae = 0;
    
    svm_parameter newparam = *param;
    newparam.probability = 0;
    svm_cross_validation(prob,&newparam,nr_fold,ymv);
    for(i=0;i<prob->l;i++)
    {
        ymv[i]=prob->y[i]-ymv[i];
        mae += fabs(ymv[i]);
    }
    mae /= prob->l;
    double std=sqrt(2*mae*mae);
    int count=0;
    mae=0;
    for(i=0;i<prob->l;i++)
        if (fabs(ymv[i]) > 5*std)
            count=count+1;
        else
            mae+=fabs(ymv[i]);
    mae /= (prob->l-count);
    info("Prob. model for test data: target value = predicted value + z,\nz: Laplace distribution e^(-|z|/sigma)/(2sigma),sigma= %g\n",mae);
    free(ymv);
    return mae;
}


// label: label name, start: begin of each class, count: #data of classes, perm: indices to the original data
// perm, length l, must be allocated before calling this subroutine
static void svm_group_classes(const svm_problem *prob, int *nr_class_ret, int **label_ret, int **start_ret, int **count_ret, int *perm)
{
    int l = prob->l;
    int max_nr_class = 16;
    int nr_class = 0;
    int *label = Malloc(int,max_nr_class);
    int *count = Malloc(int,max_nr_class);
    int *data_label = Malloc(int,l);
    int i;
    
    for(i=0;i<l;i++)
    {
        int this_label = (int)prob->y[i];
        int j;
        for(j=0;j<nr_class;j++)
        {
            if(this_label == label[j])
            {
                ++count[j];
                break;
            }
        }
        data_label[i] = j;
        if(j == nr_class)
        {
            if(nr_class == max_nr_class)
            {
                max_nr_class *= 2;
                label = (int *)realloc(label,max_nr_class*sizeof(int));
                count = (int *)realloc(count,max_nr_class*sizeof(int));
            }
            label[nr_class] = this_label;
            count[nr_class] = 1;
            ++nr_class;
        }
    }
    
    //
    // Labels are ordered by their first occurrence in the training set.
    // However, for two-class sets with -1/+1 labels and -1 appears first,
    // we swap labels to ensure that internally the binary SVM has positive data corresponding to the +1 instances.
    //
    if (nr_class == 2 && label[0] == -1 && label[1] == 1)
    {
        swap(label[0],label[1]);
        swap(count[0],count[1]);
        for(i=0;i<l;i++)
        {
            if(data_label[i] == 0)
                data_label[i] = 1;
            else
                data_label[i] = 0;
        }
    }
    
    int *start = Malloc(int,nr_class);
    start[0] = 0;
    for(i=1;i<nr_class;i++)
        start[i] = start[i-1]+count[i-1];
    for(i=0;i<l;i++)
    {
        perm[start[data_label[i]]] = i;
        ++start[data_label[i]];
    }
    start[0] = 0;
    for(i=1;i<nr_class;i++)
        start[i] = start[i-1]+count[i-1];
    
    *nr_class_ret = nr_class;
    *label_ret = label;
    *start_ret = start;
    *count_ret = count;
    free(data_label);
}

//
// Interface functions
//
svm_model *svm_train(const svm_problem *prob, const svm_parameter *param)
{
    svm_model *model = Malloc(svm_model,1);
    model->param = *param;
    model->free_sv = 0;	// XXX
    
    if(param->svm_type == ONE_CLASS ||
       param->svm_type == EPSILON_SVR ||
       param->svm_type == NU_SVR)
    {
        // regression or one-class-svm
        model->nr_class = 2;
        model->label = NULL;
        model->nSV = NULL;
        model->probA = NULL; model->probB = NULL;
        model->sv_coef = Malloc(double *,1);
        
        if(param->probability &&
           (param->svm_type == EPSILON_SVR ||
            param->svm_type == NU_SVR))
        {
            model->probA = Malloc(double,1);
            model->probA[0] = svm_svr_probability(prob,param);
        }
        
        decision_function f = svm_train_one(prob,param,0,0);
        model->rho = Malloc(double,1);
        model->rho[0] = f.rho;
        
        int nSV = 0;
        int i;
        for(i=0;i<prob->l;i++)
            if(fabs(f.alpha[i]) > 0) ++nSV;
        model->l = nSV;
        model->SV = Malloc(svm_node *,nSV);
        model->sv_coef[0] = Malloc(double,nSV);
        model->sv_indices = Malloc(int,nSV);
        int j = 0;
        for(i=0;i<prob->l;i++)
            if(fabs(f.alpha[i]) > 0)
            {
                model->SV[j] = prob->x[i];
                model->sv_coef[0][j] = f.alpha[i];
                model->sv_indices[j] = i+1;
                ++j;
            }
        
        free(f.alpha);
    }
    else
    {
        // classification
        int l = prob->l;
        int nr_class;
        int *label = NULL;
        int *start = NULL;
        int *count = NULL;
        int *perm = Malloc(int,l);
        
        // group training data of the same class
        svm_group_classes(prob,&nr_class,&label,&start,&count,perm);
        if(nr_class == 1)
            info("WARNING: training data in only one class. See README for details.\n");
        
        svm_node **x = Malloc(svm_node *,l);
        int i;
        for(i=0;i<l;i++)
            x[i] = prob->x[perm[i]];
        
        // calculate weighted C
        
        double *weighted_C = Malloc(double, nr_class);
        for(i=0;i<nr_class;i++)
            weighted_C[i] = param->C;
        for(i=0;i<param->nr_weight;i++)
        {
            int j;
            for(j=0;j<nr_class;j++)
                if(param->weight_label[i] == label[j])
                    break;
            if(j == nr_class)
                fprintf(stderr,"WARNING: class label %d specified in weight is not found\n", param->weight_label[i]);
            else
                weighted_C[j] *= param->weight[i];
        }
        
        // train k*(k-1)/2 models
        
        bool *nonzero = Malloc(bool,l);
        for(i=0;i<l;i++)
            nonzero[i] = false;
        decision_function *f = Malloc(decision_function,nr_class*(nr_class-1)/2);
        
        double *probA=NULL,*probB=NULL;
        if (param->probability)
        {
            probA=Malloc(double,nr_class*(nr_class-1)/2);
            probB=Malloc(double,nr_class*(nr_class-1)/2);
        }
        
        int p = 0;
        for(i=0;i<nr_class;i++)
            for(int j=i+1;j<nr_class;j++)
            {
                svm_problem sub_prob;
                int si = start[i], sj = start[j];
                int ci = count[i], cj = count[j];
                sub_prob.l = ci+cj;
                sub_prob.x = Malloc(svm_node *,sub_prob.l);
                sub_prob.y = Malloc(double,sub_prob.l);
                int k;
                for(k=0;k<ci;k++)
                {
                    sub_prob.x[k] = x[si+k];
                    sub_prob.y[k] = +1;
                }
                for(k=0;k<cj;k++)
                {
                    sub_prob.x[ci+k] = x[sj+k];
                    sub_prob.y[ci+k] = -1;
                }
                
                if(param->probability)
                    svm_binary_svc_probability(&sub_prob,param,weighted_C[i],weighted_C[j],probA[p],probB[p]);
                
                f[p] = svm_train_one(&sub_prob,param,weighted_C[i],weighted_C[j]);
                for(k=0;k<ci;k++)
                    if(!nonzero[si+k] && fabs(f[p].alpha[k]) > 0)
                        nonzero[si+k] = true;
                for(k=0;k<cj;k++)
                    if(!nonzero[sj+k] && fabs(f[p].alpha[ci+k]) > 0)
                        nonzero[sj+k] = true;
                free(sub_prob.x);
                free(sub_prob.y);
                ++p;
            }
        
        // build output
        
        model->nr_class = nr_class;
        
        model->label = Malloc(int,nr_class);
        for(i=0;i<nr_class;i++)
            model->label[i] = label[i];
        
        model->rho = Malloc(double,nr_class*(nr_class-1)/2);
        for(i=0;i<nr_class*(nr_class-1)/2;i++)
            model->rho[i] = f[i].rho;
        
        if(param->probability)
        {
            model->probA = Malloc(double,nr_class*(nr_class-1)/2);
            model->probB = Malloc(double,nr_class*(nr_class-1)/2);
            for(i=0;i<nr_class*(nr_class-1)/2;i++)
            {
                model->probA[i] = probA[i];
                model->probB[i] = probB[i];
            }
        }
        else
        {
            model->probA=NULL;
            model->probB=NULL;
        }
        
        int total_sv = 0;
        int *nz_count = Malloc(int,nr_class);
        model->nSV = Malloc(int,nr_class);
        for(i=0;i<nr_class;i++)
        {
            int nSV = 0;
            for(int j=0;j<count[i];j++)
                if(nonzero[start[i]+j])
                {
                    ++nSV;
                    ++total_sv;
                }
            model->nSV[i] = nSV;
            nz_count[i] = nSV;
        }
        
        info("Total nSV = %d\n",total_sv);
        
        model->l = total_sv;
        model->SV = Malloc(svm_node *,total_sv);
        model->sv_indices = Malloc(int,total_sv);
        p = 0;
        for(i=0;i<l;i++)
            if(nonzero[i])
            {
                model->SV[p] = x[i];
                model->sv_indices[p++] = perm[i] + 1;
            }
        
        int *nz_start = Malloc(int,nr_class);
        nz_start[0] = 0;
        for(i=1;i<nr_class;i++)
            nz_start[i] = nz_start[i-1]+nz_count[i-1];
        
        model->sv_coef = Malloc(double *,nr_class-1);
        for(i=0;i<nr_class-1;i++)
            model->sv_coef[i] = Malloc(double,total_sv);
        
        p = 0;
        for(i=0;i<nr_class;i++)
            for(int j=i+1;j<nr_class;j++)
            {
                // classifier (i,j): coefficients with
                // i are in sv_coef[j-1][nz_start[i]...],
                // j are in sv_coef[i][nz_start[j]...]
                
                int si = start[i];
                int sj = start[j];
                int ci = count[i];
                int cj = count[j];
                
                int q = nz_start[i];
                int k;
                for(k=0;k<ci;k++)
                    if(nonzero[si+k])
                        model->sv_coef[j-1][q++] = f[p].alpha[k];
                q = nz_start[j];
                for(k=0;k<cj;k++)
                    if(nonzero[sj+k])
                        model->sv_coef[i][q++] = f[p].alpha[ci+k];
                ++p;
            }
        
        free(label);
        free(probA);
        free(probB);
        free(count);
        free(perm);
        free(start);
        free(x);
        free(weighted_C);
        free(nonzero);
        for(i=0;i<nr_class*(nr_class-1)/2;i++)
            free(f[i].alpha);
        free(f);
        free(nz_count);
        free(nz_start);
    }
    return model;
}

// Stratified cross validation
void svm_cross_validation(const svm_problem *prob, const svm_parameter *param, int nr_fold, double *target)
{
    int i;
    int *fold_start;
    int l = prob->l;
    int *perm = Malloc(int,l);
    int nr_class;
    if (nr_fold > l)
    {
        nr_fold = l;
        fprintf(stderr,"WARNING: # folds > # data. Will use # folds = # data instead (i.e., leave-one-out cross validation)\n");
    }
    fold_start = Malloc(int,nr_fold+1);
    // stratified cv may not give leave-one-out rate
    // Each class to l folds -> some folds may have zero elements
    if((param->svm_type == C_SVC ||
        param->svm_type == NU_SVC) && nr_fold < l)
    {
        int *start = NULL;
        int *label = NULL;
        int *count = NULL;
        svm_group_classes(prob,&nr_class,&label,&start,&count,perm);
        
        // random shuffle and then data grouped by fold using the array perm
        int *fold_count = Malloc(int,nr_fold);
        int c;
        int *index = Malloc(int,l);
        for(i=0;i<l;i++)
            index[i]=perm[i];
        for (c=0; c<nr_class; c++)
            for(i=0;i<count[c];i++)
            {
                int j = i+rand()%(count[c]-i);
                swap(index[start[c]+j],index[start[c]+i]);
            }
        for(i=0;i<nr_fold;i++)
        {
            fold_count[i] = 0;
            for (c=0; c<nr_class;c++)
                fold_count[i]+=(i+1)*count[c]/nr_fold-i*count[c]/nr_fold;
        }
        fold_start[0]=0;
        for (i=1;i<=nr_fold;i++)
            fold_start[i] = fold_start[i-1]+fold_count[i-1];
        for (c=0; c<nr_class;c++)
            for(i=0;i<nr_fold;i++)
            {
                int begin = start[c]+i*count[c]/nr_fold;
                int end = start[c]+(i+1)*count[c]/nr_fold;
                for(int j=begin;j<end;j++)
                {
                    perm[fold_start[i]] = index[j];
                    fold_start[i]++;
                }
            }
        fold_start[0]=0;
        for (i=1;i<=nr_fold;i++)
            fold_start[i] = fold_start[i-1]+fold_count[i-1];
        free(start);
        free(label);
        free(count);
        free(index);
        free(fold_count);
    }
    else
    {
        for(i=0;i<l;i++) perm[i]=i;
        for(i=0;i<l;i++)
        {
            int j = i+rand()%(l-i);
            swap(perm[i],perm[j]);
        }
        for(i=0;i<=nr_fold;i++)
            fold_start[i]=i*l/nr_fold;
    }
    
    for(i=0;i<nr_fold;i++)
    {
        int begin = fold_start[i];
        int end = fold_start[i+1];
        int j,k;
        struct svm_problem subprob;
        
        subprob.l = l-(end-begin);
        subprob.x = Malloc(struct svm_node*,subprob.l);
        subprob.y = Malloc(double,subprob.l);
        
        k=0;
        for(j=0;j<begin;j++)
        {
            subprob.x[k] = prob->x[perm[j]];
            subprob.y[k] = prob->y[perm[j]];
            ++k;
        }
        for(j=end;j<l;j++)
        {
            subprob.x[k] = prob->x[perm[j]];
            subprob.y[k] = prob->y[perm[j]];
            ++k;
        }
        struct svm_model *submodel = svm_train(&subprob,param);
        if(param->probability &&
           (param->svm_type == C_SVC || param->svm_type == NU_SVC))
        {
            double *prob_estimates=Malloc(double,svm_get_nr_class(submodel));
            for(j=begin;j<end;j++)
                target[perm[j]] = svm_predict_probability(submodel,prob->x[perm[j]],prob_estimates);
            free(prob_estimates);
        }
        else
            for(j=begin;j<end;j++)
                target[perm[j]] = svm_predict(submodel,prob->x[perm[j]]);
        svm_free_and_destroy_model(&submodel);
        free(subprob.x);
        free(subprob.y);
    }
    free(fold_start);
    free(perm);
}


int svm_get_svm_type(const svm_model *model)
{
    return model->param.svm_type;
}

int svm_get_nr_class(const svm_model *model)
{
    return model->nr_class;
}

void svm_get_labels(const svm_model *model, int* label)
{
    if (model->label != NULL)
        for(int i=0;i<model->nr_class;i++)
            label[i] = model->label[i];
}

void svm_get_sv_indices(const svm_model *model, int* indices)
{
    if (model->sv_indices != NULL)
        for(int i=0;i<model->l;i++)
            indices[i] = model->sv_indices[i];
}

int svm_get_nr_sv(const svm_model *model)
{
    return model->l;
}

double svm_get_svr_probability(const svm_model *model)
{
    if ((model->param.svm_type == EPSILON_SVR || model->param.svm_type == NU_SVR) &&
        model->probA!=NULL)
        return model->probA[0];
    else
    {
        fprintf(stderr,"Model doesn't contain information for SVR probability inference\n");
        return 0;
    }
}

double svm_predict_values(const svm_model *model, const svm_node *x, double* dec_values)
{
    int i;
    if(model->param.svm_type == ONE_CLASS ||
       model->param.svm_type == EPSILON_SVR ||
       model->param.svm_type == NU_SVR)
    {
        double *sv_coef = model->sv_coef[0];
        double sum = 0;
        for(i=0;i<model->l;i++)
            sum += sv_coef[i] * Kernel::k_function(x,model->SV[i],model->param);
        sum -= model->rho[0];
        *dec_values = sum;
        
        if(model->param.svm_type == ONE_CLASS)
            return (sum>0)?1:-1;
        else
            return sum;
    }
    else
    {
        int nr_class = model->nr_class;
        int l = model->l;
        
        double *kvalue = Malloc(double,l);
        for(i=0;i<l;i++)
            kvalue[i] = Kernel::k_function(x,model->SV[i],model->param);
        
        int *start = Malloc(int,nr_class);
        start[0] = 0;
        for(i=1;i<nr_class;i++)
            start[i] = start[i-1]+model->nSV[i-1];
        
        int *vote = Malloc(int,nr_class);
        for(i=0;i<nr_class;i++)
            vote[i] = 0;
        
        int p=0;
        for(i=0;i<nr_class;i++)
            for(int j=i+1;j<nr_class;j++)
            {
                double sum = 0;
                int si = start[i];
                int sj = start[j];
                int ci = model->nSV[i];
                int cj = model->nSV[j];
                
                int k;
                double *coef1 = model->sv_coef[j-1];
                double *coef2 = model->sv_coef[i];
                for(k=0;k<ci;k++)
                    sum += coef1[si+k] * kvalue[si+k];
                for(k=0;k<cj;k++)
                    sum += coef2[sj+k] * kvalue[sj+k];
                sum -= model->rho[p];
                dec_values[p] = sum;
                
                if(dec_values[p] > 0)
                    ++vote[i];
                else
                    ++vote[j];
                p++;
            }
        
        int vote_max_idx = 0;
        for(i=1;i<nr_class;i++)
            if(vote[i] > vote[vote_max_idx])
                vote_max_idx = i;
        
        free(kvalue);
        free(start);
        free(vote);
        return model->label[vote_max_idx];
    }
}

double svm_predict(const svm_model *model, const svm_node *x)
{
    int nr_class = model->nr_class;
    double *dec_values;
    if(model->param.svm_type == ONE_CLASS ||
       model->param.svm_type == EPSILON_SVR ||
       model->param.svm_type == NU_SVR)
        dec_values = Malloc(double, 1);
    else
        dec_values = Malloc(double, nr_class*(nr_class-1)/2);
    double pred_result = svm_predict_values(model, x, dec_values);
    free(dec_values);
    return pred_result;
}

double svm_predict_probability(
                               const svm_model *model, const svm_node *x, double *prob_estimates)
{
    if ((model->param.svm_type == C_SVC || model->param.svm_type == NU_SVC) &&
        model->probA!=NULL && model->probB!=NULL)
    {
        int i;
        int nr_class = model->nr_class;
        double *dec_values = Malloc(double, nr_class*(nr_class-1)/2);
        svm_predict_values(model, x, dec_values);
        
        double min_prob=1e-7;
        double **pairwise_prob=Malloc(double *,nr_class);
        for(i=0;i<nr_class;i++)
            pairwise_prob[i]=Malloc(double,nr_class);
        int k=0;
        for(i=0;i<nr_class;i++)
            for(int j=i+1;j<nr_class;j++)
            {
                pairwise_prob[i][j]=min(max(sigmoid_predict(dec_values[k],model->probA[k],model->probB[k]),min_prob),1-min_prob);
                pairwise_prob[j][i]=1-pairwise_prob[i][j];
                k++;
            }
        multiclass_probability(nr_class,pairwise_prob,prob_estimates);
        
        int prob_max_idx = 0;
        for(i=1;i<nr_class;i++)
            if(prob_estimates[i] > prob_estimates[prob_max_idx])
                prob_max_idx = i;
        for(i=0;i<nr_class;i++)
            free(pairwise_prob[i]);
        free(dec_values);
        free(pairwise_prob);
        return model->label[prob_max_idx];
    }
    else
        return svm_predict(model, x);
}

static const char *svm_type_table[] =
{
    "c_svc","nu_svc","one_class","epsilon_svr","nu_svr",NULL
};

static const char *kernel_type_table[]=
{
    "linear","polynomial","rbf","sigmoid","precomputed",NULL
};

int svm_save_model(const char *model_file_name, const svm_model *model)
{
    FILE *fp = fopen(model_file_name,"w");
    if(fp==NULL) return -1;
    
    char *old_locale = strdup(setlocale(LC_ALL, NULL));
    setlocale(LC_ALL, "C");
    
    const svm_parameter& param = model->param;
    
    fprintf(fp,"svm_type %s\n", svm_type_table[param.svm_type]);
    fprintf(fp,"kernel_type %s\n", kernel_type_table[param.kernel_type]);
    
    if(param.kernel_type == POLY)
        fprintf(fp,"degree %d\n", param.degree);
    
    if(param.kernel_type == POLY || param.kernel_type == RBF || param.kernel_type == SIGMOID)
        fprintf(fp,"gamma %g\n", param.gamma);
    
    if(param.kernel_type == POLY || param.kernel_type == SIGMOID)
        fprintf(fp,"coef0 %g\n", param.coef0);
    
    int nr_class = model->nr_class;
    int l = model->l;
    fprintf(fp, "nr_class %d\n", nr_class);
    fprintf(fp, "total_sv %d\n",l);
    
    {
        fprintf(fp, "rho");
        for(int i=0;i<nr_class*(nr_class-1)/2;i++)
            fprintf(fp," %g",model->rho[i]);
        fprintf(fp, "\n");
    }
    
    if(model->label)
    {
        fprintf(fp, "label");
        for(int i=0;i<nr_class;i++)
            fprintf(fp," %d",model->label[i]);
        fprintf(fp, "\n");
    }
    
    if(model->probA) // regression has probA only
    {
        fprintf(fp, "probA");
        for(int i=0;i<nr_class*(nr_class-1)/2;i++)
            fprintf(fp," %g",model->probA[i]);
        fprintf(fp, "\n");
    }
    if(model->probB)
    {
        fprintf(fp, "probB");
        for(int i=0;i<nr_class*(nr_class-1)/2;i++)
            fprintf(fp," %g",model->probB[i]);
        fprintf(fp, "\n");
    }
    
    if(model->nSV)
    {
        fprintf(fp, "nr_sv");
        for(int i=0;i<nr_class;i++)
            fprintf(fp," %d",model->nSV[i]);
        fprintf(fp, "\n");
    }
    
    fprintf(fp, "SV\n");
    const double * const *sv_coef = model->sv_coef;
    const svm_node * const *SV = model->SV;
    
    for(int i=0;i<l;i++)
    {
        for(int j=0;j<nr_class-1;j++)
            fprintf(fp, "%.16g ",sv_coef[j][i]);
        
        const svm_node *p = SV[i];
        
        if(param.kernel_type == PRECOMPUTED)
            fprintf(fp,"0:%d ",(int)(p->value));
        else
            while(p->index != -1)
            {
                fprintf(fp,"%d:%.8g ",p->index,p->value);
                p++;
            }
        fprintf(fp, "\n");
    }
    
    setlocale(LC_ALL, old_locale);
    free(old_locale);
    
    if (ferror(fp) != 0 || fclose(fp) != 0) return -1;
    else return 0;
}

static char *line = NULL;
static int max_line_len;

static char* readline(FILE *input)
{
    int len;
    
    if(fgets(line,max_line_len,input) == NULL)
        return NULL;
    
    while(strrchr(line,'\n') == NULL)
    {
        max_line_len *= 2;
        line = (char *) realloc(line,max_line_len);
        len = (int) strlen(line);
        if(fgets(line+len,max_line_len-len,input) == NULL)
            break;
    }
    return line;
}

//
// FSCANF helps to handle fscanf failures.
// Its do-while block avoids the ambiguity when
// if (...)
//    FSCANF();
// is used
//
#define FSCANF(_stream, _format, _var) do{ if (fscanf(_stream, _format, _var) != 1) return false; }while(0)
bool read_model_header(FILE *fp, svm_model* model)
{
    svm_parameter& param = model->param;
    char cmd[81];
    while(1)
    {
        FSCANF(fp,"%80s",cmd);
        
        if(strcmp(cmd,"svm_type")==0)
        {
            FSCANF(fp,"%80s",cmd);
            int i;
            for(i=0;svm_type_table[i];i++)
            {
                if(strcmp(svm_type_table[i],cmd)==0)
                {
                    param.svm_type=i;
                    break;
                }
            }
            if(svm_type_table[i] == NULL)
            {
                fprintf(stderr,"unknown svm type.\n");
                return false;
            }
        }
        else if(strcmp(cmd,"kernel_type")==0)
        {
            FSCANF(fp,"%80s",cmd);
            int i;
            for(i=0;kernel_type_table[i];i++)
            {
                if(strcmp(kernel_type_table[i],cmd)==0)
                {
                    param.kernel_type=i;
                    break;
                }
            }
            if(kernel_type_table[i] == NULL)
            {
                fprintf(stderr,"unknown kernel function.\n");
                return false;
            }
        }
        else if(strcmp(cmd,"degree")==0)
            FSCANF(fp,"%d",&param.degree);
        else if(strcmp(cmd,"gamma")==0)
            FSCANF(fp,"%lf",&param.gamma);
        else if(strcmp(cmd,"coef0")==0)
            FSCANF(fp,"%lf",&param.coef0);
        else if(strcmp(cmd,"nr_class")==0)
            FSCANF(fp,"%d",&model->nr_class);
        else if(strcmp(cmd,"total_sv")==0)
            FSCANF(fp,"%d",&model->l);
        else if(strcmp(cmd,"rho")==0)
        {
            int n = model->nr_class * (model->nr_class-1)/2;
            model->rho = Malloc(double,n);
            for(int i=0;i<n;i++)
                FSCANF(fp,"%lf",&model->rho[i]);
        }
        else if(strcmp(cmd,"label")==0)
        {
            int n = model->nr_class;
            model->label = Malloc(int,n);
            for(int i=0;i<n;i++)
                FSCANF(fp,"%d",&model->label[i]);
        }
        else if(strcmp(cmd,"probA")==0)
        {
            int n = model->nr_class * (model->nr_class-1)/2;
            model->probA = Malloc(double,n);
            for(int i=0;i<n;i++)
                FSCANF(fp,"%lf",&model->probA[i]);
        }
        else if(strcmp(cmd,"probB")==0)
        {
            int n = model->nr_class * (model->nr_class-1)/2;
            model->probB = Malloc(double,n);
            for(int i=0;i<n;i++)
                FSCANF(fp,"%lf",&model->probB[i]);
        }
        else if(strcmp(cmd,"nr_sv")==0)
        {
            int n = model->nr_class;
            model->nSV = Malloc(int,n);
            for(int i=0;i<n;i++)
                FSCANF(fp,"%d",&model->nSV[i]);
        }
        else if(strcmp(cmd,"SV")==0)
        {
            while(1)
            {
                int c = getc(fp);
                if(c==EOF || c=='\n') break;
            }
            break;
        }
        else
        {
            fprintf(stderr,"unknown text in model file: [%s]\n",cmd);
            return false;
        }
    }
    
    return true;
    
}

svm_model *svm_load_model(const char *model_file_name)
{
    FILE *fp = fopen(model_file_name,"rb");
    if(fp==NULL) return NULL;
    
    char *old_locale = strdup(setlocale(LC_ALL, NULL));
    setlocale(LC_ALL, "C");
    
    // read parameters
    
    svm_model *model = Malloc(svm_model,1);
    model->rho = NULL;
    model->probA = NULL;
    model->probB = NULL;
    model->sv_indices = NULL;
    model->label = NULL;
    model->nSV = NULL;
    
    // read header
    if (!read_model_header(fp, model))
    {
        fprintf(stderr, "ERROR: fscanf failed to read model\n");
        setlocale(LC_ALL, old_locale);
        free(old_locale);
        free(model->rho);
        free(model->label);
        free(model->nSV);
        free(model);
        return NULL;
    }
    
    // read sv_coef and SV
    
    int elements = 0;
    long pos = ftell(fp);
    
    max_line_len = 1024;
    line = Malloc(char,max_line_len);
    char *p,*endptr,*idx,*val;
    
    while(readline(fp)!=NULL)
    {
        p = strtok(line,":");
        while(1)
        {
            p = strtok(NULL,":");
            if(p == NULL)
                break;
            ++elements;
        }
    }
    elements += model->l;
    
    fseek(fp,pos,SEEK_SET);
    
    int m = model->nr_class - 1;
    int l = model->l;
    model->sv_coef = Malloc(double *,m);
    int i;
    for(i=0;i<m;i++)
        model->sv_coef[i] = Malloc(double,l);
    model->SV = Malloc(svm_node*,l);
    svm_node *x_space = NULL;
    if(l>0) x_space = Malloc(svm_node,elements);
    
    int j=0;
    for(i=0;i<l;i++)
    {
        readline(fp);
        model->SV[i] = &x_space[j];
        
        p = strtok(line, " \t");
        model->sv_coef[0][i] = strtod(p,&endptr);
        for(int k=1;k<m;k++)
        {
            p = strtok(NULL, " \t");
            model->sv_coef[k][i] = strtod(p,&endptr);
        }
        
        while(1)
        {
            idx = strtok(NULL, ":");
            val = strtok(NULL, " \t");
            
            if(val == NULL)
                break;
            x_space[j].index = (int) strtol(idx,&endptr,10);
            x_space[j].value = strtod(val,&endptr);
            
            ++j;
        }
        x_space[j++].index = -1;
    }
    free(line);
    
    setlocale(LC_ALL, old_locale);
    free(old_locale);
    
    if (ferror(fp) != 0 || fclose(fp) != 0)
        return NULL;
    
    model->free_sv = 1;	// XXX
    return model;
}

void svm_free_model_content(svm_model* model_ptr)
{
    if(model_ptr->free_sv && model_ptr->l > 0 && model_ptr->SV != NULL)
        free((void *)(model_ptr->SV[0]));
    if(model_ptr->sv_coef)
    {
        for(int i=0;i<model_ptr->nr_class-1;i++)
            free(model_ptr->sv_coef[i]);
    }
    
    free(model_ptr->SV);
    model_ptr->SV = NULL;
    
    free(model_ptr->sv_coef);
    model_ptr->sv_coef = NULL;
    
    free(model_ptr->rho);
    model_ptr->rho = NULL;
    
    free(model_ptr->label);
    model_ptr->label= NULL;
    
    free(model_ptr->probA);
    model_ptr->probA = NULL;
    
    free(model_ptr->probB);
    model_ptr->probB= NULL;
    
    free(model_ptr->sv_indices);
    model_ptr->sv_indices = NULL;
    
    free(model_ptr->nSV);
    model_ptr->nSV = NULL;
}

void svm_free_and_destroy_model(svm_model** model_ptr_ptr)
{
    if(model_ptr_ptr != NULL && *model_ptr_ptr != NULL)
    {
        svm_free_model_content(*model_ptr_ptr);
        free(*model_ptr_ptr);
        *model_ptr_ptr = NULL;
    }
}

void svm_destroy_param(svm_parameter* param)
{
    free(param->weight_label);
    free(param->weight);
}

const char *svm_check_parameter(const svm_problem *prob, const svm_parameter *param)
{
    // svm_type
    
    int svm_type = param->svm_type;
    if(svm_type != C_SVC &&
       svm_type != NU_SVC &&
       svm_type != ONE_CLASS &&
       svm_type != EPSILON_SVR &&
       svm_type != NU_SVR)
        return "unknown svm type";
    
    // kernel_type, degree
    
    int kernel_type = param->kernel_type;
    if(kernel_type != LINEAR &&
       kernel_type != POLY &&
       kernel_type != RBF &&
       kernel_type != SIGMOID &&
       kernel_type != PRECOMPUTED)
        return "unknown kernel type";
    
    if(param->gamma < 0)
        return "gamma < 0";
    
    if(param->degree < 0)
        return "degree of polynomial kernel < 0";
    
    // cache_size,eps,C,nu,p,shrinking
    
    if(param->cache_size <= 0)
        return "cache_size <= 0";
    
    if(param->eps <= 0)
        return "eps <= 0";
    
    if(svm_type == C_SVC ||
       svm_type == EPSILON_SVR ||
       svm_type == NU_SVR)
        if(param->C <= 0)
            return "C <= 0";
    
    if(svm_type == NU_SVC ||
       svm_type == ONE_CLASS ||
       svm_type == NU_SVR)
        if(param->nu <= 0 || param->nu > 1)
            return "nu <= 0 or nu > 1";
    
    if(svm_type == EPSILON_SVR)
        if(param->p < 0)
            return "p < 0";
    
    if(param->shrinking != 0 &&
       param->shrinking != 1)
        return "shrinking != 0 and shrinking != 1";
    
    if(param->probability != 0 &&
       param->probability != 1)
        return "probability != 0 and probability != 1";
    
    if(param->probability == 1 &&
       svm_type == ONE_CLASS)
        return "one-class SVM probability output not supported yet";
    
    
    // check whether nu-svc is feasible
    
    if(svm_type == NU_SVC)
    {
        int l = prob->l;
        int max_nr_class = 16;
        int nr_class = 0;
        int *label = Malloc(int,max_nr_class);
        int *count = Malloc(int,max_nr_class);
        
        int i;
        for(i=0;i<l;i++)
        {
            int this_label = (int)prob->y[i];
            int j;
            for(j=0;j<nr_class;j++)
                if(this_label == label[j])
                {
                    ++count[j];
                    break;
                }
            if(j == nr_class)
            {
                if(nr_class == max_nr_class)
                {
                    max_nr_class *= 2;
                    label = (int *)realloc(label,max_nr_class*sizeof(int));
                    count = (int *)realloc(count,max_nr_class*sizeof(int));
                }
                label[nr_class] = this_label;
                count[nr_class] = 1;
                ++nr_class;
            }
        }
        
        for(i=0;i<nr_class;i++)
        {
            int n1 = count[i];
            for(int j=i+1;j<nr_class;j++)
            {
                int n2 = count[j];
                if(param->nu*(n1+n2)/2 > min(n1,n2))
                {
                    free(label);
                    free(count);
                    return "specified nu is infeasible";
                }
            }
        }
        free(label);
        free(count);
    }
    
    return NULL;
}

int svm_check_probability_model(const svm_model *model)
{
    return ((model->param.svm_type == C_SVC || model->param.svm_type == NU_SVC) &&
            model->probA!=NULL && model->probB!=NULL) ||
    ((model->param.svm_type == EPSILON_SVR || model->param.svm_type == NU_SVR) &&
     model->probA!=NULL);
}

void svm_set_print_string_function(void (*print_func)(const char *))
{
    if(print_func == NULL)
        svm_print_string = &print_string_stdout;
    else
        svm_print_string = print_func;
}

//-----------------------------End SVM---------------------------------------------------------

void readProblem(const VVD &samples, const VD &dv, struct svm_problem *prob, struct svm_node *x_space) {
    // calculate number of elements
    size_t features = samples[0].size();
    size_t elements = 0;
    for (const VD &row : samples) {
        for (const double &v : row) {
            if (v > 0) {
                elements++;
            }
        }
        // to include -1 at each row end
        elements++;
    }
    
    Printf("Non-zero elements: %lu, all elements: %lu\n", elements, features * samples.size());
    
    prob->l = (int)samples.size();
    prob->y = Malloc(double, prob->l);
    prob->x = Malloc(struct svm_node *, prob->l);
    x_space = Malloc(struct svm_node, elements);

    // fill with values
    int j = 0;
    for (int i = 0; i < prob->l; i++) {
        // set Y
        prob->y[i] = dv[i];
        
        // set X
        prob->x[i] = &x_space[j];
        for (int k = 0; k < features; k++) {
            if (samples[i][k] != 0) {
                // only non-zero added
                x_space[j].index = k + 1;// starts from 1
                x_space[j].value = samples[i][k];
                
                ++j;
            }
            //            Printf("%i : %.2f\n", x_space[j].index, x_space[j].value);
        }
        // mark end of row
        x_space[j++].index = -1;
    }
    Printf("Problem parsed, samples: %i, features: %i\n", prob->l, features);
}

//
//---------------------------------------------------------------------------------------------
//
//=============================================================================================
const static int SAMPLE_SIZE_HOR = 32;//2;//8;//16;//32;
const static int SAMPLE_SIZE_VER = 24;//48;
const static int SAMPLE_SIZE_MULT = SAMPLE_SIZE_HOR * SAMPLE_SIZE_VER;
const static int XSAMPLES = 640 / SAMPLE_SIZE_HOR;
const static int YSAMPLES = 480 / SAMPLE_SIZE_VER;

HoG hogoperator;

const static int HOG_WX = 5;
const static int HOG_WY = 5;
const static int HOG_BIN = 10;

void extractSampleHOG(const VI &img, const int x, const int y, VD &descriptor) {
    VVD res(SAMPLE_SIZE_VER, VD(SAMPLE_SIZE_HOR, 0));
    int index = x + y * 640;
    for (int j = 0; j < SAMPLE_SIZE_VER; j++) {
        VD row(SAMPLE_SIZE_HOR, 0);
        copy(img.begin() + index, img.begin() + index + SAMPLE_SIZE_HOR, row.begin());
        if (index + SAMPLE_SIZE_HOR > img.size()) {
            Assert(false, "Out of bounds");
        }
//        for (int i = 0; i < SAMPLE_SIZE_HOR; i++) {
//            row[i] = img[index + i] / 16777216.0;// grayscale
//        }
        index += 640;
        // add row
        res[j] = row;
    }
    
    // calculate HOG
    hogoperator.wx = HOG_WX;
    hogoperator.wy = HOG_WY;
    hogoperator.nbin = HOG_BIN;
    hogoperator.HOGdescriptor(res, descriptor);
}


void extractLabeledROISamples(const VI &img, const int ooiX, const int ooiY, VVD &features, VD &dv) {
    bool hasOOI = (ooiX > 0 && ooiY > 0);
    int sCount = 0, positives = 0;
    if (hasOOI) {
        // sample region of interest
        int startX = ooiX - SAMPLE_SIZE_HOR / 2;
        if (startX < 0) {
            startX += abs(startX);
        } else if(startX + SAMPLE_SIZE_HOR >= 640) {
            startX = 640 - SAMPLE_SIZE_HOR;
        }
        int startY = ooiY - SAMPLE_SIZE_VER / 2;
        if (startY < 0) {
            startY += abs(startY);
        } else if (startY + SAMPLE_SIZE_VER >= 480) {
            startY = 480 - SAMPLE_SIZE_VER;
        }
        VD sample;
        extractSampleHOG(img, startX, startY, sample);
        features.push_back(sample);
        dv.push_back(1.0);
        sCount++;
        positives++;
    } else {
        // false positives
        for (int ys = 0; ys < YSAMPLES; ys++) {
            for (int xs = 0; xs < XSAMPLES; xs++) {
                VD sample;
                extractSampleHOG(img, xs * SAMPLE_SIZE_HOR, ys * SAMPLE_SIZE_VER, sample);
                
                //            VD sample = extractSample(img, xs * SAMPLE_SIZE_HOR, ys * SAMPLE_SIZE_VER);
                features.push_back(sample);
                // not found
                dv.push_back(-1.0);
                sCount++;
            }
        }
    }
    Printf("Extracted %i samples with %i posititves\n", sCount, positives);
}

void extractROISamples(const VI &img, const int dx, const int dy, VVD &features) {
    int xcount = 640 / dx;
    int ycount = 480 / dy;
    VVD res(ycount, VD(xcount, 0));
    int sCount = 0;
    for (int ys = 0; ys < ycount; ys++) {
        int starty = ys * dy;
        if (starty + SAMPLE_SIZE_VER > 480) {
            starty = 480 - SAMPLE_SIZE_VER;
        }
        
        for (int xs = 0; xs < xcount; xs++) {
            VD sample;
            int startx = xs * dx;
            if (startx + SAMPLE_SIZE_HOR > 640) {
                startx = 640 - SAMPLE_SIZE_HOR;
            }
            extractSampleHOG(img, startx, starty, sample);
            features.push_back(sample);
            
            sCount++;
            
            if (startx + SAMPLE_SIZE_HOR == 640) {
                // outside
                break;
            }
        }
        
        if (starty + SAMPLE_SIZE_VER == 480) {
            // outside
            break;
        }
    }
    
    Printf("Sampled: %i regions from: %i\n", sCount, xcount * ycount);
}

pair<int, int>findMaximum(const VD &values, const int dx, const int dy) {
    int xcount = 640 / dx;
    int ycount = 480 / dy;
    int x = -1, y = -1, index = 0;
    int roiIndex = -1;
    double maxLabel = 0;
    for (int ys = 0; ys < ycount; ys++) {
        int starty = ys * dy;
        if (starty + SAMPLE_SIZE_VER > 480) {
            starty = 480 - SAMPLE_SIZE_VER;
        }
        
        for (int xs = 0; xs < xcount; xs++) {
            VD sample;
            int startx = xs * dx;
            if (startx + SAMPLE_SIZE_HOR > 640) {
                startx = 640 - SAMPLE_SIZE_HOR;
            }
            if (values[index] > maxLabel) {
                maxLabel = values[index];
                roiIndex = index;
                
                x = startx + SAMPLE_SIZE_HOR / 2;
                y = starty + SAMPLE_SIZE_VER / 2;
            }
            // increment
            index++;
            
            if (startx + SAMPLE_SIZE_HOR == 640) {
                // outside
                break;
            }
        }
        
        if (starty + SAMPLE_SIZE_VER == 480) {
            // outside
            break;
        }
    }
    if (roiIndex >= 0 && maxLabel > 0) {
        Printf("ROI at [%i, %i], index: %i, label: %.4f\n", x, y, roiIndex, maxLabel);
    } else {
        Printf("ROI not found---------------\n");
    }
    return pair<int, int>(x, y);
}


pair<int, int>detectRoi(const struct svm_model *model, const struct svm_problem *prob) {
    Printf("Start ROI with model: [classes: %i]\n", model->nr_class);
    int rows = prob->l;
    VD labels(rows, 0);
    int roiIndex = -1;
    double maxLabel = -1;
    for (int i = 0; i < rows; i++) {
        double* dec_values = (double *) malloc(1 * sizeof(double));
        double label = svm_predict_values(model, prob->x[i], dec_values);
        double decision = dec_values[0];
        labels[i] = label;
        if (decision > maxLabel) {
            maxLabel = decision;
            roiIndex = i;
        }
        free(dec_values);
    }
    int x = -1;
    int y = -1;
    if (roiIndex >= 0 && maxLabel > 0) {
        y = SAMPLE_SIZE_VER * floor(roiIndex / XSAMPLES);
        x = SAMPLE_SIZE_HOR * (roiIndex % XSAMPLES);
        
        y += SAMPLE_SIZE_VER / 2;
        x += SAMPLE_SIZE_HOR / 2;
        
        Printf("ROI at [%i, %i], index: %i, label: %i\n", x, y, roiIndex, maxLabel);
    }
    return pair<int, int>(x, y);
}

VD extractSample(const VI &img, const int x, const int y) {
    VD res(SAMPLE_SIZE_MULT, 0);
    int index = x + y * 640;
    int resIndex = 0;
    for (int j = 0; j < SAMPLE_SIZE_VER; j++) {
        for (int i = 0; i < SAMPLE_SIZE_HOR; i++) {
            res[resIndex++] = img[index + i] / 16777216.0;
        }
        index += 640;
    }
    return res;
}

void extractLabeledSamples(const VI &img, const int ooiX, const int ooiY, VVD &features, VD &dv) {
    bool hasOOI = (ooiX > 0 && ooiY > 0);
    int sCount = 0, positives = 0;
    for (int ys = 0; ys < YSAMPLES; ys++) {
        for (int xs = 0; xs < XSAMPLES; xs++) {
            VD sample;
            extractSampleHOG(img, xs * SAMPLE_SIZE_HOR, ys * SAMPLE_SIZE_VER, sample);
            
//            VD sample = extractSample(img, xs * SAMPLE_SIZE_HOR, ys * SAMPLE_SIZE_VER);
            features.push_back(sample);
            sCount++;
            if (hasOOI && (ooiX >= xs * SAMPLE_SIZE_HOR && ooiX < xs * SAMPLE_SIZE_HOR + SAMPLE_SIZE_HOR
                           && ooiY >= ys * SAMPLE_SIZE_VER && ooiY < ys * SAMPLE_SIZE_VER + SAMPLE_SIZE_VER)) {
                dv.push_back(1.0);
                positives++;
            } else {
                // not found
                dv.push_back(-1.0);
            }
        }
    }
    Printf("Extracted %i samples with %i posititves\n", sCount, positives);
}


void extractSamples(const VI &img, const int ooiX, const int ooiY, VVD &features, VD &dv) {
    bool hasOOI = (ooiX > 0 && ooiY > 0);
    int sCount = 0;
    for (int ys = 0; ys < YSAMPLES; ys++) {
        for (int xs = 0; xs < XSAMPLES; xs++) {
            VD sample;
            extractSampleHOG(img, xs * SAMPLE_SIZE_HOR, ys * SAMPLE_SIZE_VER, sample);
            features.push_back(sample);
            sCount++;
            if (hasOOI) {
                double dxv = (abs(ooiX - xs * SAMPLE_SIZE_HOR - SAMPLE_SIZE_HOR / 2) + 1);
                double xv = 640.0 / dxv;
                double yv = 480.0 / (abs(ooiY - ys * SAMPLE_SIZE_VER - SAMPLE_SIZE_VER / 2) + 1);
                dv.push_back((xv * yv) / 100);
            } else {
                // not found
                dv.push_back(-1.0);
            }
        }
    }
    Printf("Extracted %i samples\n", sCount);
}

double sampleRegion(const VI &img, const int x, const int y) {
    double res = 0;
    int index = x + y * 640;
    for (int j = 0; j < SAMPLE_SIZE_VER; j++) {
        for (int i = 0; i < SAMPLE_SIZE_HOR; i++) {
            res += img[index + i] / 16777216.0;
        }
        index += 640;
    }
    res /= SAMPLE_SIZE_MULT;
    return res;
}

void sampleRegions(const VI &img, const int ooiX, const int ooiY, VVD &features, VD &dv) {
    bool hasOOI = (ooiX > 0 && ooiY > 0);
    int sCount = 0;
    VD samples(YSAMPLES * XSAMPLES, 0);
    for (int ys = 0; ys < YSAMPLES; ys++) {
        for (int xs = 0; xs < XSAMPLES; xs++) {
            double sample = sampleRegion(img, xs * SAMPLE_SIZE_HOR, ys * SAMPLE_SIZE_VER);
            samples[sCount++] = sample;
        }
    }
    if (hasOOI) {
        double xv = round(100.0 * ooiX / 640.0);
        double yv = ooiY / 480.0;
        double v = (xv + yv) / 100.0;
        dv.push_back(v);
    } else {
        // not found
        dv.push_back(-1.0);
    }
    features.push_back(samples);
    
    //    Printf("Sampled %i regions\n", sCount);
}

pair<int, int>extractCoordinates(const double value) {
    int x = -1, y = -1;
    if (value > 0) {
        double v = value * 100.0;
        x = 640.0 * floor(v) / 100.0;
        y = 480.0 * (v - floor(v));
    }
    Printf("Extracted x: %i, y: %i\n", x, y);
    return pair<int, int>(x, y);
}

void sampleImageByHoG(const VI &img, const int ooiX, const int ooiY, VVD &features, VD &dv) {
    VVD res(480, VD(640, 0));
    int index = 0;
    for (int j = 0; j < 480; j++) {
        VD row(640, 0);
        copy(img.begin() + index, img.begin() + index + 640, row.begin());
        index += 640;
        // add row
        res[j] = row;
    }
    
    // calculate HOG
    VD descriptor;
    hogoperator.wx = HOG_WX;
    hogoperator.wy = HOG_WY;
    hogoperator.nbin = HOG_BIN;
    hogoperator.HOGdescriptor(res, descriptor);
    
    
    bool hasOOI = (ooiX > 0 && ooiY > 0);
    if (hasOOI) {
        double xv = round(100.0 * ooiX / 640.0);
        double yv = ooiY / 480.0;
        double v = (xv + yv) / 100.0;
        dv.push_back(v);
    } else {
        // not found
        dv.push_back(-1.0);
    }
    features.push_back(descriptor);
    
    //    Printf("Sampled %i regions\n", sCount);
}

pair<int, int>findMaximum(const VD &values) {
    int x = -1, y = -1, index = 0;
    int roiIndex = -1;
    double maxLabel = 0;
    for (int ys = 0; ys < YSAMPLES; ys++) {
        for (int xs = 0; xs < XSAMPLES; xs++) {
            if (values[index] > maxLabel) {
                maxLabel = values[index];
                roiIndex = index;
            }
            // increment
            index++;
        }
    }
    if (roiIndex >= 0 && maxLabel > 0) {
        y = SAMPLE_SIZE_VER * floor(roiIndex / XSAMPLES);
        x = SAMPLE_SIZE_HOR * (roiIndex % XSAMPLES);
        
        y += SAMPLE_SIZE_VER / 2;
        x += SAMPLE_SIZE_HOR / 2;
        
        Printf("ROI at [%i, %i], index: %i, label: %.4f\n", x, y, roiIndex, maxLabel);
    }
    return pair<int, int>(x, y);
}

void save_problem(const char *filename, const struct svm_problem *problem);
//
// ----------------------------
//
class RobotVisionTracker {
    VVD trainLeftFeatures;
    VD trainLeftDV;
    VVD trainRightFeatures;
    VD trainRightDV;
    
    RF_config conf;
    RF_Regression rfLeft;
    RF_Regression rfRight;
    
    int ooiCount = 0;
    int noOoiCount = 0;
    
    svm_model *modelLeft = NULL;
    svm_model *modelRight = NULL;
    
    svm_parameter param;
    
    svm_problem prob;
    svm_node *x_space = NULL;

public:
    RobotVisionTracker() {
        // just to be sure
        Assert(640 % SAMPLE_SIZE_HOR == 0, "Wrong horizontal sample");
        Assert(480 % SAMPLE_SIZE_VER == 0, "Wrong vertical sample");
        
        // init parameters
        param.svm_type = C_SVC;//
        param.kernel_type = RBF;
        param.degree = 3;
        param.gamma = 0;	// 1/num_features
        param.coef0 = 0;
        param.nu = 0.5;
        param.cache_size = 100;
        param.C = 1;
        param.eps = 1e-3;
        param.p = 0.1;
        param.shrinking = 1;
        param.probability = 0;//1;
        param.nr_weight = 0;
        param.weight_label = NULL;
        param.weight = NULL;
    }
    
    int training(const int videoIndex, const int frameIndex, const VI &imageDataLeft, const VI &imageDataRight, const int leftX, const int leftY, const int rightX, const int rightY) {
        
        Printf("+Adding training data: %i : %i, left[%i, %i], right[%i, %i]\n", videoIndex, frameIndex, leftX, leftY, rightX, rightY);
        
        // collect test data
        extractLabeledROISamples(imageDataLeft, leftX, leftY, trainLeftFeatures, trainLeftDV);
        extractLabeledROISamples(imageDataRight, rightX, rightY, trainRightFeatures, trainRightDV);
        
        
        
        if (leftX < 0 || leftY < 0) {
            noOoiCount++;
        } else {
            ooiCount++;
        }
        
        if (videoIndex == 1 /*&& frameIndex == 5*/) {
            return 1;
        }
        
        return 0;
    }
    
    VI testing(const int videoIndex, const int frameIndex, const VI &imageDataLeft, const VI &imageDataRight) {
        Printf("Test: %i : %i\n", videoIndex, frameIndex);
        
        // do left
        //
        VVD testFeatures;
        VD testDV;
//        extractLabeledSamples(imageDataLeft, -1, -1, testFeatures, testDV);
        extractROISamples(imageDataLeft, SAMPLE_SIZE_HOR / 4, SAMPLE_SIZE_VER / 4, testFeatures);
        
        VD res = rfLeft.predict(testFeatures, conf);
        pair<int, int> left = findMaximum(res, SAMPLE_SIZE_HOR / 4, SAMPLE_SIZE_VER / 4);
        
        // do right
        //
        testFeatures.clear();
        testDV.clear();
//        extractLabeledSamples(imageDataRight, -1, -1, testFeatures, testDV);
        extractROISamples(imageDataRight, SAMPLE_SIZE_HOR / 4, SAMPLE_SIZE_VER / 4, testFeatures);
        
        res = rfRight.predict(testFeatures, conf);
        pair<int, int> right = findMaximum(res, SAMPLE_SIZE_HOR / 4, SAMPLE_SIZE_VER / 4);
        
        VI result = {left.first, left.second, right.first, right.second};
        return result;
    }
    
    int doneTraining() {
        Printf("Frames with OOI: %i, without OOI: %i\n", ooiCount, noOoiCount);
        
        conf.nTree = 300;//500;
        conf.mtry = 60;
        //        conf.nodesize = 500;
        
        // do train
        rfLeft.train(trainLeftFeatures, trainLeftDV, conf);
        rfRight.train(trainRightFeatures, trainRightDV, conf);
        
        // release memory - waste of time - skipping
        //        trainLeftFeatures.clear();
        //        trainLeftDV.clear();
        //        trainRightFeatures.clear();
        //        trainRightDV.clear();
        
        return 0;
    }
    /*
    int training(const int videoIndex, const int frameIndex, const VI &imageDataLeft, const VI &imageDataRight, const int leftX, const int leftY, const int rightX, const int rightY) {
        
        Printf("Train: %i : %i, left[%i, %i], right[%i, %i]\n", videoIndex, frameIndex, leftX, leftY, rightX, rightY);
        
        // collect test data
        extractLabeledROISamples(imageDataLeft, leftX, leftY, trainLeftFeatures, trainLeftDV);
        extractLabeledROISamples(imageDataRight, rightX, rightY, trainRightFeatures, trainRightDV);
        
        
        
        if (leftX < 0 || leftY < 0) {
            noOoiCount++;
        } else {
            ooiCount++;
        }
        
//        if (videoIndex == 0 && frameIndex == 2) {
//            return 1;
//        }
        
        return 0;
    }
    
    VI testing(const int videoIndex, const int frameIndex, const VI &imageDataLeft, const VI &imageDataRight) {
        Printf("Test: %i : %i\n", videoIndex, frameIndex);
        
        // do left
        //
        VVD testFeaturesLeft;
        VD testDV;
        extractLabeledSamples(imageDataLeft, -1, -1, testFeaturesLeft, testDV);
        
        readProblem(testFeaturesLeft, testDV, &prob, x_space);
        pair<int, int>left = detectRoi(modelLeft, &prob);
        free(prob.y);
        free(prob.x);
        free(x_space);
        
        // do right
        //
        testDV.clear();
        VVD testFeaturesRight;
        extractLabeledSamples(imageDataRight, -1, -1, testFeaturesRight, testDV);
        
        readProblem(testFeaturesRight, testDV, &prob, x_space);
        pair<int, int>right = detectRoi(modelRight, &prob);
        free(prob.y);
        free(prob.x);
        free(x_space);
        
        VI result = {left.first, left.second, right.first, right.second};
        return result;
    }
    
    int doneTraining() {
        
        if (param.gamma == 0) {
            param.gamma = 1.0 / trainLeftFeatures[0].size();
        }
        
        Printf("Frames with OOI: %i, without OOI: %i, gamma: %.5f\n", ooiCount, noOoiCount, param.gamma);
        
        // do train left
        readProblem(trainLeftFeatures, trainLeftDV, &prob, x_space);
#ifdef SAVE_DATA
        save_problem("/Users/yaric/left_data.libsvm", &prob);
#endif
        
        if (!modelLeft) {
            svm_free_and_destroy_model(&modelLeft);
        }
        modelLeft = svm_train(&prob, &param);
        
        free(prob.y);
        free(prob.x);
        free(x_space);
        trainLeftFeatures.clear();
        trainLeftDV.clear();
        
        Printf("Left model: [classes: %i]\n", modelLeft->nr_class);
        int *labels = modelLeft->label;
        if (labels != NULL) {
            Printf("Left labels: ");
            for (int i = 0; i < sizeof(labels)/sizeof(labels[0]); i++) {
                Printf("%i, ", labels[i]);
            }
        }
        Printf("\n==========================\n");
        
        // do train right
        readProblem(trainRightFeatures, trainRightDV, &prob, x_space);
#ifdef SAVE_DATA
        save_problem("/Users/yaric/right_data.libsvm", &prob);
#endif
        
        if (!modelRight) {
            svm_free_and_destroy_model(&modelRight);
        }
        modelRight = svm_train(&prob, &param);
        
        free(prob.y);
        free(prob.x);
        free(x_space);
        trainRightFeatures.clear();
        trainRightDV.clear();
        
        Printf("Right model: [classes: %i]\n", modelRight->nr_class);
        labels = modelRight->label;
        if (labels != NULL) {
            Printf("Right labels: ");
            for (int i = 0; i < sizeof(labels)/sizeof(labels[0]); i++) {
                Printf("%i, ", labels[i]);
            }
        }
        Printf("\n==========================\n");
        
#ifdef SAVE_MODELS
        svm_save_model("/Users/yaric/left_model.svmmodel", modelLeft);
        svm_save_model("/Users/yaric/right_model.svmmodel", modelRight);
#endif
        
        // load models
     
//         modelLeft = svm_load_model("/Users/yaric/left_model.svmmodel");
//         modelRight = svm_load_model("/Users/yaric/right_model.svmmodel");
     
        return 0;
    }*/
};

void save_problem(const char *save_filename, const struct svm_problem *prob) {
    int l = prob->l;
    FILE *fp_save = fopen(save_filename,"w");
    if(fp_save==NULL) {
        fprintf(stderr,"can't open file %s\n", save_filename);
        exit(1);
    }
    for (int i = 0; i < l; i++) {
        // write line
        fprintf(fp_save, "%.16g ", prob->y[i]);
        const svm_node *p = prob->x[i];
        while (p->index != -1) {
            fprintf(fp_save,"%d:%.8g ", p->index, p->value);
            p++;
        }
        fprintf(fp_save, "\n");
    }
    
    fclose(fp_save);
}


#endif
