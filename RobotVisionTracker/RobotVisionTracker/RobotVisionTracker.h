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
        // alloc on heap - may be too big
        double* X_test = (double *) calloc(n_size * p_size, sizeof(double));
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
        free(X_test);
        
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

//-----------------------------------LBP-------------------------------------------------------

/* Compare a value pointed to by 'ptr' to the 'center' value and
 * increment pointer. Comparison is made by masking the most
 * significant bit of an integer (the sign) and shifting it to an
 * appropriate position. */
#define compab_mask_inc(ptr,shift) { value |= ((unsigned int)(*center - *ptr - 1) & 0x80000000) >> (31-shift); ptr++; }
/* Compare a value 'val' to the 'center' value. */
#define compab_mask(val,shift) { value |= ((unsigned int)(*center - (val) - 1) & 0x80000000) >> (31-shift); }
/* Predicate 1 for the 3x3 neighborhood */
#define predicate 1
/* The number of bits */
#define bits 8

typedef struct
{
    int x,y;
} integerpoint;

typedef struct
{
    double x,y;
} doublepoint;

integerpoint points[bits];
doublepoint offsets[bits];

/*
 * Get a bilinearly interpolated value for a pixel.
 */
inline double interpolate_at_ptr(int* upperLeft, int i, int columns)
{
    double dx = 1-offsets[i].x;
    double dy = 1-offsets[i].y;
    return
    *upperLeft*dx*dy +
    *(upperLeft+1)*offsets[i].x*dy +
    *(upperLeft+columns)*dx*offsets[i].y +
    *(upperLeft+columns+1)*offsets[i].x*offsets[i].y;
}

/*
 * Calculate the point coordinates for circular sampling of the neighborhood.
 */
void calculate_points(void)
{
    double step = 2 * M_PI / bits, tmpX, tmpY;
    int i;
    for (i=0;i<bits;i++)
    {
        tmpX = predicate * cos(i * step);
        tmpY = predicate * sin(i * step);
        points[i].x = (int)tmpX;
        points[i].y = (int)tmpY;
        offsets[i].x = tmpX - points[i].x;
        offsets[i].y = tmpY - points[i].y;
        if (offsets[i].x < 1.0e-10 && offsets[i].x > -1.0e-10) /* rounding error */
            offsets[i].x = 0;
        if (offsets[i].y < 1.0e-10 && offsets[i].y > -1.0e-10) /* rounding error */
            offsets[i].y = 0;
        
        if (tmpX < 0 && offsets[i].x != 0)
        {
            points[i].x -= 1;
            offsets[i].x += 1;
        }
        if (tmpY < 0 && offsets[i].y != 0)
        {
            points[i].y -= 1;
            offsets[i].y += 1;
        }
    }
}

/*
 * Calculate the LBP histogram for an integer-valued image. This is an
 * optimized version of the basic 8-bit LBP operator. Note that this
 * assumes 4-byte integers. In some architectures, one must modify the
 * code to reflect a different integer size.
 *
 * img: the image data, an array of rows*columns integers arranged in
 * a horizontal raster-scan order
 * rows: the number of rows in the image
 * columns: the number of columns in the image
 * result: an array of 256 integers. Will hold the 256-bin LBP histogram.
 * interpolated: if != 0, a circular sampling of the neighborhood is
 * performed. Each pixel value not matching the discrete image grid
 * exactly is obtained using a bilinear interpolation. You must call
 * calculate_points (only once) prior to using the interpolated version.
 * return value: result
 */
void lbp_histogram(int* img, int rows, int columns, int* result, int interpolated) {
    int leap = columns*predicate;
    /*Set up a circularly indexed neighborhood using nine pointers.*/
    int
    *p0 = img,
    *p1 = p0 + predicate,
    *p2 = p1 + predicate,
    *p3 = p2 + leap,
    *p4 = p3 + leap,
    *p5 = p4 - predicate,
    *p6 = p5 - predicate,
    *p7 = p6 - leap,
    *center = p7 + predicate;
    unsigned int value;
    int pred2 = predicate << 1;
    int r,c;
    
    memset(result,0,sizeof(int)*256); /* Clear result histogram */
    
    if (!interpolated)
    {
        for (r=0;r<rows-pred2;r++)
        {
            for (c=0;c<columns-pred2;c++)
            {
                value = 0;
                
                /* Unrolled loop */
                compab_mask_inc(p0,0);
                compab_mask_inc(p1,1);
                compab_mask_inc(p2,2);
                compab_mask_inc(p3,3);
                compab_mask_inc(p4,4);
                compab_mask_inc(p5,5);
                compab_mask_inc(p6,6);
                compab_mask_inc(p7,7);
                center++;
                
                result[value]++; /* Increase histogram bin value */
            }
            p0 += pred2;
            p1 += pred2;
            p2 += pred2;
            p3 += pred2;
            p4 += pred2;
            p5 += pred2;
            p6 += pred2;
            p7 += pred2;
            center += pred2;
        }
    }
    else
    {
        p0 = center + points[5].x + points[5].y * columns;
        p2 = center + points[7].x + points[7].y * columns;
        p4 = center + points[1].x + points[1].y * columns;
        p6 = center + points[3].x + points[3].y * columns;
        
        for (r=0;r<rows-pred2;r++)
        {
            for (c=0;c<columns-pred2;c++)
            {
                value = 0;
                
                /* Unrolled loop */
                compab_mask_inc(p1,1);
                compab_mask_inc(p3,3);
                compab_mask_inc(p5,5);
                compab_mask_inc(p7,7);
                
                /* Interpolate corner pixels */
                compab_mask((int)(interpolate_at_ptr(p0,5,columns)+0.5),0);
                compab_mask((int)(interpolate_at_ptr(p2,7,columns)+0.5),2);
                compab_mask((int)(interpolate_at_ptr(p4,1,columns)+0.5),4);
                compab_mask((int)(interpolate_at_ptr(p6,3,columns)+0.5),6);
                p0++;
                p2++;
                p4++;
                p6++;
                center++;
                
                result[value]++;
            }
            p0 += pred2;
            p1 += pred2;
            p2 += pred2;
            p3 += pred2;
            p4 += pred2;
            p5 += pred2;
            p6 += pred2;
            p7 += pred2;
            center += pred2;
        }
    }
}


//=============================================================================================
typedef enum _FeatureDescriptor {
    HOG,
    LBP
}FeatureDescriptor;

// holds number of detections
int trueRoiCount;
// holds number of
int falseRoiCounts;

const static int SAMPLE_SIZE_HOR = 20;//32;//2;//8;//16;//32;
const static int SAMPLE_SIZE_VER = 20;//24;//48;
const static int SAMPLE_SIZE_MULT = SAMPLE_SIZE_HOR * SAMPLE_SIZE_VER;
const static int XSAMPLES = 640 / SAMPLE_SIZE_HOR;
const static int YSAMPLES = 480 / SAMPLE_SIZE_VER;

HoG hogoperator;

const static int HOG_WX = 5;
const static int HOG_WY = 5;
const static int HOG_BIN = 10;
FeatureDescriptor descrType = HOG;//LBP;//

void extractSampleLBP(const VI &img, const int x, const int y, VD &descriptor) {
    VI sample(SAMPLE_SIZE_MULT, 0);
    int index = x + y * 640;
    for (int j = 0; j < SAMPLE_SIZE_VER; j++) {
        copy(img.begin() + index, img.begin() + index + SAMPLE_SIZE_HOR, sample.begin() + j * SAMPLE_SIZE_HOR);
        if (index + SAMPLE_SIZE_HOR > img.size()) {
            Assert(false, "Out of bounds");
        }

        index += 640;
    }
    int *result = (int *)calloc(256, sizeof(int));
    lbp_histogram(&sample[0], SAMPLE_SIZE_HOR, SAMPLE_SIZE_VER, result, 0);
    for (int i = 0; i < 256; i++) {
        descriptor.push_back(result[i] / 256.0);
    }
}

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

void extractSampleDescriptor(const VI &img, const int x, const int y, VD &descriptor) {
    if (descrType == HOG) {
        extractSampleHOG(img, x, y, descriptor);
    } else {
        extractSampleLBP(img, x, y, descriptor);
    }
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
        extractSampleDescriptor(img, startX, startY, sample);
        features.push_back(sample);
        dv.push_back(1.0);
        sCount++;
        positives++;
        
        trueRoiCount++;
    } else {
        // false positives
        for (int ys = 0; ys < YSAMPLES; ys++) {
            for (int xs = 0; xs < XSAMPLES; xs++) {
                if (xs % 2 == 0) {
                    // count each second
                    continue;
                }
                VD sample;
                extractSampleDescriptor(img, xs * SAMPLE_SIZE_HOR, ys * SAMPLE_SIZE_VER, sample);
                
                //            VD sample = extractSample(img, xs * SAMPLE_SIZE_HOR, ys * SAMPLE_SIZE_VER);
                features.push_back(sample);
                // not found
                dv.push_back(-1.0);
                sCount++;
                
                falseRoiCounts++;
            }
        }
    }
    Printf("Extracted %i samples with %i posititves\n", sCount, positives);
}

void extractSamplesByDescriptors(const VI &img, const int dx, const int dy, VVD &features) {
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
            extractSampleDescriptor(img, startx, starty, sample);
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

pair<int, int>findMaximumBySlidingWindow(const VD &values, const int dx, const int dy) {
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
                // look for the best match
                maxLabel = values[index];
                roiIndex = index;
                
                x = startx + SAMPLE_SIZE_HOR / 2;
                y = starty + SAMPLE_SIZE_VER / 2;
            }
            // increment
            index++;
            
            if (startx + SAMPLE_SIZE_HOR == 640 || index >= values.size()) {
                // outside
                break;
            }
        }
        
        if (starty + SAMPLE_SIZE_VER == 480 || index >= values.size()) {
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

void extractLabeledSamples(const VI &img, const int ooiX, const int ooiY, VVD &features, VD &dv) {
    bool hasOOI = (ooiX > 0 && ooiY > 0);
    int sCount = 0, positives = 0;
    for (int ys = 0; ys < YSAMPLES; ys++) {
        for (int xs = 0; xs < XSAMPLES; xs++) {
            VD sample;
            extractSampleDescriptor(img, xs * SAMPLE_SIZE_HOR, ys * SAMPLE_SIZE_VER, sample);
            
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

public:
    RobotVisionTracker() {
        // just to be sure
        Assert(640 % SAMPLE_SIZE_HOR == 0, "Wrong horizontal sample");
        Assert(480 % SAMPLE_SIZE_VER == 0, "Wrong vertical sample");
        
        trueRoiCount = 0;
        falseRoiCounts = 0;
        
        // initiate LPB
//        calculate_points();
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
        
//        if (videoIndex == 1 /*&& frameIndex == 5*/) {
//            return 1;
//        }
        
        return 0;
    }
    
    VI testing(const int videoIndex, const int frameIndex, const VI &imageDataLeft, const VI &imageDataRight) {
        Printf("Test: %i : %i\n", videoIndex, frameIndex);
        
        int dx = SAMPLE_SIZE_HOR / 4;
        int dy = SAMPLE_SIZE_VER / 4;
        
        Printf("Sliding step dx: %i, dy: %i\n", dx, dy);
        
        // do left
        //
        VVD testFeaturesLeft;
//        VD testDV;
//        extractLabeledSamples(imageDataLeft, -1, -1, testFeaturesLeft, testDV);
        extractSamplesByDescriptors(imageDataLeft, dx, dy, testFeaturesLeft);
        
        VD res = rfLeft.predict(testFeaturesLeft, conf);
        pair<int, int> left = findMaximumBySlidingWindow(res, dx, dy);
        
        // do right
        //
        VVD testFeaturesRight;
//        extractLabeledSamples(imageDataRight, -1, -1, testFeatures, testDV);
        extractSamplesByDescriptors(imageDataRight, dx, dy, testFeaturesRight);
        
        res = rfRight.predict(testFeaturesRight, conf);
        pair<int, int> right = findMaximumBySlidingWindow(res, dx, dy);
        
        VI result = {left.first, left.second, right.first, right.second};
        return result;
    }
    
    int doneTraining() {
        Printf("Frames with OOI: %i, without OOI: %i, true ROI number: %i, false ROI number: %i\n", ooiCount, noOoiCount, trueRoiCount, falseRoiCounts);
        
        conf.nTree = 200;//300;//
        conf.mtry = 60;//80;
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
};



#endif
