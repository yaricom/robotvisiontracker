
#ifndef RobotVisionTracker_AsteroidRejector_h
#define RobotVisionTracker_AsteroidRejector_h
#include <bits/stdtr1c++.h>
#include <sys/time.h>
//#include <emmintrin.h>

using namespace std;

#define INLINE   inline __attribute__ ((always_inline))
#define NOINLINE __attribute__ ((noinline))

#define ALIGNED __attribute__ ((aligned(16)))

#define likely(x)   __builtin_expect(!!(x),1)
#define unlikely(x) __builtin_expect(!!(x),0)

#define SSELOAD(a)     _mm_load_si128((__m128i*)&a)
#define SSESTORE(a, b) _mm_store_si128((__m128i*)&a, b)

#define FOR(i,a,b)  for(int i=(a);i<(b);++i)
#define REP(i,a)    FOR(i,0,a)
#define ZERO(m)     memset(m,0,sizeof(m))
#define ALL(x)      x.begin(),x.end()
#define PB          push_back
#define S           size()
#define LL          long long
#define ULL         unsigned long long
#define LD          long double
#define MP          make_pair
#define X           first
#define Y           second
#define VC          vector
#define PII         pair <int, int>
#define VI          VC < int >
#define VVI         VC < VI >
#define VVVI        VC < VVI >
#define VPII        VC < PII >
#define VD          VC < double >
#define VVD         VC < VD >
#define VS          VC < string >
#define VVS         VC < VS >
#define DB(a)       cerr << #a << ": " << (a) << endl;

template<class T> void print(VC < T > v) {cerr << "[";if (v.S) cerr << v[0];FOR(i, 1, v.S) cerr << ", " << v[i];cerr << "]" << endl;}
template<class T> string i2s(T x) {ostringstream o; o << x; return o.str();}
VS splt(string s, char c = ' ') {VS all; int p = 0, np; while (np = s.find(c, p), np >= 0) {if (np != p) all.PB(s.substr(p, np - p)); p = np + 1;} if (p < s.S) all.PB(s.substr(p)); return all;}

double getTime() {
    timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + 1e-6 * tv.tv_usec;
}

int THREADS_NO = 1;
bool PRINT_IMPORTANCES = false;
bool CHOOSE_IMPORTANCES = false;

struct RNG {
    unsigned int MT[624];
    int index;
    
    RNG(int seed = 1) {
        init(seed);
    }
    
    void init(int seed = 1) {
        MT[0] = seed;
        FOR(i, 1, 624) MT[i] = (1812433253UL * (MT[i-1] ^ (MT[i-1] >> 30)) + i);
        index = 0;
    }
    
    void generate() {
        const unsigned int MULT[] = {0, 2567483615UL};
        REP(i, 227) {
            unsigned int y = (MT[i] & 0x8000000UL) + (MT[i+1] & 0x7FFFFFFFUL);
            MT[i] = MT[i+397] ^ (y >> 1);
            MT[i] ^= MULT[y&1];
        }
        FOR(i, 227, 623) {
            unsigned int y = (MT[i] & 0x8000000UL) + (MT[i+1] & 0x7FFFFFFFUL);
            MT[i] = MT[i-227] ^ (y >> 1);
            MT[i] ^= MULT[y&1];
        }
        unsigned int y = (MT[623] & 0x8000000UL) + (MT[0] & 0x7FFFFFFFUL);
        MT[623] = MT[623-227] ^ (y >> 1);
        MT[623] ^= MULT[y&1];
    }
    
    unsigned int rand() {
        if (index == 0) {
            generate();
        }
        
        unsigned int y = MT[index];
        y ^= y >> 11;
        y ^= y << 7  & 2636928640UL;
        y ^= y << 15 & 4022730752UL;
        y ^= y >> 18;
        index = index == 623 ? 0 : index + 1;
        return y;
    }
    
    INLINE int next() {
        return rand();
    }
    
    INLINE int next(int x) {
        return rand() % x;
    }
    
    INLINE int next(int a, int b) {
        return a + (rand() % (b - a));
    }
    
    INLINE double nextDouble() {
        return (rand() + 0.5) * (1.0 / 4294967296.0);
    }
};

static RNG rng;


struct Detection {
    int ID;
    string detNo;
    int frameNo;
    int sourceNo;
    double date;
    double ra;
    double dec;
    double x;
    double y;
    double mag;
    double fwhm;
    double elong;
    double theta;
    double rmse;
    double delta;
    int truth;
    
    Detection(string s) {
        VS v = splt(s.c_str(), ' ');
        
        ID = atoi(v[0].c_str());
        detNo = v[1];
        frameNo = atoi(v[2].c_str());
        sourceNo = atoi(v[3].c_str());
        date = atof(v[4].c_str());
        ra = atof(v[5].c_str());
        dec = atof(v[6].c_str());
        x = atof(v[7].c_str());
        y = atof(v[8].c_str());
        mag = atof(v[9].c_str());
        fwhm = atof(v[10].c_str());
        elong = atof(v[11].c_str());
        theta = atof(v[12].c_str());
        rmse = atof(v[13].c_str());
        delta = atof(v[14].c_str());
        truth = v.S > 15 ? atoi(v[15].c_str()) : -1;
    }
};

struct TreeNode {
    int level;
    int feature;
    double value;
    double result;
    int left;
    int right;
    
    TreeNode() {
        level = -1;
        feature = -1;
        value = 0;
        result = 0;
        left = -1;
        right = -1;
    }
};

struct RandomForestConfig {
    VI randomFeatures = {5};
    VI randomPositions = {2};
    int featuresIgnored = 0;
    int groupFeature = -1; //NOT IMPLEMENTED
    int maxLevel = 100;
    int maxNodeSize = 1;
    int threadsNo = 1;
    double bagSize = 1.0;
    double timeLimit = 0;
    bool useBootstrapping = true;
    bool computeImportances = false;
    bool computeOOB = false; //NOT IMPLEMENTED
};

class DecisionTree {
public:
    VC < TreeNode > nodes;
    VD importances;
    
    DecisionTree() { }
    
    template < class T > DecisionTree(VC < VC < T > > &features, VC < T > &results, RandomForestConfig &config, int seed) {
        RNG r(seed);
        
        if (config.computeImportances) {
            importances = VD(features[0].S);
        }
        
        VI chosenGroups(features.S);
        REP(i, (int)(features.S * config.bagSize)) chosenGroups[r.next(features.S)]++;
        
        int bagSize = 0;
        REP(i, features.S) if (chosenGroups[i]) bagSize++;
        
        VI bag(bagSize);
        VI weight(features.S);
        
        int pos = 0;
        
        REP(i, features.S) {
            weight[i] = config.useBootstrapping ? chosenGroups[i] : min(1, chosenGroups[i]);
            if (chosenGroups[i]) bag[pos++] = i;
        }
        
        TreeNode root;
        root.level = 0;
        root.left = 0;
        root.right = pos;
        nodes.PB(root);
        
        VI stack;
        stack.PB(0);
        
        while (stack.S) {
            bool equal = true;
            
            int curNode = stack[stack.S - 1];
            stack.pop_back();
            
            int bLeft = nodes[curNode].left;
            int bRight = nodes[curNode].right;
            int bSize = bRight - bLeft;
            
            int totalWeight = 0;
            T totalSum = 0;
            FOR(i, bLeft, bRight) {
                totalSum += results[bag[i]] * weight[bag[i]];
                totalWeight += weight[bag[i]];
            }
            
            assert(bSize > 0);
            
            FOR(i, bLeft + 1, bRight) if (results[bag[i]] != results[bag[i - 1]]) {
                equal = false;
                break;
            }
            
            if (equal || bSize <= config.maxNodeSize || nodes[curNode].level >= config.maxLevel) {
                nodes[curNode].result = totalSum / totalWeight;
                continue;
            }
            
            int bestFeature = -1;
            int bestLeft = 0;
            int bestRight = 0;
            T bestValue = 0;
            T bestMSE = 1e99;
            
            const int randomFeatures = config.randomFeatures[min(nodes[curNode].level, (int)config.randomFeatures.S - 1)];
            REP(i, randomFeatures) {
                
                int featureID = config.featuresIgnored + r.next(features[0].S - config.featuresIgnored);
                
                T vlo, vhi;
                vlo = vhi = features[bag[bLeft]][featureID];
                FOR(j, bLeft + 1, bRight) {
                    vlo = min(vlo, features[bag[j]][featureID]);
                    vhi = max(vhi, features[bag[j]][featureID]);
                }
                if (vlo == vhi) continue;
                
                const int randomPositions = config.randomPositions[min(nodes[curNode].level, (int)config.randomPositions.S - 1)];
                REP(j, randomPositions) {
                    T splitValue = features[bag[bLeft + r.next(bSize)]][featureID];
                    if (splitValue == vlo) {
                        j--;
                        continue;
                    }
                    
                    T sumLeft = 0;
                    int totalLeft = 0;
                    int weightLeft = 0;
                    FOR(k, bLeft, bRight) {
                        int p = bag[k];
                        if (features[p][featureID] < splitValue) {
                            sumLeft += results[p] * weight[p];
                            weightLeft += weight[p];
                            totalLeft++;
                        }
                    }
                    
                    T sumRight = totalSum - sumLeft;
                    int weightRight = totalWeight - weightLeft;
                    int totalRight = bSize - totalLeft;
                    
                    if (totalLeft == 0 || totalRight == 0)
                        continue;
                    
                    T meanLeft = sumLeft / weightLeft;
                    T meanRight = sumRight / weightRight;
                    T error = 0;
                    //MCE
                    // error = (1 - meanLeft) * (1 - meanLeft) * (1 - meanLeft) * sumLeft +
                    // meanLeft * meanLeft * meanLeft * (weightLeft - sumLeft) +
                    // (1 - meanRight) * (1 - meanRight) * (1 - meanRight) * sumRight +
                    // meanRight * meanRight * meanRight * (weightRight - sumRight);
                    //MQE
                    error = (1 - meanLeft) * (1 - meanLeft) * (1 - meanLeft) * (1 - meanLeft) * sumLeft +
                    meanLeft * meanLeft * meanLeft * meanLeft * (weightLeft - sumLeft) +
                    (1 - meanRight) * (1 - meanRight) * (1 - meanRight) * (1 - meanRight) * sumRight +
                    meanRight * meanRight * meanRight * meanRight * (weightRight - sumRight);
                    
                    //MSE
                    // error = (1 - meanLeft) * (1 - meanLeft) * sumLeft +
                    // meanLeft * meanLeft * (weightLeft - sumLeft) +
                    // (1 - meanRight) * (1 - meanRight) * sumRight +
                    // meanRight * meanRight * (weightRight - sumRight);
                    
                    
                    if (error < bestMSE) {
                        bestMSE = error;
                        bestValue = splitValue;
                        bestFeature = featureID;
                        bestLeft = totalLeft;
                        bestRight = totalRight;
                        if (error == 0) goto outer;
                    }
                }
            }
        outer:
            
            if (bestLeft == 0 || bestRight == 0) {
                nodes[curNode].result = totalSum / totalWeight;
                continue;
            }
            
            if (config.computeImportances) {
                importances[bestFeature] += 1. * (bRight - bLeft) / features.S;
            }
            
            T mean = totalSum / totalWeight;
            
            T nextValue = -1e99;
            FOR(i, bLeft, bRight) if (features[bag[i]][bestFeature] < bestValue) nextValue = max(nextValue, features[bag[i]][bestFeature]);
            // assert(nextValue > -1e90);
            
            TreeNode left;
            TreeNode right;
            
            left.level = right.level = nodes[curNode].level + 1;
            nodes[curNode].feature = bestFeature;
            nodes[curNode].value = (bestValue + nextValue) / 2.0;
            if (!(nodes[curNode].value > nextValue)) nodes[curNode].value = bestValue;
            nodes[curNode].left = nodes.S;
            nodes[curNode].right = nodes.S + 1;
            
            int bMiddle = bRight;
            FOR(i, bLeft, bMiddle) {
                if (features[bag[i]][nodes[curNode].feature] >= nodes[curNode].value) {
                    swap(bag[i], bag[--bMiddle]);
                    i--;
                    continue;
                }
            }
            
            assert(bestLeft == bMiddle - bLeft);
            assert(bestRight == bRight - bMiddle);
            
            left.left = bLeft;
            left.right = bMiddle;
            right.left = bMiddle;
            right.right = bRight;
            
            stack.PB(nodes.S);
            stack.PB(nodes.S + 1);
            
            nodes.PB(left);
            nodes.PB(right);
            
        }
        
        nodes.shrink_to_fit();
    }
    
    template < class T > double estimate(VC < T > &features) {
        TreeNode *pNode = &nodes[0];
        while (true) {
            if (pNode->feature < 0) return pNode->result;
            pNode = &nodes[features[pNode->feature] < pNode->value ? pNode->left : pNode->right];
        }
    }
};

class RandomForest {
public:
    
    VC < DecisionTree > trees;
    VD importances;
    
    void clear() {
        trees.clear();
        trees.shrink_to_fit();
    }
    
    template < class T > void train(VC < VC < T > > &features, VC < T > &results, RandomForestConfig &config, int treesNo) {
        double startTime = getTime();
        if (config.threadsNo == 1) {
            REP(i, treesNo) {
                if (config.timeLimit && getTime() - startTime > config.timeLimit) break;
                trees.PB(DecisionTree(features, results, config, rng.next()));
            }
        } else {
            thread *threads = new thread[config.threadsNo];
            mutex mutex;
            REP(i, config.threadsNo)
            threads[i] = thread([&] {
                while (true) {
                    auto tree = DecisionTree(features, results, config, rng.next());
                    mutex.lock();
                    if (trees.S < treesNo)
                        trees.PB(tree);
                    bool done = trees.S >= treesNo || config.timeLimit && getTime() - startTime > config.timeLimit;
                    mutex.unlock();
                    if (done) break;
                }
            });
            REP(i, config.threadsNo) threads[i].join();
            delete[] threads;
        }
        
        if (config.computeImportances) {
            importances = VD(features[0].S);
            for (DecisionTree tree : trees)
                REP(i, importances.S)
                importances[i] += tree.importances[i];
            double sum = 0;
            REP(i, importances.S) sum += importances[i];
            REP(i, importances.S) importances[i] /= sum;
        }
    }
    
    template < class T > double estimate(VC < T > &features) {
        assert(trees.S);
        
        double sum = 0;
        REP(i, trees.S) sum += trees[i].estimate(features);
        return sum / trees.S;
    }
    
};

VVD trainFeatures;
VVD testFeatures;

template < class T > void removeElements(VC < T > &v, VI &pos) {
    sort(ALL(pos));
    reverse(ALL(pos));
    REP(i, pos.S) {
        swap(v[pos[i]], v[v.S-1]);
        v.pop_back();
    }
}

VD extractFeatures(VI &img, VS &detections) {
    VC < Detection > vdet;
    REP(i, detections.S) vdet.PB(Detection(detections[i]));
    
    int N = detections.S;
    
    VD feat;
    feat.PB(vdet[0].ID);
    feat.PB(vdet[0].truth);
    
    VVD av;
    
    if (false && img.S) {
        VC < double > colors;
        VC < double > contrast;
        
        REP(i,N) {
            double sum = 0;
            REP(j,64*64) sum += img[i*64*64+j];
            colors.PB(sum);
        }
        REP(i,N) {
            double sum = 0;
            REP(x,64) REP(y, 64) {
                if (y < 63) sum += abs(img[i*64*64+x*64+y]-img[i*64*64+x*64+y+1]);
                if (x < 63) sum += abs(img[i*64*64+x*64+y]-img[i*64*64+x*64+y+64]);
            }
            contrast.PB(sum);
        }
        sort(ALL(colors));
        sort(ALL(contrast));
        
        // REP(i,N) feat.PB(colors[i]);
        // REP(i,N) feat.PB(contrast[i]);
    }
    
    for (auto det : vdet) {
        VD v;
        v.PB(det.date);
        // v.PB(det.sourceNo);
        v.PB(det.ra);
        v.PB(det.dec);
        v.PB(det.x);
        v.PB(det.y);
        v.PB(det.mag);
        v.PB(det.fwhm);
        v.PB(det.elong);
        v.PB(det.theta);
        v.PB(det.rmse);
        v.PB(det.delta);
        av.PB(v);
    }
    
    int lastSource = 0;
    int sources = 0;
    int wrongSource = -1;
    REP(i, av.S) {
        if (vdet[i].sourceNo > 0) {
            sources++;
            lastSource = i;
        } else {
            wrongSource = i;
        }
    }
    assert(sources >= 3);
    feat.PB(wrongSource);
    feat.PB(wrongSource);
    feat.PB(wrongSource);
    // feat.PB(sources);
    // feat.PB(sources);
    
    double minx = +1e30, maxx = -1e30, miny = +1e30, maxy = -1e30;
    double minra = +1e30, maxra = -1e30, mindec = +1e30, maxdec = -1e30;
    double mindate = +1e30, maxdate = -1e30;
    REP(i, N) {
        minx = min(minx, vdet[i].x);
        maxx = max(maxx, vdet[i].x);
        miny = min(miny, vdet[i].y);
        maxy = max(maxy, vdet[i].y);
        minra = min(minra, vdet[i].ra);
        maxra = max(maxra, vdet[i].ra);
        mindec = min(mindec, vdet[i].dec);
        maxdec = max(maxdec, vdet[i].dec);
        mindate = min(mindate, vdet[i].date);
        maxdate = max(maxdate, vdet[i].date);
    }
    
    // feat.PB(sqrt((maxx-minx)*(maxx-minx)+(maxy-miny)*(maxy-miny))/(maxdate-mindate));
    feat.PB(sqrt((maxra-minra)*(maxra-minra)+(maxdec-mindec)*(maxdec-mindec))/(maxdate-mindate));
    // feat.PB(sqrt((maxra-minra)*(maxra-minra)+(maxdec-mindec)*(maxdec-mindec))/(maxdate-mindate));
    // feat.PB(sqrt((maxra-minra)*(maxra-minra)+(maxdec-mindec)*(maxdec-mindec))/(maxdate-mindate));
    
    
    if (sources == 3) {
        REP(i, av.S) if (vdet[i].sourceNo == 0) {
            REP(j, av[0].S) {
                double sum = 0;
                REP(k, av.S) if (k != i) sum += av[k][j];
                av[i][j] = sum / 3;
            }
            vdet[i].sourceNo = 1;
            sources++;
        }
    }
    
    REP(i, av[0].S) {
        double sum = 0;
        VD v;
        REP(j, av.S) {
            if (vdet[j].sourceNo == 0) continue;
            v.PB(av[j][i]);
            sum += av[j][i];
        }
        sort(ALL(v));
        double avg = sum / 4;
        double stddev = 0;
        REP(j, 4) stddev += (v[j] - avg) * (v[j] - avg);
        stddev = sqrt(stddev / 4);
        
        feat.PB(avg);
        
        feat.PB(v[0]);
        feat.PB(v[1]);
        feat.PB(v[2]);
        feat.PB(v[3]);
        
        feat.PB(av[0][i] - avg);
        // feat.PB(av[1][i] - avg);
        // feat.PB(av[2][i] - avg);
        feat.PB(av[3][i] - avg);
        
        feat.PB(v[0] - avg);
        // feat.PB(v[1] - avg);
        // feat.PB(v[2] - avg);
        feat.PB(v[3] - avg);
        
        // feat.PB(v[1] - v[0]);
        // feat.PB(v[2] - v[1]);
        // feat.PB(v[3] - v[2]);
        feat.PB(v[3] - v[0]);
        
        feat.PB(v[0] == 0 ? 0 : v[3] / v[0]);
        feat.PB(avg == 0 ? 0 : v[3] / avg);
        feat.PB(avg == 0 ? 0 : v[0] / avg);
        
        feat.PB(stddev);
    }
    
    // VPII pairs = {MP(6, 5), MP(7, 5), MP(7, 6), MP(9, 5), MP(9, 6), MP(9, 7), MP(10, 9)};
    // for (PII p : pairs) {
    // int i = p.X;
    // int j = p.Y;
    // double sum0 = 0;
    // double sum1 = 0;
    // REP(k, av.S) {
    // if (vdet[k].sourceNo == 0) continue;
    // sum0 += av[k][i];
    // sum1 += av[k][j];
    // }
    // sum0 /= sources;
    // sum1 /= sources;
    // feat.PB(sum0 * sum1);
    // feat.PB(sum0 / sum1);
    // }
    
    VI rem = {13, 15, 16, 17, 18, 19, 48, 49, 50, 51, 52, 58, 60, 62, 63, 64, 65, 66, 72, 73, 74, 81, 82, 95, 96, 100, 101, 102, 109, 114, 115, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 137, 138, 139, 140, 141, 142, 143, 144, 145, 151, 152, 153, 154, 155, 156, 157, 158, 159};
    // REP(i, rem.S) rem[i] += 4;
    removeElements(feat, rem);
    
    return feat;
}

void chooseImportances(VD &imp) {
    string s = "";
    FOR(i, 2, imp.S) {
        if (imp[i] < 0.0035) {
            if (s.S) s += ", ";
            s += i2s(i);
        }
    }
    s = "{" + s + "};";
    cerr << s << endl;
}

void printImportances(VD &imp) {
    int p = 0;
    p += 2;
    
    REP(i, 5) cerr << imp[p++] << ' ';
    cerr << endl;
    REP(i, 3) cerr << imp[p++] << ' ';
    cerr << endl;
    REP(j, 11) {
        REP(i, 9) fprintf(stderr, "%.5f ", imp[p+i]);
        cerr << endl;
        p += 9;
    }	
    
    // VPII pairs = {MP(6, 5), MP(7, 5), MP(7, 6), MP(9, 5), MP(9, 6), MP(9, 7), MP(10, 9)};
    // for (PII r : pairs) {
    // int j = r.X;
    // int k = r.Y;
    // REP(i, 2) fprintf(stderr, "%.5f ", imp[p+i]);
    // cerr << j << ' ' << k << endl;
    // p += 2;
    // }
    
    assert(p == imp.S);
    
}


VVD extractAllFeatures(VI &img, VS &detections) {
    unordered_map < int, VI > detID;
    REP(i, detections.S) {
        int id = atoi(splt(detections[i], ' ')[0].c_str());
        detID[id].PB(i);
    }
    
    VVD rv;
    
    for (auto p : detID) {
        VI nimg;
        nimg.reserve(64*64*p.Y.S);
        VS ndetections;
        for (int i : p.Y) {
            ndetections.PB(detections[i]);
            if (img.S) REP(j, 64*64) nimg.PB(img[i*64*64+j]);
        }
        VD nfeat = extractFeatures(nimg, ndetections);
        rv.PB(nfeat);
    }
    
    return rv;
}

class AsteroidRejector {public:	
    
    AsteroidRejector() {
        rng.init(1);
        trainFeatures.clear();
        testFeatures.clear();
    }
    
    int trainCall = 0;
    int testCall = 0;
    
    int trainingData(VI &imageData, VS &detections) {
        trainCall++;
        VVD newFeatures = extractAllFeatures(imageData, detections);
        trainFeatures.insert(trainFeatures.end(), ALL(newFeatures));
        return 0;
    }
    
    int testingData(VI &imageData, VS &detections) {
        testCall++;
        VVD newFeatures = extractAllFeatures(imageData, detections);
        testFeatures.insert(testFeatures.end(), ALL(newFeatures));
        return 0;
    }
    
    VI getAnswer() {
        VD truth(trainFeatures.S);
        REP(i, trainFeatures.S) truth[i] = trainFeatures[i][1];
        
        RandomForestConfig rfCfg;
        rfCfg.randomFeatures = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16};
        rfCfg.randomPositions = {2, 2, 3, 3, 4};
        rfCfg.featuresIgnored = 2;
        rfCfg.maxNodeSize = 1;
        rfCfg.bagSize = 3.0;
        rfCfg.timeLimit = 60 * 30;
        rfCfg.threadsNo = THREADS_NO;
        rfCfg.computeImportances = true;
        
        double startTime = getTime();
        RandomForest RF;
        
        RF.train(trainFeatures, truth, rfCfg, 4000);
        if (PRINT_IMPORTANCES) printImportances(RF.importances);
        if (CHOOSE_IMPORTANCES) chooseImportances(RF.importances);
        // cerr << "Train Time: " << (getTime() - startTime) << endl;
        
        VC < pair < double, int > > order;
        REP(i, testFeatures.S) order.PB(MP(RF.estimate(testFeatures[i]), (int)(testFeatures[i][0]+1e-9)));
        // REP(i, testFeatures.S) order.PB(MP(0, (int)(testFeatures[i][0]+1e-9)));
        sort(ALL(order));
        
        VI rv;
        REP(i, order.S) rv.PB(order[i].Y);
        // reverse(ALL(rv));
        
        return rv;
    }
    
};


#endif
