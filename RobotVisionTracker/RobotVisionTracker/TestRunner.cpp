//
//  TestRunner.cpp
//  RobotVisionTracker
//
//  Created by Iaroslav Omelianenko on 3/30/15.
//  Copyright (c) 2015 Nologin. All rights reserved.
//

#include <stdio.h>
#include <ctime>
#include "lodepng/lodepng.h"

//#define X_CONTEST

#ifdef X_CONTEST
#include "RobotVisionTrackerX.h"
#else
#include "RobotVisionTracker.h"
#endif

inline int SPrintf(char *buf, size_t size, const char *fmt, ...) {
    va_list args;
    va_start(args, fmt);
    int ret = vsnprintf(buf, size, fmt, args);
    va_end(args);
    return ret;
}

VS splt(std::string s, char c = ',') {
    VS all;
    int p = 0, np;
    while (np = (int)s.find(c, p), np >= 0) {
        if (np != p)
            all.push_back(s.substr(p, np - p));
        else
            all.push_back("");
        
        p = np + 1;
    }
    if (p < s.size())
        all.push_back(s.substr(p));
    return all;
}

VI imageToArray(const char* filename) {
    std::vector<unsigned char> pixels;// raw pixels
    unsigned width, height;
    
    //decode
    unsigned error = lodepng::decode(pixels, width, height, filename);
    
    //if there's an error, display it
    if(error) std::cout << "decoder error " << error << ": " << lodepng_error_text(error) << std::endl;
    
    //the pixels are now in the vector "image", 4 bytes per pixel, ordered RGBARGBA..., use it as texture, draw it, ...
    VI res(width * height, 0);
    int pos = 0;
    for (int i = 0; i < pixels.size(); i += 4) {
        int v0 = (int)(pixels[i]);
        if (v0 < 0) v0 += 256;
        int v1 = (int)(pixels[i+1]);
        if (v1 < 0) v1 += 256;
        int v2 = (int)(pixels[i+2]);
        if (v2 < 0) v2 += 256;
        res[pos++] = v0 | (v1<<8) | (v2<<16);
    }
    return res;
}

struct VideoFrame {
    VI leftImg;
    VI rightImg;
    
    void Load(const string &fileNameLeft, const string &fileNameRight) {
        leftImg = imageToArray(fileNameLeft.c_str());
        rightImg = imageToArray(fileNameRight.c_str());
    }
};

struct VideoAnnotations {
    VI leftX;
    VI rightX;
    VI leftY;
    VI rightY;
    int numFrames;
    
    void Load(const string &fileName, bool bLeft) {
        
        if (bLeft) {
            leftX.clear();
            leftY.clear();
            rightX.clear();
            rightY.clear();
        }
        
        //
        // load data
        //
        std::ifstream datafile (fileName);
        if (!datafile.is_open()) {
            std::cerr << "Error in opening file: " << fileName << std::endl;
            return;
        }
        
        Printf("Loading annotation file: %s\n", fileName.c_str());
        
        // read annotations
        std::string sname;
        while (! datafile.eof() ) {
            getline(datafile, sname);
            if (strlen(sname.c_str()) == 0) break;
            
            int num = atoi(sname.substr(0, 5).c_str());
            while (leftX.size() <= num) {
                leftX.push_back(-1);
                leftY.push_back(-1);
                rightX.push_back(-1);
                rightY.push_back(-1);
            }
            std::string scnt;
            getline(datafile, scnt);
            std::string spnt;
            getline(datafile, spnt);
            VS items = splt(spnt, ' ');
            if (bLeft) {
                leftX[num] = atof(items[0].c_str());
                leftY[num] = atof(items[1].c_str());
            } else {
                rightX[num] = atof(items[0].c_str());
                rightY[num] = atof(items[1].c_str());
            }
        }
        numFrames = (int)leftX.size();
    }
};

class RobotVisionTrackerVis {
    const static int NUM_OF_TESTING_FRAMES = 50;
    const static int delay = 100;
    VS trainingVideos;
    VS testingVideos;
    
#ifdef X_CONTEST
    RobotVisionTrackerX *task;
#else
    RobotVisionTracker *task;
#endif
    
public:
    string folder;
#ifdef X_CONTEST
    RobotVisionTrackerVis (RobotVisionTrackerX *task) :task(task) {
#else
    RobotVisionTrackerVis (RobotVisionTracker *task) :task(task) {
#endif
        for (int i = 1; i < 16; i++) {
            trainingVideos.push_back("Closeup-Car" + to_string(i));
        }
        for (int i = 16; i < 21; i++) {
            testingVideos.push_back("Closeup-Car" + to_string(i));
        }
    }
    double doExec() {
        VideoFrame aFrame;
        VideoAnnotations aAnnotations;
        
        Printf("Training with %lu videos\n", trainingVideos.size());
        
        const char *fmt = "%05d.png";
        int sz = std::snprintf(nullptr, 0, fmt, std::sqrt(2));
        char buf[sz + 1]; // note +1 for null terminator
        
        int videoIndex = 0;
        for (string sVideoName : trainingVideos) {
            // load annotations
            aAnnotations.Load(folder + sVideoName + "_Left_annotation_pt.txt", true);
            aAnnotations.Load(folder + sVideoName + "_Right_annotation_pt.txt", false);

            bool bExitTraining = false;
            for (int f = 1; f < aAnnotations.numFrames; f++) {
                std::snprintf(&buf[0], sizeof(buf), fmt, f);
                
                string leftFilename = folder + sVideoName + "_Left_" + buf;
                string rightFilename = folder + sVideoName + "_Right_" + buf;
                
                aFrame.Load(leftFilename, rightFilename);
                
                int ret = task->training(videoIndex, f, aFrame.leftImg, aFrame.rightImg,
                                         aAnnotations.leftX[f], aAnnotations.leftY[f], aAnnotations.rightX[f], aAnnotations.rightY[f]);
                if (ret == 1) {
                    bExitTraining = true; // stop receiving training images
                    break;
                }
            }
            videoIndex++;
            if (bExitTraining) break;
        }
        
        // call doneTraining function
        task->doneTraining();
        
        // do test
        VideoAnnotations aTestAnnotations[5];
        videoIndex = 0;
        
        MT_RNG rnd;
        int N = (int)testingVideos.size();
        
        string sTestVideoName[5];
        int frameStart[5];
        for (int i = 0; i < 5; i++) {
            int iN = rnd.nextInt(N);
            sTestVideoName[i] = testingVideos[iN];
            // load annotations for video
            aTestAnnotations[i].Load(folder + sTestVideoName[i] + "_Left_annotation_pt.txt", true);
            aTestAnnotations[i].Load(folder + sTestVideoName[i] + "_Right_annotation_pt.txt", false);
            // pick random starting frame
            frameStart[i] = max(1, 1 + rnd.nextInt(aTestAnnotations[i].numFrames - NUM_OF_TESTING_FRAMES / 5));
            Printf("Testing with video '%s' starting at frame: %d\n", sTestVideoName[i].c_str(), frameStart[i]);
        }

        int gtfLeftX[NUM_OF_TESTING_FRAMES];
        int gtfLeftY[NUM_OF_TESTING_FRAMES];
        int gtfRightX[NUM_OF_TESTING_FRAMES];
        int gtfRightY[NUM_OF_TESTING_FRAMES];
        int userLeftX[NUM_OF_TESTING_FRAMES];
        int userLeftY[NUM_OF_TESTING_FRAMES];
        int userRightX[NUM_OF_TESTING_FRAMES];
        int userRightY[NUM_OF_TESTING_FRAMES];
        
        double score = 0;
        clock_t processingTime = 0;

        for (int f = 0; f < NUM_OF_TESTING_FRAMES; f++) {
            userLeftX[f] = userLeftY[f] = userRightX[f] = userRightY[f] = -1;
        }
        for (int f = 0; f < NUM_OF_TESTING_FRAMES; f++) {
            int videoIndex = f / 10;
            int frameIndex = f % 10;
            int index = frameIndex + frameStart[videoIndex];
            
            std::snprintf(&buf[0], sizeof(buf), fmt, index);
            string leftFilename = folder + sTestVideoName[videoIndex] + "_Left_" + buf;
            string rightFilename = folder + sTestVideoName[videoIndex] + "_Right_" + buf;
            aFrame.Load(leftFilename, rightFilename);
            
            gtfLeftX[f] = aTestAnnotations[videoIndex].leftX[index];
            gtfLeftY[f] = aTestAnnotations[videoIndex].leftY[index];
            gtfRightX[f] = aTestAnnotations[videoIndex].rightX[index];
            gtfRightY[f] = aTestAnnotations[videoIndex].rightY[index];
            
            clock_t startTime = clock();
            
            VI ret = task->testing(videoIndex, frameIndex, aFrame.leftImg, aFrame.rightImg);
            
            // store values
            userLeftX[f] = ret[0];
            userLeftY[f] = ret[1];
            userRightX[f] = ret[2];
            userRightY[f] = ret[3];
            
            clock_t endTime = clock();
            
            processingTime += endTime - startTime;
        }
        
        processingTime = 1000 * ((float)processingTime) / CLOCKS_PER_SEC;
        
        // calculate score
        correctOutOfBounds(gtfLeftX, gtfLeftY);
        correctOutOfBounds(gtfRightX, gtfRightY);
        
        double l10 = countFrameCorrect(gtfLeftX, gtfLeftY, userLeftX, userLeftY, 10 * 10);
        double l20 = countFrameCorrect(gtfLeftX, gtfLeftY, userLeftX, userLeftY, 20 * 20);
        double l50 = countFrameCorrect(gtfLeftX, gtfLeftY, userLeftX, userLeftY, 50 * 50);
        double r10 = countFrameCorrect(gtfRightX, gtfRightY, userRightX, userRightY, 10 * 10);
        double r20 = countFrameCorrect(gtfRightX, gtfRightY, userRightX, userRightY, 20 * 20);
        double r50 = countFrameCorrect(gtfRightX, gtfRightY, userRightX, userRightY, 50 * 50);
        
        Printf("R[10] = %.2f %.2f\n", l10, r10);
        Printf("R[20] = %.2f %.2f\n", l20, r20);
        Printf("R[50] = %.2f %.2f\n", l50, r50);
        
        double accuracyScore = 10000.0 * (50.0 * (l10 + r10) + 35.0 * (l20 + r20) + 15.0 * (l50 + r50));
        
        Printf("Accuracy Score = %.1f\n", accuracyScore);
        Printf("Processing Time = %d ms\n", processingTime);
        
        double timeMultiplier = 1.0;
        double T = 0.001 * processingTime;
        if (T > 3.33) timeMultiplier = 1.3536 - 0.2939 * log(T);
        if (T > 100.0) timeMultiplier = 0;
        timeMultiplier += 1.0;
        
        Printf("Time Multiplier = %.2f\n", timeMultiplier);
        
        score = timeMultiplier * accuracyScore;
        
        return score;
    }
    
    
private:
    void correctOutOfBounds(int *userX, int *userY) {
        for (int i = 0; i < sizeof(userX)/sizeof(userX[0]); i++) {
            if (userX[i] < 0 || userX[i] > 639 || userY[i] < 0 || userY[i] > 479) {
                userX[i] = -1;
                userY[i] = -1;
            }
        }
    }
    
    double countFrameCorrect(const int *gtfX, const int *gtfY, const int *userX, const int *userY, const int radi) {
        double p = 0;
        for (int i = 0; i < NUM_OF_TESTING_FRAMES; i++) {
            if (gtfX[i] < 0 && userX[i] < 0) {
                // correctly detected that object is not in view
                p += 1.0 / NUM_OF_TESTING_FRAMES;
                continue;
            }
            if (gtfX[i] < 0 && userX[i] >= 0) {
                // incorrectly detected an object
                continue;
            }
            if (gtfX[i] >= 0 && userX[i] < 0) {
                // incorrectly detected no object
                continue;
            }
            double dx = (gtfX[i] - userX[i]);
            double dy = (gtfY[i] - userY[i]);
            if (dx*dx+dy*dy <= radi) p += 1.0 / NUM_OF_TESTING_FRAMES;
        }
        return p;
    }
    
};

int main(int argc, const char * argv[]) {
    if (argc < 1) {
        printf("Usage: folder\n");
        return 0;
    }
#ifdef X_CONTEST
    RobotVisionTrackerX task;
#else
    RobotVisionTracker task;
#endif
    RobotVisionTrackerVis runner(&task);
    runner.folder = argv[1];
    
    double score = runner.doExec();
    fprintf(stderr, "Score = %.1f\n", score);
}
