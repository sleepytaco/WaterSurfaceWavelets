#include "Eigen/Dense"
#include "profilebuffer.h"
#include <iostream>
#include <unordered_set>
#include <deque>
#include <set>
#include <vector>

#pragma once

using namespace Eigen;
using namespace std;

// class to represent the discretized 4D amplitude function; this amplitude A is the simulation variable
class Amplitude
{
public:
    Amplitude(int xSamples, int thetaSamples, int kSamples);

    void setXMinMax(float xMin, float xMax);
    void setYMinMax(float yMin, float yMax);
    void setThetaMinMax(float thetaMin, float thetaMax);
    void setKMinMax(float kMin, float kMax);

    Vector2d idxToPos(int i, int j);
    Vector2d posToIdxSpace(Vector2d pos); // returns (a, b) in "index space" eg. a = i + 0.123, b = j = 0.456

    double interpolateAmplitude(Vector2d idxSpacePos); // interpolate 4 nearest points

    // eqn 17 in paper
    double advectionSpeed(double waveNumber);
    Vector2d advectionPos(Vector2d pos, double dt, double theta, double waveNumber);
    void advectionStep(double dt);

     // eqn 16 in paper
    double getInterpolatedAmplitude(Vector2d x, double theta, double waveNumber);

    void precomputeProfileBuffers(double time);

    double waterHeight(Vector2d pos);

    void timeStep(double dt);

private:
    // index into x, y, theta, k
    float dx, dy;

    // samples
    int dimXY = 1024;
    int dimTheta = 16;
    int dimK = 1;

    float dXY = (xMax - xMin) / dimXY; // the "width" of each cell

    int w;
    int h;
    int k;
    int theta;

    int numWaveNumberSamples = 1;
    int numThetaSamples = 16; // number of samples for integrating water height

    double wavelengthMin = 0.02; // note: if these are changed so should the values in profilebuffer.h, these should be moved to a config file
    double wavelengthMax = 13.0;

    double m_time = 10.0; // accumulate time across timesteps, for some reason first few seconds are cursed so we skip to 10s

    vector<float> amplitudeGrid;
// TODO: make Grid4d class
//    Grid4d currAmplitude;
//    Grid4d newAmplitude;

    int gridIndex(size_t i1, size_t i2, size_t i3, size_t i4){ return i1 + w * (i2 + h * (i3 * k + i4));};

    // simulation domain range
    float xMin; float xMax;
    float yMin; float yMax;
    float thetaMin; float thetaMax;
    float kMin; float kMax;

    ProfileBuffer m_profileBuffer;
};

