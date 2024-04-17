#include "Eigen/Dense"

#include "grid.h"

#include "profilebuffer.h"

#include <iostream>
#include <unordered_set>
#include <deque>
#include <set>
#include <vector>

#pragma once

using namespace Eigen;

// class to represent the discretized 4D amplitude function; this amplitude A is the simulation variable
class Amplitude
{
public:
    Amplitude();
    Amplitude(int xSamples, int thetaSamples, int kSamples);


    void setXMinMax(float min, float max) { xMin = min; xMax = max; dXY = (xMax - xMin) / dimXY;};
    void setYMinMax(float min, float max) { yMin = min; yMax = max; dXY = (yMax - yMin) / dimXY;};;
    void setThetaMinMax(float min, float max) { thetaMin = min; thetaMax = max; dTheta = (thetaMax - thetaMin) / dimTheta;};;
    void setKMinMax(float min, float max) { kMin = min; kMax = max; dK = (kMax - kMin) / dimK;};;

    Vector2d idxToPos(int i, int j);
    Vector2d posToIdxSpace(Vector2d pos); // returns (a, b) in "index space" eg. a = i + 0.123, b = j = 0.456

    double interpolateAmplitude(Vector2d idxSpacePos, int thetaIdx); // interpolate 4 nearest points

    // eqn 17 in paper
    double advectionSpeed(double waveNumber);
    Vector2d advectionPos(Vector2d pos, double dt, double theta, double waveNumber);
    void advectionStep(double dt);

     // eqn 16 in paper
    double interpolateAmplitude4d(Vector2d x, double theta, double waveNumber);

    void precomputeProfileBuffers(double time);

    double waterHeight(Vector2d pos);

    void timeStep(double dt);

private:

    Grid m_currentAmplitude;
    Grid m_newAmplitude;


    double m_time = 0.0; // accumulate time across timesteps

    // number of samples for 4D amplitude grid
    int dimXY = 1024;
    int dimTheta = 16;
    int dimK = 1;


    // simulation domain range
    // double dx, dy; // assuming this same as dXY
    double xMin; double xMax; double dXY; // the "width" of each cell XY grid cell
    double yMin; double yMax;
    double thetaMin; double thetaMax;  double dTheta;
    double kMin; double kMax; double dK;

    double wavelengthMin = 0.02; // note: if these are changed so should the values in profilebuffer.h, these should be moved to a config file
    double wavelengthMax = 13.0;

    // number of samples for integrating water height
    int numWaveNumberSamples = 1;
    int numThetaSamples = 16;

    ProfileBuffer m_profileBuffer;
};

