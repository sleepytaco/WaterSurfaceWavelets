#include "Eigen/Dense"

#include "grid.h"
#include "config.h"
#include "profilebuffer.h"

#include <iostream>
#include <optional>
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

//    void setXMinMax(float min, float max) { xMin = min; xMax = max; dXY = (xMax - xMin) / dimXY;};
//    void setYMinMax(float min, float max) { yMin = min; yMax = max; dXY = (yMax - yMin) / dimXY;};
//    void setThetaMinMax(float min, float max) { thetaMin = min; thetaMax = max; dTheta = (thetaMax - thetaMin) / dimTheta;};
//    void setKMinMax(float min, float max) { kMin = min; kMax = max; dK = (kMin - kMax) / 2;};

    Vector2d idxToPos(int i, int j);
    Vector2d posToIdxSpace(Vector2d pos); // returns (a, b) in "index space" eg. a = i + 0.123, b = j = 0.456

    double interpolateAmplitude(Vector2d idxSpacePos, int thetaIdx); // interpolate 4 nearest points

    // eqn 17 in paper
    double advectionSpeed(double waveNumber);

    double advectionAccel(double waveNumber);

    Vector2d advectionPos(Vector2d pos, double dt, double theta, double waveNumber);

    void boundaryReflection(Vector2d& advPos, int& thetaIdx);

    void advectionStep(double dt);

    double diffusionStep(double dt, double A, int xIdx, int yIdx, int thetaIdx, double waveNumber);

    void amplitudeSpread(double dt);

     // eqn 16 in paper
    double interpolateAmplitude4d(Vector2d x, double theta, double waveNumber);

    void precomputeProfileBuffers(double time);

    Vector3d waterHeight(Vector2d pos);

    void timeStep(double dt);

    void init_ls(std::vector<Vector3f>& terrain);

    float gamma(double waveNumber, double theta);

    float delta(double waveNumber, double theta);

    void toggleWind() { m_applyWindSim = !m_applyWindSim; };

    Grid m_currentAmplitude;
    Grid m_newAmplitude;

private:


    double bilerp(Vector2d x, double theta, double waveNumber);

    double catmullRom(std::vector<double>& segments, double adv_t);

    double spatialDiffusion(double dt, Vector2d idxPos, int xIdx, int yIdx, int thetaIdx, double waveNumber);


    //checking boundaries

    std::optional<Vector2d> boundary_check(Vector2d pos);



    double m_time = 0.0; // accumulate time across timesteps
    double m_timeToChangeWind = 3.0;
    double m_timeAcc = 0.0;

    bool m_applyWindSim = false;

    // number of samples for 4D amplitude grid
    const int dimXY = config.dimXY;
    const int dimTheta = config.dimTheta;
    const int dimK = config.dimK;

    // number of samples for integrating water height
    int numWaveNumberSamples = config.numWaveNumberSamples;
    int numThetaSamples = config.numThetaSamples;

    // simulation domain range
    // double dx, dy; // assuming this same as dXY
    double xMin=config.xMin; double xMax=config.xMax; double dXY=config.dXY; // the "width" of each cell XY grid cell
    double yMin=config.yMin; double yMax=config.yMax;
    double thetaMin=config.thetaMin; double thetaMax=config.thetaMax;  double dTheta=config.dTheta;
    double kMin=config.kMin; double kMax=config.kMax; double dK=config.dK;

    std::vector<Vector3f> terrain;

    ProfileBuffer m_profileBuffer;
};

