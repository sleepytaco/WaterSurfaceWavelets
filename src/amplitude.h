#include "Eigen/Dense"
#include "grid.h"
#include <iostream>
#include <unordered_set>
#include <deque>
#include <set>

using namespace Eigen;

// class to represent the discretized 4D amplitude function; this amplitude A is the simulation variable
class Amplitude
{
public:
    Amplitude();


    void setXMinMax(float xMin, float xMax);
    void setYMinMax(float yMin, float yMax);
    void setThetaMinMax(float thetaMin, float thetaMax);
    void setKMinMax(float kMin, float kMax);


    float getInterpolatedAmplitudeVal(Vector2f x, Vector2f k); // find the relevant
private:
    // index into x, y, theta, k

    Grid m_grid;


    // simulation domain range
    float xMin; float xMax;
    float yMin; float yMax;
    float thetaMin; float thetaMax;
    float kMin; float kMax;
};

