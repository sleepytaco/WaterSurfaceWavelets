#include "Eigen/Dense"
#include <iostream>
#include <unordered_set>
#include <deque>
#include <set>

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

    float getAmplitudeVal(int a, int b, int c) {return amplitudeGrid[a][a][b][c];}; // gets the amplitude value stored in the grid
    float getInterpolatedAmplitudeVal(Vector2f x, Vector2f k); // find the relevant
private:
    // index into x, y, theta, k
    vector<vector<vector<vector<float>>>> amplitudeGrid;

    // simulation domain range
    float xMin; float xMax;
    float yMin; float yMax;
    float thetaMin; float thetaMax;
    float kMin; float kMax;
};

