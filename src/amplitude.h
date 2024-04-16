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

    float getAmplitudeVal(Vector2i a, int b, int c); // gets the amplitude value stored in the grid
    float getInterpolatedAmplitudeVal(Vector2f x, Vector2f k); // find the relevant
private:
    // index into x, y, theta, k

    int w;
    int h;
    int k;
    int theta;


    vector<float> amplitudeGrid;


    int gridIndex(size_t i1, size_t i2, size_t i3, size_t i4){ return i1 + w * (i2 + h * (i3 * k + i4));};

    // simulation domain range
    float xMin; float xMax;
    float yMin; float yMax;
    float thetaMin; float thetaMax;
    float kMin; float kMax;
};

