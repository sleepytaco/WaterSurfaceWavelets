#include "amplitude.h"

// the hat function is centered at centerVal and is zero outside minVal and maxVal
double hatFunction(int idx, double x, int dimX, double dX) {

    double centerVal = idx * dX; // value to center the hat (basis) function at
    double minVal = centerVal - dX;
    double maxVal = centerVal + dX;

    if (idx == 0) { // on the left boundary of sim domain the hat function is just the right half of the hat triangle
        if (x < 0) {
            return 0;
        } else {  // right half of the hat function
            return (maxVal - x) / (maxVal - centerVal);
        }
    } else if (idx == dimX) { // on the right boundary of sim domain the hat function is just the left half of the hat triangle
        if (x >= dimX*dX) {
            return 0;
        } else {  // left half of the hat function
            return (x - minVal) / (centerVal - minVal);
        }
    } else {
        if (x < minVal || x > maxVal) { // the hat function is 0 beyond a certain range
            return 0;
        } else if  (x < centerVal) {  // left half of the hat function
            return (x - minVal) / (centerVal - minVal);
        } else {  // right half of the hat function
            return (maxVal - x) / (maxVal - centerVal);
        }
    }

    print("forbidden area in hatFunction"); // have this here just in case none of the conditions above are executed (which should not happen)
}

// this reflects eqn 16 in the paper (the one using hat basis funcs to approximate A(x, theta, k))
double Amplitude::interpolateAmplitude4d(Vector2d x, double theta, double waveNumber) {
    double approxA = 0;
    // print("interpolate 4d");

    // this takes wayyy too slow to run as expected
    for (int i=0; i<=dimXY; ++i) { // a
        for (int j=0; j<=dimXY; ++j) { // a
            for (int b=0; b<=dimTheta; ++b) { // b
                for (int c=0; c<=dimK; ++c) { // c
                   // see eqn 16 for phi_a, nu_b, psi_c (which are basis functions)
                   double phi_a = hatFunction(i, x.x(), dimXY, dXY) * hatFunction(j, x.y(), dimXY, dXY);
                   double nu_b = hatFunction(b, theta, dimTheta, dTheta) ;
                   double psi_c = hatFunction(c, waveNumber, dimK, dK);

                   approxA += m_currentAmplitude.get(i, j, b, c) * phi_a * nu_b * psi_c;
                }
            }
        }
    }

    // TODO: working on narrowing down the search as most of the times phi_a, nu_b, psi_c are going to be 0s, except around
    // the inputs x, theta, waveNumber values

    return approxA;
}
