#include "amplitude.h"

// TODO: potentially fix this to reflect eqn 16 in the paper (one with basis funcs)
double Amplitude::interpolateAmplitude4d(Vector2d x, double theta, double waveNumber) {
    Vector2d idxSpacePos = posToIdxSpace(x);
    double idxSpaceX = idxSpacePos.x();
    double idxSpaceY = idxSpacePos.y();
    double idxSpaceTheta = (float)dimTheta * theta / (2.0 * M_PI);

    int idxXMin = floor(idxSpaceX);
    int idxYMin = floor(idxSpaceY);
    int idxThetaMin = floor(idxSpaceTheta);

    std::vector<double> amplitudes;
    std::vector<double> weights;

    // Doing a similar thing to what we discussed with amplitude interpolation for advection for now
    // Aka: perform interpolation between the nearest 8 points
    for (int i = idxXMin; i <= idxXMin + 1; i++) {
        for (int j = idxYMin; j <= idxYMin + 1; j++) {
            for (int k = idxThetaMin; k <= idxThetaMin + 1; k++) {
                double amplitude = m_currentAmplitude.get(i, j, k, 0);
                double weight = (Vector3d(i, j, k) - Vector3d(idxSpaceX, idxSpaceY, idxSpaceTheta)).norm();
                amplitudes.push_back(amplitude);
                weights.push_back(weight);
            }
        }
    }

    double weightSum = 0;
    for (double weight: weights) {
        weightSum += weight;
    }
    for (double& weight: weights) {
        weight /= weightSum;
        assert(weight >= 0 && weight <= 1);
    }

    double interpolatedAmplitude = 0;
    for (int i = 0; i < amplitudes.size(); i++) {
        interpolatedAmplitude += amplitudes[i] * weights[i];
    }

    // part that I'm least sure about but I'm not sure what the basis function for wave number is supposed to be besides this
    //    interpolatedAmplitude *= m_profileBuffer.spectrum(m_profileBuffer.dispersion(waveNumber));

    return interpolatedAmplitude;
}

// TODO: turn this simple into cubic spatial interpolation instead of the simple adj amplitude samples weighted avg
double Amplitude::interpolateAmplitude(Vector2d idxSpacePos, int thetaIdx) {

    // obtain floor of x and y coordinates
    int i = floor(idxSpacePos[0]);
    int j = floor(idxSpacePos[1]);

    // obtain relevant values adjacent to idxSpacePos
    double topLeft     = m_currentAmplitude.get(i  , j  , thetaIdx, 0);
    double topRight    = m_currentAmplitude.get(i+1, j  , thetaIdx, 0);
    double bottomLeft  = m_currentAmplitude.get(i  , j+1, thetaIdx, 0);
    double bottomRight = m_currentAmplitude.get(i+1, j+1, thetaIdx, 0);

    // calculate relative weightings based on distance
    double weightTopLeft     = (idxSpacePos - Vector2d(i  , j  )).norm();
    double weightTopRight    = (idxSpacePos - Vector2d(i+1, j  )).norm();
    double weightBottomLeft  = (idxSpacePos - Vector2d(i  , j+1)).norm();
    double weightBottomRight = (idxSpacePos - Vector2d(i+1, j+1)).norm();

    // calculate total sum of weights, so we can normalize
    double normalizeFactor = weightTopLeft + weightTopRight + weightBottomLeft + weightBottomRight;

    // return weighted average of adjacent amplitudes
    return weightTopLeft     / normalizeFactor * topLeft +
           weightTopRight    / normalizeFactor * topRight +
           weightBottomLeft  / normalizeFactor * bottomLeft +
           weightBottomRight / normalizeFactor * bottomRight;
}

double Amplitude::bilerp(Vector2d x, double theta, double waveNumber){
    // obtain floor of x and y coordinates
    int i = floor(x[0]);
    int j = floor(x[1]);

    double topLeft     = m_currentAmplitude.get(i  , j  , theta, 0);
    double topRight    = m_currentAmplitude.get(i+1, j  , theta, 0);
    double bottomLeft  = m_currentAmplitude.get(i  , j+1, theta, 0);
    double bottomRight = m_currentAmplitude.get(i+1, j+1, theta, 0);


    Vector4d f = Vector4d(topLeft, topRight, bottomLeft, bottomRight);

    Matrix4d coords;
    coords<< (i + 1) * (j + 1), -(i + 1) * j, - i * (j + 1), i * j,
        -(j + 1), j , j + 1, -j,
        -(i + 1), i + 1, i, -i,
        1, -1, -1, 1;

    Vector4d weights = coords * f;
    double outp = weights(0) + weights(1) * x(0) + weights(2) * x(1) + weights(3) * x(0) * x(1);
    return outp;
}

double Amplitude::catmullRom(std::vector<double>& segments, double adv_t){
    double d_1 = (segments[2] - segments[0])/2.0;
    double d_2 = (segments[3] - segments[1])/2.0;
    double dq = segments[2] - segments[1];

    double interpolated = segments[1] + d_1 * adv_t + (3 * dq - 2 * d_1 - d_2) * adv_t * adv_t + (-2 * dq + d_1 + d_2) * adv_t * adv_t *adv_t;
    return interpolated;
}
