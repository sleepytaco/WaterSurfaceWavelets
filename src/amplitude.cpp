#include "amplitude.h"


Amplitude::Amplitude() {

}
Amplitude::Amplitude(int xSamples, int thetaSamples, int kSamples) : dimXY(xSamples), dimTheta(thetaSamples), dimK(kSamples) {
    m_currentAmplitude = Grid();
    m_newAmplitude = Grid();
    m_profileBuffer = ProfileBuffer();
}

// TODO
Vector2d Amplitude::idxToPos(int i, int j) {
    return Vector2d(0, 0);
}
// TODO
Vector2d Amplitude::posToIdxSpace(Vector2d pos) {
    return Vector2d(0, 0);
}

// assuming deep water dispersion: this is the equation for omega' if omega = sqrt(gk)
// TODO: this might be too simple and need to factor in additional terms
double Amplitude::advectionSpeed(double waveNumber) {
    return sqrt(9.8) * 0.5 / sqrt(waveNumber);
}

Vector2d Amplitude::advectionPos(Vector2d pos, double dt, double theta, double waveNumber) {
    Vector2d waveDirection = Vector2d(cos(theta), sin(theta));
    return pos - dt * advectionSpeed(waveNumber) * waveDirection;
}

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
                double amplitude = m_currentAmplitude(i, j, k, 0);
                double weight = (Vector3d((double)i, (double)j, (double)k) - Vector3d(idxSpaceX, idxSpaceY, idxSpaceTheta)).norm();
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
    interpolatedAmplitude *= m_profileBuffer.spectrum(m_profileBuffer.dispersion(waveNumber));

    return interpolatedAmplitude;
}

// advectionStep moves the current amplitudeGrid forward in time
void Amplitude::advectionStep(double dt) {
    for (int i=0; i<dimXY; ++i) { // a
        for (int j=0; j<dimXY; ++j) { // a
            for (int theta=0; theta<dimTheta; ++theta) { // b
                Vector2d x_a = idxToPos(i, j); // x_a = (x, y)
                double theta_b = theta*dTheta;
                double k_c = dK;

                Vector2d advPos = advectionPos(x_a, dt, theta_b, k_c); // pos "back in time" ---- x_jump, y_jump
                Vector2d idxSpaceAdvPos = posToIdxSpace(advPos);

                // fill in the newAmplitudeGrid with interpolated amplitude values at (x_jump, y_jump)
                m_newAmplitude(i, j, theta, 0) += interpolateAmplitude(idxSpaceAdvPos, theta);
            }
        }
    }

    std::swap(m_newAmplitude, m_currentAmplitude);
}

void Amplitude::precomputeProfileBuffers(double time) {
    m_profileBuffer.precompute(time);
}

double Amplitude::waterHeight(Vector2d pos) {
    double totalHeight = 0;

    for (int b = 1; b <= numThetaSamples; b++) {
        double theta = 2.0 * M_PI * (double)b / numThetaSamples;
        Vector2d waveDirection = Vector2d(cos(theta), sin(theta));
        double p = waveDirection.dot(pos);

        for (int c = 1; c <= numWaveNumberSamples; c++) {
            double fraction = (double)c / numWaveNumberSamples;
            double wavelength = wavelengthMax * fraction + wavelengthMin * (1 - fraction); // not 100% sure on this
            double waveNumber = 2.0 * M_PI / wavelength;
            totalHeight += interpolateAmplitude4d(pos, theta, waveNumber) * m_profileBuffer.getValueAt(p); // no shot this works first time. check here when things inevitably break
        }
    }

    return totalHeight;
}

void Amplitude::timeStep(double dt) {
    m_time += dt;
    advectionStep(dt);
    precomputeProfileBuffers(m_time);
}
