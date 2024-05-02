#include "amplitude.h"
#include <random>


Amplitude::Amplitude() {
    std::cout << "amplitude constructor" << std::endl;
    double lower_bound = -20;
    double upper_bound = 20;
    std::uniform_real_distribution<double> unif(lower_bound,upper_bound);
    std::default_random_engine re;

    m_currentAmplitude = Grid();
    m_newAmplitude = Grid();
    for (int i=0; i<=dimXY; ++i) { // a
        for (int j=0; j<=dimXY; ++j) { // a
            for (int theta=0; theta<=dimTheta; ++theta) { // b
                // uncomment to init a sqaure with 0 amplitude in the center of the grid
//                Vector2d x_a = idxToPos(i, j); // x_a = (x, y)
//                if ((x_a.x() >= 500 && x_a.x() <= 3500) && (x_a.y() >= 500 && x_a.y() <= 3500)) {
//                    //m_currentAmplitude.get(i, j, theta, 0) = unif(re); 20 * sin((i + j) / 1);
//                    continue;
//                }
                if ((i >= 50 && i <= 80 && j >= 50 && j <= 80))
//                if (j >= 100)
                    m_currentAmplitude.get(i, j, theta, 0) = 30; unif(re) * sin((i + j) / 2);
            }
        }
    }
    m_profileBuffer = ProfileBuffer();
}

Vector2d Amplitude::idxToPos(int i, int j) {
    double x = xMin + (i+0.5) * dXY;
    double y = yMin + (j+0.5) * dXY;
    return Vector2d(x, y);
}

Vector2d Amplitude::posToIdxSpace(Vector2d pos) {
    double idxSpaceX = (pos.x() - xMin) / dXY-0.5;
    double idxSpaceY = (pos.y() - yMin) / dXY-0.5;

//    if (idxSpaceX > dimXY) idxSpaceX -= dimXY;
//    if (idxSpaceY > dimXY) idxSpaceY -= dimXY;
//    if (idxSpaceX < 0) idxSpaceX += dimXY;
//    if (idxSpaceY < 0) idxSpaceY += dimXY;
//    if (idxSpaceX > dimXY) idxSpaceX = dimXY;
//    if (idxSpaceY > dimXY) idxSpaceY = dimXY;
//    if (idxSpaceX < 0) idxSpaceX = 0;
//    if (idxSpaceY < 0) idxSpaceY = 0;

    return Vector2d(idxSpaceX, idxSpaceY);
}


// Returns the new theta index after reflection
// Note: implementation details says to split the amplitude to the two nearest thetas to the reflected directions
// but not doing this now because it would require a more complicated structure
void Amplitude::boundaryReflection(Vector2d& advPos, int& thetaIdx) {
    // Only apply boundary conditions if on boundary of grid
    // TODO: maybe also apply boudnary conditions if water is hitting terrain
    Vector2d boundaryNormal = Vector2d(0, 0);
    double boundaryDistance = 0;
    if (advPos.x() < xMin + dXY) {
        boundaryNormal = Vector2d(1, 0);
        boundaryDistance = abs(advPos.x() - xMin);
    }
    else if (advPos.x() > xMax - dXY) {
        boundaryNormal = Vector2d(-1, 0);
        boundaryDistance = abs(advPos.x() - xMax);
    }
    else if (advPos.y() < yMin + dXY) {
        boundaryNormal = Vector2d(0, 1);
        boundaryDistance = abs(advPos.y() - yMin);
    }
    else if (advPos.y() > yMax - dXY) {
        boundaryNormal = Vector2d(0, -1);
        boundaryDistance = abs(advPos.y() - yMax);
    }
    else {
        return;
    }

    double theta = thetaIdx * dTheta;
    Vector2d waveDirection = Vector2d(cos(theta), sin(theta));

    Vector2d advPos_refl = advPos + 2 * boundaryDistance * boundaryNormal;
//    if (advPos_refl != advPos)
//        std::cout << advPos.x() << "," << advPos.y() << " -> " << advPos_refl.x() << "," << advPos_refl.y() << std::endl;
    advPos = advPos_refl;
    // assert(advPos.x() >= xMin && advPos.x() <= xMax && advPos.y() >= yMin && advPos.y() <= yMax);

    Vector2d k_refl = waveDirection - 2 * (boundaryNormal.dot(waveDirection)) * boundaryNormal;
    double theta_refl = atan2(k_refl.y(), k_refl.x());
//    std::cout << theta * 180 / M_PI << " -> " << theta_refl * 180 / M_PI << " " << boundaryNormal.x() << "," << boundaryNormal.y() << std::endl;
    if (theta_refl < 0) theta_refl += 2*M_PI;
    int thetaIdx_refl = floor(theta_refl / dTheta);
    thetaIdx = thetaIdx_refl;
}

// assuming deep water dispersion: this is the equation for omega' if omega = sqrt(gk)
// TODO: this might be too simple and need to factor in additional terms
double Amplitude::advectionSpeed(double waveNumber) {
    return 0.5 * sqrt(config.g / waveNumber) * 10;
//    double val = config.g/(2 * sqrt(config.g * waveNumber + config.sigma * pow(waveNumber, 3)));
//    return val;
}

double Amplitude::advectionAccel(double waveNumber) {
    //return sqrt(9.8) * 0.5 / sqrt(waveNumber);
    double val =-pow(config.g, 2)/(4 * pow(config.g * waveNumber + config.sigma * pow(waveNumber, 3), 1.5));
    return val;
}

Vector2d Amplitude::advectionPos(Vector2d pos, double dt, double theta, double waveNumber) {
    Vector2d waveDirection = Vector2d(cos(theta), sin(theta));
    Vector2d advPos = pos - dt * advectionSpeed(waveNumber) * waveDirection; // pos "back in time" ---- x_jump, y_jump
    return advPos;
}

// advectionStep moves the current amplitudeGrid forward in time
void Amplitude::advectionStep(double dt) {

    #pragma omp parallel for collapse(2)
    for (int i=0; i<=dimXY; ++i) { // a
        for (int j=0; j<=dimXY; ++j) { // a

            Vector2d x_a = idxToPos(i, j); // x_a = (x, y)
            double k_c = dK; // wave number

            for (int theta=0; theta<dimTheta; ++theta) { // b

                double theta_b = theta*dTheta;
                Vector2d advPos = advectionPos(x_a, dt, theta_b, k_c);// pos "back in time" ---- x_jump, y_jump
                int theta_refl = theta;
                boundaryReflection(advPos, theta_refl);
                Vector2d idxSpaceAdvPos = posToIdxSpace(advPos);
                double A = diffusionStep(dt, idxSpaceAdvPos, i, j, theta_refl, k_c);
//                std::cout << A << std::endl;
                m_newAmplitude.get(i, j, theta, 0) = A;

            }
        }
    }
    std::swap(m_newAmplitude, m_currentAmplitude);
}

double Amplitude::diffusionStep(double dt, Vector2d idxSpaceAdvPos, int xIdx, int yIdx, int thetaIdx, double waveNumber) {
    // fill in the newAmplitudeGrid with interpolated amplitude values at (x_jump, y_jump)
    double A = bilerp(idxSpaceAdvPos, thetaIdx, waveNumber);
    //calulate diffusion here
    double dissapation = 0; //2*1e-6*k_c*k_c*A; // including this term from paper stabalizes the water to a sheet pretty quickly
    double gamma = 0.025 * advectionSpeed(waveNumber) * pow(dTheta, 2)/dXY;
    double prevThetaIdx = thetaIdx == 0 ? 15 : thetaIdx - 1;
    double nextThetaIdx = thetaIdx == 15 ? 0 : thetaIdx + 1;
    double prev = m_currentAmplitude.get(xIdx, yIdx, prevThetaIdx, 0);
    double next = m_currentAmplitude.get(xIdx, yIdx, nextThetaIdx, 0);
    double s_d = (1 - 2 * gamma * dt/pow(dTheta, 2)) * A + gamma * dt/pow(dTheta, 2) * (prev + next);
    return s_d - dissapation;
}

void Amplitude::precomputeProfileBuffers(double time) {
    m_profileBuffer.precompute(time);
}

Vector3d Amplitude::waterHeight(Vector2d pos) {
    Vector3d total = Vector3d(0, 0, 0);
    double totalHeight = 0;

    for (int b = 0; b < numThetaSamples; b++) {
        double theta = 2.0 * M_PI * (double)b / (double)numThetaSamples;
        Vector2d waveDirection = Vector2d(cos(theta), sin(theta));
        double p = waveDirection.dot(pos);
        Vector2d profile = m_profileBuffer.getValueAt(p);

        for (int c = 1; c <= numWaveNumberSamples; c++) {

            double fraction = (double)c / (double)numWaveNumberSamples;
            double wavelength = config.wavelengthMax * fraction + config.wavelengthMin * (1 - fraction); // not 100% sure on this
            double waveNumber = 2.0 * M_PI / wavelength;

//            double dWaveNumber = ((2.0 * M_PI)/config.wavelengthMax - (2.0 * M_PI)/config.wavelengthMin) / numWaveNumberSamples;
//            double waveNumber = c * dWaveNumber;

            // Gerstner Waves: https://people.computing.clemson.edu/~jtessen/reports/papers_files/coursenotes2004.pdf
            Vector2d profileXZ = waveDirection * profile.x();
            Vector3d profilePos = bilerp(posToIdxSpace(pos), b, waveNumber) * Vector3d(profileXZ.x(), profile.y(), profileXZ.y());
            Vector2d XZScalar = -waveDirection / waveNumber;
            profilePos = profilePos.cwiseProduct(Vector3d(XZScalar.x(), 1, XZScalar.y()));
            total += profilePos; // no shot this works first time. check here when things inevitably break
        }
    }

    return total;
}

void Amplitude::timeStep(double dt) {
    m_time += dt;
    advectionStep(dt);
    precomputeProfileBuffers(m_time);
//    print("sim time elapsed: " + std::to_string(m_time));
}
