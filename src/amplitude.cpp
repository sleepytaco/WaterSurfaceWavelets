#include "amplitude.h"
#include <random>


Amplitude::Amplitude() {
    std::cout << "amplitude constructor" << std::endl;

    m_currentAmplitude = Grid();
    m_newAmplitude = Grid();
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
    return 0.5 * sqrt(config.g / waveNumber);
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

                double A = spatialDiffusion(dt, x_a, i, j, theta, k_c);

//                std::cout << A << std::endl;
                m_newAmplitude.set(i, j, theta, 0, A);

            }
        }
    }
    std::swap(m_newAmplitude, m_currentAmplitude);
    m_currentAmplitude.applyWind();
}



void Amplitude::precomputeProfileBuffers(double time) {
    m_profileBuffer.precompute(time);
}

std::tuple<Vector3d, Vector3d> Amplitude::waterHeight(Vector2d pos) {
    Vector3d total = Vector3d(0, 0, 0);
    double totalHeight = 0;
    Vector3d totalXDeriv = Vector3d::Zero();
    Vector3d totalYDeriv = Vector3d::Zero();

    for (int b = 0; b < numThetaSamples; b++) {
        double theta = 2.0 * M_PI * (double)b / (double)numThetaSamples;
        Vector2d waveDirection = Vector2d(cos(theta), sin(theta));
        double p = waveDirection.dot(pos) * config.pScalar;
        Vector4d profile = m_profileBuffer.getValueAt(p);

        for (int c = 1; c <= numWaveNumberSamples; c++) {

            double fraction = (double)c / (double)numWaveNumberSamples;
            double wavelength = config.wavelengthMax * fraction + config.wavelengthMin * (1 - fraction); // not 100% sure on this
            double waveNumber = 2.0 * M_PI / wavelength;

//            double dWaveNumber = ((2.0 * M_PI)/config.wavelengthMax - (2.0 * M_PI)/config.wavelengthMin) / numWaveNumberSamples;
//            double waveNumber = c * dWaveNumber;

            // Gerstner Waves: https://people.computing.clemson.edu/~jtessen/reports/papers_files/coursenotes2004.pdf
            Vector2d profileXZ = waveDirection * profile.x();
            double amplitude = bilerp(posToIdxSpace(pos), b, waveNumber);
            Vector3d profilePos = amplitude * Vector3d(waveDirection.x() * profile.x(), profile.y(), waveDirection.y() * profile.z());
            Vector2d XZScalar = -waveDirection / waveNumber;
            profilePos = profilePos.cwiseProduct(Vector3d(XZScalar.x(), 1, XZScalar.y()));
            total += profilePos; // no shot this works first time. check here when things inevitably break
            totalXDeriv += amplitude * waveDirection.x() * Vector3d(profile.z(), profile.w(), 0);
            totalYDeriv += amplitude * waveDirection.y() * Vector3d(0, profile.w(), profile.z());
        }
    }

    Vector3d normal = totalXDeriv.cross(totalYDeriv);
    if (normal == Vector3d::Zero()) {
        normal = Vector3d(0, 1, 0);
    } else {
//        std::cout << pos.x() << "," << pos.y() << " => " << normal.x() << "," << normal.y() << " " << normal.z() << std::endl;
    }
    normal = normal.normalized();
    if (normal.y() < -0.001) normal *= -1;
    if (normal.y() < 0.001) normal = Vector3d(0, 1, 0);
    return {total, normal};
}

void Amplitude::timeStep(double dt) {
    m_time += dt;
    m_timeAcc += dt;
    advectionStep(dt);
    precomputeProfileBuffers(m_time);

    static std::uniform_real_distribution<double> changeTime(2,5);
    static std::default_random_engine re;
    if (m_timeAcc > m_timeToChangeWind) {
        m_timeAcc = 0;
        m_timeToChangeWind = changeTime(re);
        m_currentAmplitude.changeWind();
    }
//    print("sim time elapsed: " + std::to_string(m_time));
}
