#include "amplitude.h"
#include <random>

void print(auto&& a) {std::cout << a << std::endl;}

Amplitude::Amplitude() {
    std::cout << "amplitude constructor" << std::endl;
    double lower_bound = -30;
    double upper_bound = 30;
    std::uniform_real_distribution<double> unif(lower_bound,upper_bound);
    std::default_random_engine re;

    m_currentAmplitude = Grid();
    m_newAmplitude = Grid();
    for (int i=0; i<dimXY; ++i) { // a
        for (int j=0; j<dimXY; ++j) { // a
            for (int theta=0; theta<dimTheta; ++theta) { // b
                // uncomment to init a sqaure with 0 amplitude in the center of the grid
//                Vector2d x_a = idxToPos(i, j); // x_a = (x, y)
//                if ((x_a.x() >= 500 && x_a.x() <= 3500) && (x_a.y() >= 500 && x_a.y() <= 3500)) {
//                    //m_currentAmplitude.get(i, j, theta, 0) = unif(re); 20 * sin((i + j) / 1);
//                    continue;
//                }

                m_currentAmplitude.get(i, j, theta, 0) = unif(re) * sin((i + j) / 2);
            }
        }
    }
    m_profileBuffer = ProfileBuffer();
}

Vector2d Amplitude::idxToPos(int i, int j) {
    double x = xMin + i * dXY;
    double y = yMin + j * dXY;
    return Vector2d(x, y);
}

Vector2d Amplitude::posToIdxSpace(Vector2d pos) {
    double idxSpaceX = (pos.x() - xMin) / dXY;
    double idxSpaceY = (pos.y() - yMin) / dXY;

    if (idxSpaceX > dimXY) idxSpaceX -= dimXY;
    if (idxSpaceY > dimXY) idxSpaceY -= dimXY;
    if (idxSpaceX < 0) idxSpaceX += dimXY;
    if (idxSpaceY < 0) idxSpaceY += dimXY;

    return Vector2d(idxSpaceX, idxSpaceY);
}

// assuming deep water dispersion: this is the equation for omega' if omega = sqrt(gk)
// TODO: this might be too simple and need to factor in additional terms
double Amplitude::advectionSpeed(double waveNumber) {
    //return sqrt(9.8) * 0.5 / sqrt(waveNumber);
    double val = config.g/(2 * sqrt(config.g * waveNumber + config.sigma * pow(waveNumber, 3)));
    return val;
}

double Amplitude::advectionAccel(double waveNumber) {
    //return sqrt(9.8) * 0.5 / sqrt(waveNumber);
    double val =-pow(config.g, 2)/(4 * pow(config.g * waveNumber + config.sigma * pow(waveNumber, 3), 1.5));
    return val;
}

Vector2d Amplitude::advectionPos(Vector2d pos, double dt, double theta, double waveNumber) {
    Vector2d waveDirection = Vector2d(cos(theta), sin(theta));    
    Vector2d advPos = pos - dt * advectionSpeed(waveNumber) * waveDirection; // pos "back in time" ---- x_jump, y_jump

    // boundary condition check - if adv pos leaves the simulation domain
    if (advPos.x() < xMin || advPos.x() > xMax-dXY || advPos.y() < yMin || advPos.y() > yMax-dXY) {

        //print(advPos.transpose());
        Vector2d boundaryNormal(0, 0); // figure out the boundary normal based on dist to grid boundaries

        // assuming the grid is made by lines x=0, x=4000, y=0, y=4000
        double topDist = abs(advPos.y() - yMin);
        double rightDist = abs(advPos.x() - xMax);
        double bottomDist = abs(advPos.y() - yMax);
        double leftDist = abs(advPos.x() - xMin);
        double minDist = std::min({leftDist, rightDist, topDist, bottomDist});

        if (minDist == topDist) boundaryNormal = Vector2d(0, -1); // top boundary - normal points down
        else if (minDist == rightDist) boundaryNormal = Vector2d(-1, 0); // right boundary - normal points left
        else if (minDist == bottomDist) boundaryNormal = Vector2d(0, 1); // bottom boundary - normal points up
        else if (minDist == leftDist) boundaryNormal = Vector2d(1, 0); // left boundary - normal points right

        // if (boundaryNormal.dot(waveVector) < 0) waveVector = -waveVector;
        Vector2d k_refl = waveDirection - 2 * (boundaryNormal.dot(waveDirection)) * boundaryNormal;

        Vector2d pos_new = advPos + k_refl;
        k_refl.normalize();
        advPos = pos_new - dt * advectionSpeed(waveNumber) * k_refl;

        // assert that the new advection pos is within the sim domain --- it fails tho...
        // assert(advPos.x() >= xMin && advPos.x() <= xMax && advPos.y() >= yMin && advPos.y() <= yMax);
    }

    return advPos;
}

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

// advectionStep moves the current amplitudeGrid forward in time
void Amplitude::advectionStep(double dt) {

    #pragma omp parallel for collapse(2)
    for (int i=0; i<dimXY; ++i) { // a
        for (int j=0; j<dimXY; ++j) { // a
            for (int theta=0; theta<dimTheta; ++theta) { // b
                Vector2d x_a = idxToPos(i, j); // x_a = (x, y)
                double theta_b = theta*dTheta;
                double k_c = dK; // wave number

                Vector2d advPos = advectionPos(x_a, dt, theta_b, k_c);// pos "back in time" ---- x_jump, y_jump
                Vector2d idxSpaceAdvPos = posToIdxSpace(advPos);

                // fill in the newAmplitudeGrid with interpolated amplitude values at (x_jump, y_jump)
                double A = interpolateAmplitude(idxSpaceAdvPos, theta);
                double dissapation = 0; //2*1e-6*k_c*k_c*A; // including this term from paper stabalizes the water to a sheet pretty quickly
                m_newAmplitude.get(i, j, theta, 0) = A - dissapation;
            }
        }
    }

    std::swap(m_newAmplitude, m_currentAmplitude);
}

void Amplitude::precomputeProfileBuffers(double time) {
    m_profileBuffer.precompute(time);
}

Vector3d Amplitude::waterHeight(Vector2d pos) {
    Vector3d total = Vector3d(0, 0, 0);
    double totalHeight = 0;

    for (int b = 1; b <= numThetaSamples; b++) {
        double theta = 2.0 * M_PI * (double)b / (double)numThetaSamples;
        Vector2d waveDirection = Vector2d(cos(theta), sin(theta));
        double p = waveDirection.dot(pos);
        Vector2d profile = m_profileBuffer.getValueAt(p);

        for (int c = 1; c <= numWaveNumberSamples; c++) {
            double fraction = (double)c / (double)numWaveNumberSamples;
            double wavelength = config.wavelengthMax * fraction + config.wavelengthMin * (1 - fraction); // not 100% sure on this
            double waveNumber = 2.0 * M_PI / wavelength;

            // Gerstner Waves: https://people.computing.clemson.edu/~jtessen/reports/papers_files/coursenotes2004.pdf
            Vector2d profileXZ = waveDirection * profile.x();
            Vector3d profilePos = interpolateAmplitude(pos, theta) * Vector3d(profileXZ.x(), profile.y(), profileXZ.y());
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
    std::cout << "sim time elapsed: " << m_time << std::endl;
}
