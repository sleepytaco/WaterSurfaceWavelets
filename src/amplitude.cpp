#include "amplitude.h"
#include <random>

void print(auto&& a) {std::cout << a << std::endl;}

Amplitude::Amplitude() {
    std::cout << "amplitude constructor" << std::endl;
    double lower_bound = -2;
    double upper_bound = 2;
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

//    if (idxSpaceX > dimXY) idxSpaceX -= dimXY;
//    if (idxSpaceY > dimXY) idxSpaceY -= dimXY;
//    if (idxSpaceX < 0) idxSpaceX += dimXY;
//    if (idxSpaceY < 0) idxSpaceY += dimXY;
    if (idxSpaceX > dimXY) idxSpaceX = dimXY;
    if (idxSpaceY > dimXY) idxSpaceY = dimXY;
    if (idxSpaceX < 0) idxSpaceX = 0;
    if (idxSpaceY < 0) idxSpaceY = 0;

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
    if (advPos.x() < xMin) {
        boundaryNormal = Vector2d(1, 0);
        boundaryDistance = abs(advPos.x() - xMin);
    }
    else if (advPos.x() > xMax) {
        boundaryNormal = Vector2d(-1, 0);
        boundaryDistance = abs(advPos.x() - xMax);
    }
    else if (advPos.y() < yMin) {
        boundaryNormal = Vector2d(0, 1);
        boundaryDistance = abs(advPos.y() - yMin);
    }
    else if (advPos.y() > yMax) {
        boundaryNormal = Vector2d(0, -1);
        boundaryDistance = abs(advPos.y() - yMax);
    }
    else {
        return;
    }

    double theta = thetaIdx * dTheta;
    Vector2d waveDirection = Vector2d(cos(theta), sin(theta));

    Vector2d advPos_refl = advPos - 2 * boundaryDistance * boundaryNormal;
    advPos = advPos_refl;
    // assert(advPos.x() >= xMin && advPos.x() <= xMax && advPos.y() >= yMin && advPos.y() <= yMax);

    Vector2d k_refl = waveDirection - 2 * (boundaryNormal.dot(waveDirection)) * boundaryNormal;
    double theta_refl = atan2(k_refl.y(), k_refl.x());
    if (theta_refl < 0) theta_refl += 2*M_PI;
    int thetaIdx_refl = round(theta_refl / dTheta);
    thetaIdx = thetaIdx_refl;
}

// assuming deep water dispersion: this is the equation for omega' if omega = sqrt(gk)
// TODO: this might be too simple and need to factor in additional terms
double Amplitude::advectionSpeed(double waveNumber) {
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
    return advPos;
}

// advectionStep moves the current amplitudeGrid forward in time
void Amplitude::advectionStep(double dt) {

    #pragma omp parallel for collapse(2)
    for (int i=0; i<dimXY; ++i) { // a
        for (int j=0; j<dimXY; ++j) { // a
            for (int theta=0; theta<dimTheta; ++theta) { // b
                Vector2d x_a = idxToPos(i, j); // x_a = (x, y)
                double k_c = dK; // wave number
                double A = spacialDiffusion(dt, x_a, i, j, theta, k_c);
                m_newAmplitude.get(i, j, theta, 0) = A;
            }
        }
    }
    std::swap(m_newAmplitude, m_currentAmplitude);
}

double Amplitude::spacialDiffusion(double dt, Vector2d idxPos, int xIdx, int yIdx, int thetaIdx, double waveNumber){
    Vector2d p_0 = idxPos;
    double delta = 0.00001 * pow(dXY, 2) * pow(dK, 2) * fabs(advectionAccel(dK));

    double theta_d = thetaIdx*dTheta;

    Vector2d advPos = advectionPos(idxPos, dt, theta_d, dK);
    int theta_refl = thetaIdx;
    //boundaryReflection(advPos, thetaIdx);

    double ref_d = theta_refl * dTheta;

    Vector2d idxAdvpos = posToIdxSpace(advPos);
    Vector2d d = Vector2d(cos(ref_d), sin(ref_d));

    Vector2d p_n1 = idxPos - dXY * d;
    Vector2d p_n2 = idxPos - dXY * d * 2;
    Vector2d p_p1 = idxPos + dXY * d;
    Vector2d p_p2 = idxPos + dXY * d * 2;
    Vector2d p_p3 = idxPos + dXY * d * 3;

    std::vector<Vector2d> p{p_n2, p_n1, p_0, p_p1, p_p2, p_p3};
    std::vector<double> v;
    std::vector<double> nv;
    #pragma omp parallel
    for(int i = 0; i < p.size(); ++i){
        v.push_back(bilerp(posToIdxSpace(p[i]), theta_refl, waveNumber));
    }

    for(int i = 1; i < v.size() - 1; ++i){
        double diffused = (1 - 2 * delta * dt/pow(dXY, 2)) * v[i] + delta * dt/pow(dXY, 2) * (v[i - 1] + v[i + 1]);
        nv.push_back(diffused);
    }

    double adv_t = (advPos - p_0).norm()/dXY;
    double spatial_diffuse = catmullRom(nv, adv_t);
    return spatial_diffuse;
}

double Amplitude::diffusionStep(double dt, Vector2d idxPos, int xIdx, int yIdx, int thetaIdx, double waveNumber) {
    // fill in the newAmplitudeGrid with interpolated amplitude values at (x_jump, y_jump)
    double gamma = 0.025 * advectionSpeed(waveNumber) * pow(dTheta, 2)/dimXY;


    double A = bilerp(idxPos, thetaIdx, waveNumber);

    //calulate diffusion here
    double dissapation = 0; //2*1e-6*k_c*k_c*A; // including this term from paper stabalizes the water to a sheet pretty quickly

    double prev = m_currentAmplitude.get(xIdx, yIdx, (thetaIdx - 1)%dimTheta, 0);
    double next = m_currentAmplitude.get(xIdx, yIdx, (thetaIdx + 1)%dimTheta, 0);
    double s_d = (1 - 2 * gamma * dt/pow(dTheta, 2)) * A + gamma * dt/pow(dTheta, 2) * (prev + next);
    return s_d - dissapation;
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
//    std::cout << "sim time elapsed: " << m_time << std::endl;
}
