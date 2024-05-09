#include "amplitude.h"

double Amplitude::spatialDiffusion(double dt, Vector2d idxPos, int xIdx, int yIdx, int thetaIdx, double waveNumber){
    /*
     *     Vector2d p_0 = idxPos;
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
     */

    double delta = 0.00001 * pow(dXY, 2) * pow(dK, 2) * fabs(advectionAccel(dK));
    double theta_d = thetaIdx*dTheta;

    Vector2d advPos = advectionPos(idxPos, dt, theta_d, dK);
    int theta_refl = thetaIdx;


    double ref_d = theta_refl * dTheta;

    Vector2d d = Vector2d(cos(ref_d), sin(ref_d));

    Vector2d p_n1 = idxPos + dXY * d;
    Vector2d p_n2 = idxPos + dXY * d * 2;
    Vector2d p_p1 = idxPos - dXY * d;
    Vector2d p_p2 = idxPos - dXY * d * 2;
    Vector2d p_p3 = idxPos - dXY * d * 3;

    Vector2d idxAdvpos = posToIdxSpace(advPos);

    std::vector<Vector2d> p{p_n2, p_n1, idxPos, p_p1, p_p2, p_p3};
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


    double adv_t = (advPos - idxPos).norm()/dXY;
    double spatial_diffuse = catmullRom(nv, adv_t);
    double angle_diffuse = diffusionStep(dt, spatial_diffuse, xIdx, yIdx, thetaIdx, waveNumber);
    return angle_diffuse;
}

double Amplitude::diffusionStep(double dt, double A, int xIdx, int yIdx, int thetaIdx, double waveNumber) {
    // fill in the newAmplitudeGrid with interpolated amplitude values at (x_jump, y_jump)
    double dissapation = 0; //2*1e-6*k_c*k_c*A; // including this term from paper stabalizes the water to a sheet pretty quickly
    double gamma = 0.025 * advectionSpeed(waveNumber) * pow(dTheta, 2)/dXY;
    double prevThetaIdx = thetaIdx == 0 ? 15 : thetaIdx - 1;
    double nextThetaIdx = thetaIdx == 15 ? 0 : thetaIdx + 1;
    double prev = m_currentAmplitude.get(xIdx, yIdx, prevThetaIdx, 0);
    double next = m_currentAmplitude.get(xIdx, yIdx, nextThetaIdx, 0);
    double s_d = (1 - 2 * gamma * dt/pow(dTheta, 2)) * A + gamma * dt/pow(dTheta, 2) * (prev + next);
    return s_d - dissapation;
}
