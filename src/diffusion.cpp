#include "amplitude.h"

double Amplitude::spatialDiffusion(double dt, Vector2d idxPos, int xIdx, int yIdx, int thetaIdx, double waveNumber){
    Vector2d p_0 = idxPos;

    double theta_d = thetaIdx*dTheta;

    Vector2d advPos = advectionPos(idxPos, dt, theta_d, dK);
    int theta_refl = thetaIdx;

    Vector2d idxAdvpos = posToIdxSpace(advPos);
    double spatial_diffuse = bilerp(idxAdvpos, theta_refl, dK);
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
