//#include "levelset.h"
//#include "amplitude.h"
//LevelSet::LevelSet(LevelSetShape shape) : shape(shape) {};

//void LevelSet::boundaryReflection(Vector2d& advPos, int& thetaIdx) {
//    if ((shape = LevelSetShape::Circle)) {
//        circleReflection(advPos, thetaIdx);
//    } else if (shape == LevelSetShape::Square) {
//        squareReflection(advPos, thetaIdx);
//    }
//}

//bool Amplitude::boundaryReflection(double dt, int i, int j, int& thetaMin, int& thetaMax) {

//    int thetaIdx = 0;
//    double highestA = m_currentAmplitude.get(i, j, thetaIdx, 0);
//    // this gets the direction in which the wave at i, j is predominantly moving in
//    for (int th=1; th<config.dimTheta; ++th) {
//        double currA = m_currentAmplitude.get(i, j, th, 0);
//        if (currA > highestA) {
//            highestA = currA;
//            thetaIdx = th;
//        }
//    }

//    // advect the i, j pos forward in the wave direction corresponding to highestA
//    Vector2d x_a = idxToPos(i, j); // x_a = (x, y)
//    double k_c = dK; // wave number
//    double theta_d = thetaIdx*dTheta;
//    Vector2d waveDirection = Vector2d(cos(theta_d), sin(theta_d));
//    Vector2d advPos = x_a + dt * advectionSpeed(k_c) * waveDirection; // pos "FORWARD in time" ---- x_jump, y_jump

//    Vector2d boundaryNormal = Vector2d(0, 0);
//    double boundaryDistance = 0;
//    if (advPos.x() < xMin + dXY) {
//        boundaryNormal = Vector2d(1, 0);
//        boundaryDistance = abs(advPos.x() - xMin);
//    }
//    else if (advPos.x() > xMax - dXY) {
//        boundaryNormal = Vector2d(-1, 0);
//        boundaryDistance = abs(advPos.x() - xMax);
//    }
//    else if (advPos.y() < yMin + dXY) {
//        boundaryNormal = Vector2d(0, 1);
//        boundaryDistance = abs(advPos.y() - yMin);
//    }
//    else if (advPos.y() > yMax - dXY) {
//        boundaryNormal = Vector2d(0, -1);
//        boundaryDistance = abs(advPos.y() - yMax);
//    }
//    else {
//        return false;
//    }


//    waveDirection.norm();
//    Vector2d reflectedWaveDirection = waveDirection - 2 * (boundaryNormal.dot(waveDirection)) * boundaryNormal;
//    reflectedWaveDirection = boundaryNormal;

//    // convert reflectedWaveDirection into a discrete thetaIdx
//    double thetaX = std::atan2(reflectedWaveDirection.x(), reflectedWaveDirection.y()); // angle with x-axis
//    if (thetaX < 0) thetaX += 2 * M_PI; // ensure angle is in [0, 2*pi)
//    int thetaIndex = static_cast<int>(thetaX / (2 * M_PI / config.dimTheta)); // convert angle to bin index

//    // print("old theta idx: " + std::to_string(thetaIdx) + " | new theta idx: " + std::to_string(thetaIndex) );

//    // "spreading" the amplitude over a range around the thetaIdx (not exactly right...? as i dont do any modulo and clamp to 0/16 at boundaries)
//    int thetaRange = 4;
//    thetaMin = std::max(thetaIndex - thetaRange, 0);
//    thetaMax = std::min(thetaIndex + thetaRange, config.dimTheta-1);

////    for (int theta=0; theta<config.dimTheta; ++theta) { // b

////        if (theta>= minTheta && theta<=maxTheta) {
////            double newAmp = -30;
////            m_newAmplitude.set(i, j, theta, 0, newAmp);
////            for (int idx=1; idx<10; ++idx) {
////                m_newAmplitude.set(i+idx, j+idx, theta, 0, newAmp);
////                m_newAmplitude.set(i+idx, j, theta, 0, newAmp);
////                m_newAmplitude.set(i, j+idx, theta, 0, newAmp);
////            }
////        } else {
////            m_newAmplitude.set(i, j, theta, 0, 0);
////        }
////    }

//    return true;
//}
