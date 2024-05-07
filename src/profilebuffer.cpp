#include "profilebuffer.h"
#include <iostream>
#include <cmath>
#include <cassert>

ProfileBuffer::ProfileBuffer()
{
    std::cout << "profile buffer constructor" << std::endl;
}

// The profiles this produces are very small in magnitude (like 0.001) but the visualization looks wave-like, maybe we just need to scale it?
void ProfileBuffer::precompute(double time) {
    for (int i = 0; i < bufferSize; ++i) {
        double p = (i * bufferExtent) / bufferSize;
        FUNCTION fun = integrand(p, time);
        Vector4d value = integrate(waveNumberMin, waveNumberMax, fun);
        m_buffer[i] = value;
    }
}

Vector4d ProfileBuffer::getValueAt(double p) {
    double interpolatePos = p * bufferSize / bufferExtent; // where we want to interpolate
    int prev = (int)floor(interpolatePos);
    int next = prev + 1;
    double remainder = interpolatePos - prev;
    // https://stackoverflow.com/questions/3417183/modulo-of-negative-numbers
    int prevWrapped = (prev % bufferSize + bufferSize) % bufferSize; // wrap for if p is outside of buffer
    int nextWrapped = (next % bufferSize + bufferSize) % bufferSize;
    assert(prevWrapped >= 0 && prevWrapped < bufferSize);
    assert(nextWrapped >= 0 && nextWrapped < bufferSize);
    return m_buffer[prevWrapped] * remainder + m_buffer[nextWrapped] * (1 - remainder);
}

// The authors say they use eq. 32, but I have no idea how because eq. 32 is dependent on wave direction
// but our integral is only over wavenumber. Instead, we're using eq. 20 because it only depends on wavenumber
double ProfileBuffer::spectrum(double waveNumber) {
//    double omega = dispersion(waveNumber);
//    const static double A = 0.0081 * pow(g, 2);
//    const static double B = 0.6858 * pow(g/U, 4);
//    return A / pow(omega, 5) * exp(-B / pow(omega, 4));
    double A = pow(waveNumber, 0.21); // original pow(2, 1.5*zeta)
    double B = exp(-1.8 * pow(waveNumber, 2) / pow(config.U, 4));
    return 0.14 * sqrt(A * B);
}

Vector4d ProfileBuffer::cubicPolynomial(double p, double waveNumber, double time) {
    double omega = dispersion(waveNumber);
    double phi = spectrum(omega);
    double waveLength = 2 * M_PI / waveNumber;
    double s = p / bufferExtent;
    double weight1 = h00(s);
    double weight2 = h01(s);
    double angle1 = waveNumber * p - omega * time;
    double angle2 = waveNumber * (p + bufferExtent) - omega * time;
    double angle3 = waveNumber * (p - bufferExtent) - omega * time;
    double sinAngle1 = sin(angle1);
    double cosAngle1 = cos(angle1);
    double sinAngle2 = sin(angle2);
    double cosAngle2 = cos(angle2);
    double sinAngle3 = sin(angle3);
    double cosAngle3 = cos(angle3);
    double horizontal = -(sinAngle1 + weight1 * sinAngle2 + weight2 * sinAngle3);
    double vertical = cosAngle1 + weight1 * cosAngle2 + weight2 * cosAngle3;
    double horizontalDeriv = -waveNumber * (cosAngle1 + weight1 * cosAngle2 + weight2 * cosAngle3);
    double verticalDeriv = -waveNumber * (sinAngle1 + weight1 * sinAngle2 + weight2 * sinAngle3);
    return 0.5 * phi * waveLength * Vector4d(horizontal, vertical, horizontalDeriv, verticalDeriv);
}

FUNCTION ProfileBuffer::integrand(double p, double time) {
    return [this, p, time](double waveNumber) -> Vector4d {
        return cubicPolynomial(p, waveNumber, time);
    };
}


// Simpson's 1/3 rule: https://en.wikipedia.org/wiki/Simpson%27s_rule
Vector4d ProfileBuffer::integrate(double minBound, double maxBound, FUNCTION& fun) {
    double delta = (maxBound - minBound) / numSubintervals;
    Vector4d integral = Vector4d::Zero();
    double x = minBound;

    for (int i = 0; i < numSubintervals; ++i) {
        integral += delta / 6.0 * (fun(x) + 4 * fun(x + delta / 2) + fun(x + delta));
        x += delta;
    }

    return integral;
}

double ProfileBuffer::dispersion(double waveNumber) {
    return sqrt(g * waveNumber); // deep water
}

double ProfileBuffer::h00(double s) {
    return 2.0 * pow(s, 3) - 3.0 * pow(s, 2) + 1.0;
}

double ProfileBuffer::h01(double s) {
    return -2.0 * pow(s, 3) + 3.0 * pow(s, 2);
}

