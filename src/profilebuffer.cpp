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
        Vector2d value = integrate(waveNumberMin, waveNumberMax, fun);
        m_buffer[i] = value;
    }
}

Vector2d ProfileBuffer::getValueAt(double p) {
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
double ProfileBuffer::spectrum(double omega) {
    const static double A = 0.0081 * pow(g, 2);
    const static double B = 0.6858 * pow(g/U, 4);
    return A / pow(omega, 5) * exp(-B / pow(omega, 4));
}

FUNCTION ProfileBuffer::integrand(double p, double time) {
    return [this, p, time](double waveNumber) -> Vector2d {
        double omega = dispersion(waveNumber);
        double waveLength = 2 * M_PI / waveNumber;
        double s = p / bufferExtent;
        assert(s >= 0 && s <= 1);
        // TODO: the paper doesn't do a good job of differentiating between wavenumber and wavelength, we might have to fiddle with the below math
        double phi = spectrum(omega);
//        double phi = 1;
        double weight1 = h00(s);
        double weight2 = h01(s);
        double angle1 = waveNumber * p - omega * time;
        double angle2 = waveNumber * (p + bufferExtent) - omega * time;
        double angle3 = waveNumber * (p - bufferExtent) - omega * time;
        double part1X = -sin(angle1);
        double part2X = weight1 * -sin(angle2);
        double part3X = weight2 * -sin(angle3);
        double part1Y = cos(angle1);
        double part2Y = weight1 * cos(angle2);
        double part3Y = weight2 * cos(angle3);
//        std::cout << "evaluated at " << waveNumber << ", p = " << p << std::endl;
        return 0.5 * phi * waveLength * Vector2d(part1X + part2X + part3X, part1Y + part2Y + part3Y);
    };
}


// Simpson's 1/3 rule: https://en.wikipedia.org/wiki/Simpson%27s_rule
Vector2d ProfileBuffer::integrate(double minBound, double maxBound, FUNCTION& fun) {
    double delta = (maxBound - minBound) / numSubintervals;
    Vector2d integral = Vector2d(0, 0);
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

