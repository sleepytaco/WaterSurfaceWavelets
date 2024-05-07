#ifndef PROFILEBUFFER_H
#define PROFILEBUFFER_H

#include "config.h"
#include <functional>
#include <array>
#include "Eigen/Dense"

using namespace Eigen;

#define FUNCTION std::function<Vector4d(double)>

class ProfileBuffer
{
public:
    ProfileBuffer();

    void precompute(double time);
    Vector4d getValueAt(double p);
    double dispersion(double waveNumber);
    double spectrum(double waveNumber);

private:
    // Returns the integrand, a function that just takes k (wavelength)
    FUNCTION integrand(double p, double time);

    // Simpson's 1/3 rule: https://en.wikipedia.org/wiki/Simpson%27s_rule
    Vector4d integrate(double minBound, double maxBound, FUNCTION& fun);

    // Cubic polynomial weighting for seamless tiling (implementation details)
    double h00(double s);
    double h01(double s);
    Vector4d cubicPolynomial(double p, double waveNumber, double time);

    // Pierson-Moskowitz, eq. 20 from https://dl.acm.org/doi/pdf/10.1145/2791261.2791267

    // Parameters TODO: move all of these settings to a config file
    constexpr static const int bufferSize = config.bufferSize;
    constexpr static const double wavelengthMin = config.wavelengthMin;
    constexpr static const double wavelengthMax = config.wavelengthMax;
    constexpr static const double waveNumberMin = config.waveNumberMin; // wave number = 2pi/wavelength
    constexpr static const double waveNumberMax = config.waveNumberMax;
    constexpr static const double bufferExtent = config.bufferExtent;
    constexpr static const int numSubintervals = config.numSubintervals; // For integration
    constexpr static const double g = config.g;
    constexpr static const double U = config.U; // average wind speed

    // Members
    std::array<Vector4d, bufferSize> m_buffer;
};

#endif // PROFILEBUFFER_H
