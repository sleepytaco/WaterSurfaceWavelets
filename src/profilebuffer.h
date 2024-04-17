#ifndef PROFILEBUFFER_H
#define PROFILEBUFFER_H

#include <functional>
#include <array>

#define FUNCTION std::function<double(double)>

class ProfileBuffer
{
public:
    ProfileBuffer();

    void precompute(double time);
    double getValueAt(double p);
    void integrationTest(); // temporary
    double dispersion(double waveNumber);
    double spectrum(double omega);

private:
    // Returns the integrand, a function that just takes k (wavelength)
    FUNCTION integrand(double p, double time);

    // Simpson's 1/3 rule: https://en.wikipedia.org/wiki/Simpson%27s_rule
    double integrate(double minBound, double maxBound, FUNCTION& fun);

    // See implementation details for explanation
    double h00(double s);
    double h01(double s);

    // Pierson-Moskowitz, eq. 20 from https://dl.acm.org/doi/pdf/10.1145/2791261.2791267

    // Parameters TODO: move all of these settings to a config file
    constexpr static const int bufferSize = 4096;
    constexpr static const double wavelengthMin = 0.02;
    constexpr static const double wavelengthMax = 13.0;
    constexpr static const double waveNumberMin = 314.159; // wave number = 2pi/wavelength
    constexpr static const double waveNumberMax = 0.483;
    constexpr static const double bufferExtent = 100.0;
    constexpr static const int numSubintervals = 50; // For integration
    constexpr static const double g = 9.81;
    constexpr static const double U = 10.0; // average wind speed

    // Members
    std::array<double, bufferSize> m_buffer;
};

#endif // PROFILEBUFFER_H
