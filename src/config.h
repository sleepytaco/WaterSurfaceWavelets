#ifndef CONFIG_H
#define CONFIG_H

#include <iostream>
#include <math.h>

#pragma once

void print(auto&& a) {std::cout << a << std::endl;}

struct Config {

    // ------ Amplitude.h ----------
    // number of samples for integrating water height
    const int numWaveNumberSamples = 1;
    const int numThetaSamples = 16;

    // number of samples for 4D amplitude grid
    const int dimXY = 128; // we assume same X and Y samples
    const int dimTheta = 16;
    const int dimK = 1;

    // simulation domain range
    // sample min max simulation ranges from supplemental paper - assuming meters as units
    const double xMin=-25; const double xMax=25;
    const double yMin=-25; const double yMax=25;
    const double thetaMin=0; const double thetaMax=2*M_PI;
    const double kMin=2.0*M_PI/0.02; const double kMax=2.0*M_PI/13.0; // wavenumber (k) = 2pi / wavelength (lamdba)

    // calculate the spatial res, theta res, wavenumber res based on the sim domain ranges
    const double dXY = (xMax - xMin) / dimXY; // the "width" of each cell XY grid cell basedon the sim domain range
    const double dTheta = (thetaMax - thetaMin) / dimTheta;
    const double dK = (kMin - kMax) / 2;
    // --------------------------------


    // ------- ProfileBuffer.h --------
    constexpr static const int bufferSize = 256;
    constexpr static const double wavelengthMin = 0.02; // defined
    constexpr static const double wavelengthMax = 13.0;
    constexpr static const double waveNumberMin = 314.159; // wave number = 2pi/wavelength
    constexpr static const double waveNumberMax = 0.483;
    constexpr static const double bufferExtent = 40.0;
    constexpr static const double pScalar = 2.5;
    constexpr static const int numSubintervals = 50; // For integration
    constexpr static const double g = 9.81;
    constexpr static const double sigma = 	0.072;
    constexpr static const double U = 1.0; // average wind speed
    // ------------------------------


    // -------- Solid-fluid Coupling --------
    const float objUniformScale = 5; // scaling factor to scale the loaded meshes by
    const float objMass = 1;
    const float fluidDensity = 997; // units are kg/m^3 --- this is the density of water
    // --------------------------------------

};

// The global Config object
extern Config config;
#endif // CONFIG_H


