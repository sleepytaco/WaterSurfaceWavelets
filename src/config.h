#ifndef CONFIG_H
#define CONFIG_H

#include <iostream>
#include <math.h>

#pragma once

struct Config {

    // ------ Amplitude.h ----------

    // number of samples for integrating water height
    const int numWaveNumberSamples = 1;
    const int numThetaSamples = 16;

    // number of samples for 4D amplitude grid
    const int dimXY = 256; // we assume same X and Y samples
    const int dimTheta = 16;
    const int dimK = 1;

    // simulation domain range
    const double xMin=0; const double xMax=4000;
    const double yMin=0; const double yMax=4000;
    const double thetaMin=0; const double thetaMax=2*M_PI;
    const double kMin=2.0*M_PI/0.02; const double kMax=2.0*M_PI/13.0; // wavenumber (k) = 2pi / wavelength (lamdba)

    const double dXY = (xMax - xMin) / dimXY; // the "width" of each cell XY grid cell basedon the sim domain range
    const double dTheta = (thetaMax - thetaMin) / dimTheta;
    const double dK = (kMin - kMax) / 2;


    //    constexpr static const double wavelengthMin = 0.02; // note: if these are changed so should the values in profilebuffer.h, these should be moved to a config file
    //    constexpr static const double wavelengthMax = 13.0;
    // --------------------------------


    // ------- ProfileBuffer.h --------
    constexpr static const int bufferSize = 4096;
    constexpr static const double wavelengthMin = 0.02; // defined
    constexpr static const double wavelengthMax = 13.0;
    constexpr static const double waveNumberMin = 314.159; // wave number = 2pi/wavelength
    constexpr static const double waveNumberMax = 0.483;
    constexpr static const double bufferExtent = 40.0;
    constexpr static const int numSubintervals = 50; // For integration
    constexpr static const double g = 9.81;
    constexpr static const double U = 10.0; // average wind speed
    // ------------------------------


    // ------- Grid.h -------------
//    size_t dims = 16;
//    size_t theta = 16;
//    size_t k = 1;
    // ----------------------------

};

// The global Config object
extern Config config;
#endif // CONFIG_H


