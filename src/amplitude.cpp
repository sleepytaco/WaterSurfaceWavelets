#include "amplitude.h"

Amplitude::Amplitude(int xSamples, int thetaSamples, int kSamples) {
    amplitudeGrid.resize(this->w * this->h * this->k * this->theta);
}




float Amplitude::getAmplitudeVal(Vector2i a, int b, int c){
    return this->amplitudeGrid[this->gridIndex(a(0), a(1), b, c)];
}
