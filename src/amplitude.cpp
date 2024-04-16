#include "amplitude.h"

Amplitude::Amplitude(int xSamples, int thetaSamples, int kSamples) {
    amplitudeGrid = vector<vector<vector<vector<float>>>>(xSamples,
                        vector<vector<vector<float>>>(xSamples,
                            vector<vector<float>>(thetaSamples,
                                vector<float>(kSamples, 0))));
}
