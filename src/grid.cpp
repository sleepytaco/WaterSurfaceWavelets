#include "grid.h"

Grid::Grid()
{
    std::cout << "grid constructor" << std::endl;

    size_t gridSize = dims * dims * k * theta;
    amplitudeGrid.resize(gridSize);
    std:fill(amplitudeGrid.begin(), amplitudeGrid.end(), 0);
    this->unif = std::uniform_real_distribution<double>(lower_bound, upper_bound);
    std::uniform_real_distribution<double> init(-30, 30);

    for (int i=0; i<=dims; ++i) { // a
        for (int j=0; j<=dims; ++j) { // a
            for (int th=0; th<=theta; ++th) { // b
                // uncomment to init a sqaure with 0 amplitude in the center of the grid
                //                Vector2d x_a = idxToPos(i, j); // x_a = (x, y)
                //                if ((x_a.x() >= 500 && x_a.x() <= 3500) && (x_a.y() >= 500 && x_a.y() <= 3500)) {
                //                    //m_currentAmplitude.get(i, j, theta, 0) = unif(re); 20 * sin((i + j) / 1);
                //                    continue;
                //                }
//                if ((i >= 50 && i <= 80 && j >= 50 && j <= 80))
                //                if (j >= 100)
                //m_currentAmplitude.get(i, j, theta, 0) = unif(re) * sin((i + j) / 2);
//                this->set(i, j, th, 0, init(re));

            }
        }
    }

    applyWind();
}


double& Grid::get(Vector2i pos, int theta, int k) {
    return amplitudeGrid[gridIndex(pos(0), pos(1), theta, k)];
}

double Grid::get(int x, int y, int theta, int k) {
    if(!containment(x, y, theta, k)){
        if(theta == 5)
        return 10;
        else
        return unif(re);


    }
//    x%=dims;
//    y%=dims;
    return amplitudeGrid[gridIndex(x,  y, theta, k)];
}

void Grid::set(int x, int y, int theta, int k, double A){
    if(x >= 0 && x < dims && y >= 0 && y <dims && theta>=0 && theta < this->theta)
    amplitudeGrid[gridIndex(x,  y, theta, k)] = A;
}

void Grid::changeWind() {
    static std::uniform_real_distribution<double> dirDelta(-3, 3);
    static std::uniform_real_distribution<double> intensityDelta(-2, 2);

    windDirection += dirDelta(re);
    if (windDirection < 0) windDirection += config.dimTheta;
    else if (windDirection > config.dimTheta) windDirection -= config.dimTheta;

    windIntensity += intensityDelta(re);
    if (windIntensity < 0) windIntensity += 7;
    else if (windIntensity > 7) windIntensity -= 7;
}

void Grid::applyWind() {

    static std::uniform_real_distribution<double> thetaDelta(-1, 1);
    static std::uniform_real_distribution<double> ampScale(0, 1);

    for (int i=0; i<=dims; ++i) { // a
        for (int j=0; j<=dims; ++j) { // a
        int thIdx = round(windDirection + thetaDelta(re)); // just to get some slight randomness in the wind direction each time
        if (thIdx < 0) thIdx += config.dimTheta;
        else if (thIdx > config.dimTheta) thIdx -= config.dimTheta;
            double currA = this->get(i, j, thIdx, 0);
            this->set(i, j, thIdx, 0, currA + windIntensity*ampScale(re));
        }
    }
}
