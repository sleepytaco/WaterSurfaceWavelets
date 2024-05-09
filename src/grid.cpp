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
                this->set(i, j, th, 0, init(re));

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

void Grid::applyWind() {

    std::uniform_real_distribution<double> ampUpperBound(2, 5);
    std::uniform_real_distribution<double> thetaIdx(3, 10);
    std::uniform_real_distribution<double> ampRand(0, ampUpperBound(re));

    for (int i=0; i<=dims; ++i) { // a
        for (int j=0; j<=dims; ++j) { // a
            int thIdx = round(thetaIdx(re)); // just to get some slight randomness in the wind direction each time
            double currA = this->get(i, j, thIdx, 0);
            this->set(i, j, thIdx, 0, currA+ampRand(re));
        }
    }
}
