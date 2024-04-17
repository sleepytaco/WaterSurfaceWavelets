#include "grid.h"

Grid::Grid()
{
    std::cout << "grid constructor" << std::endl;

    size_t gridSize = dims * dims * k * theta;
    amplitudeGrid.resize(gridSize);
    std:fill(amplitudeGrid.begin(), amplitudeGrid.end(), 0);
}


double& Grid::get(Vector2i pos, int theta, int k) {
    return amplitudeGrid[gridIndex(pos(0), pos(1), theta, k)];
}

double& Grid::get(int x, int y, int theta, int k) {
    return amplitudeGrid[gridIndex(x,  y, theta, k)];
}
