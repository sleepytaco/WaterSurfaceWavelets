#include "grid.h"

Grid::Grid()
{
    size_t gridSize = dims * dims * k * theta;
    amplitudeGrid.resize(gridSize);
    gridStep.resize(gridSize);

    for(int i = 0; i < mesh_dims; ++i){
        for(int j = 0; j < mesh_dims; ++j){
            float period = 8 * M_PI * (i)/100.f;
            vertices.push_back(Vector3f(j, 10 * sin(period), i));
        }
    }

    for(int i = 0; i < mesh_dims - 1; ++i){
        for(int j = 0; j < mesh_dims- 1; ++j){
            int c1 = i + j * mesh_dims;
            int c2 = i + 1 + j * mesh_dims;
            int c3 = i + (j + 1) * mesh_dims;
            int c4 = i + 1 + (j + 1) * mesh_dims;

            triangles.push_back(Vector3i(c3, c4, c1));
            triangles.push_back(Vector3i(c4, c2, c1));
        }

    }
}


double& Grid::operator()(Vector2i pos, int theta, int k) {
    return amplitudeGrid[gridIndex(pos(0), pos(1), theta, k)];
}

double& Grid::operator()(int x, int y, int theta, int k) {
    return amplitudeGrid[gridIndex(x,  y, theta, k)];
}
