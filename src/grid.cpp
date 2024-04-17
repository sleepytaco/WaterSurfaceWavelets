#include "grid.h"

Grid::Grid()
{
    size_t gridSize = this->dims * this->dims * this->k * this->theta;
    this->amplitudeGrid.resize(gridSize);
    this->gridStep.resize(gridSize);

    for(int i = 0; i < this->mesh_dims; ++i){
        for(int j = 0; j < this->mesh_dims; ++j){
            float period = 8 * M_PI * (i)/100.f;
            this->vertices.push_back(Vector3f(j, 10 * sin(period), i));
        }
    }

    for(int i = 0; i < this->mesh_dims - 1; ++i){
        for(int j = 0; j < this->mesh_dims- 1; ++j){
            int c1 = i + j * this->mesh_dims;
            int c2 = i + 1 + j * this->mesh_dims;
            int c3 = i + (j + 1) * this->mesh_dims;
            int c4 = i + 1 + (j + 1) * this->mesh_dims;

            this->triangles.push_back(Vector3i(c3, c4, c1));
            this->triangles.push_back(Vector3i(c4, c2, c1));
        }

    }
}


float Grid::getAmplitudeVal(Vector2i a, int b, int c){
    return this->gridIndex(a(0), a(1), b, c);
}
