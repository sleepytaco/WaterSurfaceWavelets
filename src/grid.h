#ifndef GRID_H
#define GRID_H

#include "Eigen/Dense"
#include <iostream>
#include <unordered_set>
#include <deque>
#include <set>

using namespace Eigen;

class Grid
{
public:
    Grid();

    float getAmplitudeVal(Vector2i a, int b, int c); // gets the amplitude value stored in the grid

    void setAmplitude(int x, int y, int theta, int k, float val);

    void updateGrid(){this->amplitudeGrid = this->gridStep;}

    std::vector<Vector3f>& getVertices(){return this->vertices;}
    std::vector<Vector3i>& getTriangles(){return this->triangles;}


private:
    const size_t dims = 4096;
    const size_t theta = 16;
    const size_t k = 1;

    const size_t mesh_dims = 128;

    std::vector<float> amplitudeGrid;

    std::vector<float> gridStep;

    std::vector<Vector3f> vertices;
    std::vector<Vector3i> triangles;

    int gridIndex(int i1, int i2, int i3, int i4){return i1 + dims * (i2 + dims * (i3 * theta + i4));}

};

#endif // GRID_H
