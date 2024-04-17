#ifndef GRID_H
#define GRID_H

#include "Eigen/Dense"
#include <iostream>
#include <unordered_set>
#include <deque>
#include <set>
#include <vector>

using namespace Eigen;

class Grid
{
public:
    Grid();

    double& get(Vector2i pos, int theta, int k);
    double& get(int x, int y, int theta, int k);

    std::vector<Vector3f>& getVertices(){return vertices;}
    std::vector<Vector3i>& getTriangles(){return triangles;}


private:
    size_t dims = 256;
    size_t theta = 16;
    size_t k = 1;

    size_t mesh_dims = 128;

    std::vector<double> amplitudeGrid;

    std::vector<Vector3f> vertices;
    std::vector<Vector3i> triangles;

    int gridIndex(int i1, int i2, int i3, int i4){return i4 + k * (i3 + theta * (i2 + dims * i1));}
};

#endif // GRID_H
