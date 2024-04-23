#ifndef GRID_H
#define GRID_H

#include "config.h"
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
    size_t dims = config.bufferSize;
    size_t theta = config.dimTheta;
    size_t k = config.dimK;
//    size_t mesh_dims = 128;

    std::vector<double> amplitudeGrid;

    std::vector<Vector3f> vertices;
    std::vector<Vector3i> triangles;

    int gridIndex(int i1, int i2, int i3, int i4){return i4 + k * (i3 + theta * (i2 + dims * i1));}
};

#endif // GRID_H
