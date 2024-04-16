#include "simulation.h"
#include "graphics/meshloader.h"

#include <iostream>
#include <set>
#include <map>
#include <vector>

using namespace std;
using namespace Eigen;

Simulation::Simulation() {}

void Simulation::init(Eigen::Vector3f &coeffMin, Eigen::Vector3f &coeffMax)
{
    std::vector<Vector3f> vertices;
    std::vector<Vector3i> triangles;
    int size = 128;

    for(int i = 0; i < size; ++i){
        for(int j = 0; j < size; ++j){
            float period = 8 * M_PI * (i)/100.f;
            vertices.push_back(Vector3f(j, 10 * sin(period), i));
        }
    }

    for(int i = 0; i < size - 1; ++i){
        for(int j = 0; j < size - 1; ++j){
            int c1 = i + j * size;
            int c2 = i + 1 + j * size;
            int c3 = i + (j + 1) * size;
            int c4 = i + 1 + (j + 1) * size;

            triangles.push_back(Vector3i(c3, c4, c1));
            triangles.push_back(Vector3i(c4, c2, c1));
        }

    }





    m_shape.init(vertices, triangles);


    MatrixX3f all_vertices = MatrixX3f(vertices.size(), 3);
    int i = 0;
    for (unsigned long i = 0; i < vertices.size(); ++i) {
        all_vertices.row(i) = vertices[i];
    }
    coeffMin = all_vertices.colwise().minCoeff();
    coeffMax = all_vertices.colwise().maxCoeff();
}

// Move an anchored vertex, defined by its index, to targetPosition
void Simulation::move(int vertex, Vector3f targetPosition)
{
    std::vector<Eigen::Vector3f> new_vertices = m_shape.getVertices();
    const std::unordered_set<int>& anchors = m_shape.getAnchors();

    // TODO: implement ARAP here
    new_vertices[vertex] = targetPosition;

    // Here are some helpful controls for the application
    //
    // - You start in first-person camera mode
    //   - WASD to move, left-click and drag to rotate
    //   - R and F to move vertically up and down
    //
    // - C to change to orbit camera mode
    //
    // - Right-click (and, optionally, drag) to anchor/unanchor points
    //   - Left-click an anchored point to move it around
    //
    // - Minus and equal keys (click repeatedly) to change the size of the vertices

    m_shape.setVertices(new_vertices);
}
