#include "simulation.h"
#include "graphics/meshloader.h"
#include "amplitude.h"

#include <iostream>
#include <set>
#include <map>
#include <vector>

using namespace std;
using namespace Eigen;

Simulation::Simulation()
    : m_amplitude(1024, 16, 1)
{
}

void Simulation::update(double deltaTime) {
    m_amplitude.timeStep(deltaTime / 2);
    setWaterHeights();
}

void Simulation::setWaterHeights() {
    std::vector<Eigen::Vector3f> new_vertices = m_shape.getVertices();
    for (int i = 0; i < new_vertices.size(); i++) {
        Vector3f vertex = new_vertices[i];
        Vector2d xz = Vector2d(vertex.x(), vertex.z());
        new_vertices[i].y() = m_amplitude.waterHeight(xz);
    }
    m_shape.setVertices(new_vertices);
}

void Simulation::init(Eigen::Vector3f &coeffMin, Eigen::Vector3f &coeffMax)
{
    std::vector<Vector3f> vertices;
    std::vector<Vector3i> triangles;
    int size = 64;

    for(int i = 0; i < size; ++i){
        for(int j = 0; j < size; ++j){
            float period = 8 * M_PI * (i)/100.f;
//            vertices.push_back(Vector3f(j, 10 * sin(period), i));
            vertices.push_back(Vector3f(j, 0, i));

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
