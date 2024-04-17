#include "simulation.h"
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
    std::vector<Vector3f> vertices = this->m_grid.getVertices();
    std::vector<Vector3i> triangles = this->m_grid.getTriangles();



    m_shape.init(vertices, triangles);


    MatrixX3f all_vertices = MatrixX3f(vertices.size(), 3);
    int i = 0;
    for (unsigned long i = 0; i < vertices.size(); ++i) {
        all_vertices.row(i) = vertices[i];
    }
    coeffMin = all_vertices.colwise().minCoeff();
    coeffMax = all_vertices.colwise().maxCoeff();
}
