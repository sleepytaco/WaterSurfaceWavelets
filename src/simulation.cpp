#include "simulation.h"
#include "amplitude.h"
#include "config.h"
#include <iostream>
#include <set>
#include <map>
#include <vector>
#include "graphics/meshloader.h"

using namespace std;
using namespace Eigen;

Simulation::Simulation()
    : m_amplitude()
{}

void Simulation::update(double deltaTime) {
    for (int i=0; i<m_particleSystems.size(); ++i) { // step each particle system forward
        solver->RK4(m_particleSystems[i], deltaTime / 2); // use RK4 to integrate the state of the particle system forward in time
        m_fallingShapes[i]->setVertices(m_particleSystems[i]->getVertices());
    }
    m_amplitude.timeStep(deltaTime);
    setWaterHeights();
}

void Simulation::setWaterHeights() {
//    std::vector<Eigen::Vector3f> new_vertices = m_shape.getVertices();
    for (int i = 0; i < undisturbedPoints.size(); i++) {
        Vector3f vertex = undisturbedPoints[i];
        Vector2d xz = Vector2d(vertex.x() / config.dimXY * (config.xMax - config.xMin) + config.xMin, vertex.z() / config.dimXY * (config.yMax - config.yMin) + config.yMin);
//        std::cout << "Getting height at " << xz.x() << "," << xz.y() << std::endl;
        Vector3f displacement = m_amplitude.waterHeight(xz).cast<float>();
        newPoints[i].x() = vertex.x() + displacement.x();
        newPoints[i].y() = displacement.z();
        newPoints[i].z() = vertex.z() + displacement.y();
    }
    m_shape.setVertices(newPoints);
}

void Simulation::init(Eigen::Vector3f &coeffMin, Eigen::Vector3f &coeffMax)
{
    // drop objects on the water surface mesh (calling these adds System() instances to m_particleSystems)
    initFallingParticleSystem("./meshes/cube.obj", Vector3f(config.dimXY/2, 30, config.dimXY/2));
    // initFallingParticleSystem("./meshes/cube.obj", Vector3f(config.dimXY/4, 30, config.dimXY/4));

    // give all particleSystems access to the amplitude function (for solid-fluid coupling)
    for (int i=0; i<m_particleSystems.size(); ++i) {
        m_particleSystems[i]->setAmplitudeFunction(&m_amplitude);
    }

    std::vector<Vector3f> vertices;
    std::vector<Vector3i> triangles;
    int size = config.dimXY;

    for(int i = 0; i < size; ++i){
        for(int j = 0; j < size; ++j){
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

    undisturbedPoints = m_shape.getVertices();
    newPoints.resize(undisturbedPoints.size());
}

// initializes an obj to drop on the water surface mesh
void Simulation::initFallingParticleSystem(string meshPath, Vector3f startPos) {
    // meshPath: path to object that you want to drop on the water surface
    // startPos: the starting point for the object to start falling from

    vector<Vector3f> vertices; // gets filled in by loadTriMesh
    vector<Vector3i> triangles; // gets filled in by loadTriMesh

    System* sys = nullptr;

    if (MeshLoader::loadTriMesh(meshPath, vertices, triangles)) {
        sys = new System(vertices, triangles); // init a new particle system
        sys->setParticleMass(config.objMass);

        Shape* fallingShape = new Shape();
        fallingShape->init(vertices, triangles);

        Eigen::Affine3f modelMatrix = Eigen::Affine3f::Identity();
        modelMatrix.scale(Eigen::Vector3f(config.objUniformScale, config.objUniformScale, config.objUniformScale)); // need to scale the obj coz the water mesh is way too large
        modelMatrix.translation() = Eigen::Vector3f(startPos[0], startPos[1], startPos[2]);
        fallingShape->setModelMatrix(modelMatrix); // y-up axis

        sys->setFallingShape(fallingShape); // set the shape this particle system instance is trying to simulate at its core
        sys->setWaterSurfaceShape(&m_shape); // give the falling shape access to water surface mesh shape for solid-fluid coupling

        m_fallingShapes.push_back(fallingShape); // global list of all falling shapes in the scene
        m_particleSystems.push_back(sys); // global list of all particle system instances
    } else {
        assert(!"Fail to load mesh obj for the falling object D:");
    }
}

