#pragma once

#include "amplitude.h"
#include "graphics/shape.h"
#include "Eigen/StdList"
#include "Eigen/StdVector"
#include "grid.h"
#include <random>
#include "FEM/solver.h"

class Shader;

class Simulation
{
private:
    Shape m_shape; // this is the water surface mesh (scared to rename this as the ARAP boilerplate relies on this literally everywhere)
    Amplitude m_amplitude;
    std::uniform_real_distribution<double> unif;
    std::default_random_engine re;
    double lower_bound;
    double upper_bound;
    std::vector<Eigen::Vector3f> undisturbedPoints;
    std::vector<Eigen::Vector3f> newPoints;


    // ============== Solid-fluid Coupling related stuff
    std::vector<Shape*> m_fallingShapes; // global list of all falling objs in the scene
    std::vector<System*> m_particleSystems; // global list of all particle system instances for falling objs in the scene

    void initFallingParticleSystem(std::string meshPath, Vector3f startPos); // initializes a falling object in the scene and adds it to the two global lists above

    Solver* solver; // used to integrate the all m_particleSystems states forward in time (contains eulerstep, midpointstep, rk4 integrators)
    // ===========================================

public:
    Simulation();

    void init(Eigen::Vector3f &min, Eigen::Vector3f &max);
    void move(int vertex, Eigen::Vector3f pos);

    void setWaterHeights();
    void update(double deltaTime);


    // ================== If You Choose To Modify The Code Below, It's On You

    int getClosestVertex(Eigen::Vector3f start, Eigen::Vector3f ray, float threshold)
    {
        return m_shape.getClosestVertex(start, ray, threshold);
    }

    void draw(Shader *shader, GLenum mode)
    {
        m_shape.draw(shader, mode); // water surface mesh
        for (int i=0; i<m_fallingShapes.size(); ++i) { // draw all the objs thrown on the water surface mesh
            m_fallingShapes[i]->draw(shader, mode);
        }
    }

    SelectMode select(Shader *shader, int vertex)
    {
        return m_shape.select(shader, vertex);
    }

    bool selectWithSpecifiedMode(Shader *shader, int vertex, SelectMode mode)
    {
        return m_shape.selectWithSpecifiedMode(shader, vertex, mode);
    }

    bool getAnchorPos(int lastSelected, Eigen::Vector3f& pos, Eigen::Vector3f ray, Eigen::Vector3f start)
    {
        return m_shape.getAnchorPos(lastSelected, pos, ray, start);
    }
};
