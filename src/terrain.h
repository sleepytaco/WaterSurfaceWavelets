#ifndef TERRAIN_H
#define TERRAIN_H

#include "config.h"
#include "Eigen/Dense"
#include "graphics/shape.h"


using namespace Eigen;

class Terrain
{
public:
    Terrain();
    ~Terrain();

    void init();

    int resolution = config.dimXY;
    int m_lookup;

    Shape m_shape;

    void draw(Shader *shader, GLenum mode)
    {
        m_shape.draw(shader, mode);
    }

    std::vector<Vector3f>& get(){return m_vertices;}

private:
    Vector2d sample_vector(int x, int y);

    Vector3f getPosition(int x, int y);

    double getHeight(double x, double y);

    // Computes color of vertex using normal and, optionally, position
    Vector3d getColor(Vector3f normal, Vector3f position);

    // Computes the intensity of Perlin noise at some point
    double computePerlin(double x, double y);

    void generateTerrain();

    std::vector<Vector2d> randlookup;

    std::vector<Vector3f> m_vertices;


};

#endif // TERRAIN_H
