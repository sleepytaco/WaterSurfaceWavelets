#pragma once

#include <GL/glew.h>
#include <vector>
#include <unordered_set>

#define EIGEN_DONT_VECTORIZE
#define EIGEN_DISABLE_UNALIGNED_ARRAY_ASSERT
#include "Eigen/StdVector"
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Eigen::Matrix2f)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Eigen::Matrix3f)
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Eigen::Matrix3i)
#include "Eigen/Dense"

enum SelectMode
{
    None     = 0,
    Anchor   = 1,
    Unanchor = 2
};

class Shader;

class Shape
{
public:
    Shape();

    void init(const std::vector<Eigen::Vector3f> &vertices, const std::vector<Eigen::Vector3i> &triangles);
    void setVertices(const std::vector<Eigen::Vector3f> &vertices);
    void setVertices(const std::vector<Eigen::Vector3f> &vertices, const std::vector<Eigen::Vector3f> &normals);


    void setModelMatrix(const Eigen::Affine3f &model);
    void setModelMatrix(const Eigen::Matrix4f &model);
    void setupSkybox();

    void draw(Shader *shader, GLenum mode);
    SelectMode select(Shader *shader, int vertex);
    bool selectWithSpecifiedMode(Shader *shader, int vertex, SelectMode mode);
    int  getClosestVertex(Eigen::Vector3f start, Eigen::Vector3f ray, float threshold);
    bool getAnchorPos(int lastSelected, Eigen::Vector3f& pos, Eigen::Vector3f ray, Eigen::Vector3f start);

    const std::vector<Eigen::Vector3f>& getVertices();
    const std::vector<Eigen::Vector3i>& getFaces();
    const std::unordered_set<int>& getAnchors();
    Eigen::Matrix4f& getModelMatrix() { return m_modelMatrix;};

    // making an assumption that all shapes are approximated by a sphere
    Eigen::Vector3f getShapeCentroid() {return shapeCentroid;};

    void setColor(float r, float g, float b) {m_red=r; m_blue=b; m_green=g;}

private:
    GLuint m_surfaceVao;
    GLuint m_surfaceVbo;
    GLuint m_surfaceIbo;

    unsigned int m_numSurfaceVertices;
    unsigned int m_verticesSize;
    float m_red;
    float m_blue;
    float m_green;
    float m_alpha;

    std::vector<Eigen::Vector3i> m_faces;
    std::vector<Eigen::Vector3f> m_vertices;
    std::unordered_set<int>      m_anchors;

    Eigen::Matrix4f m_modelMatrix;
    int lastSelected = -1;

    // Helpers
    // making an assumption that all shapes are approximated by a sphere
    Eigen::Vector3f shapeCentroid; // avg of all vertices
    void updateShapeCentroid(std::vector<Eigen::Vector3f>& verts); // recompute the mean of all vertices and set it to be shapeCentroid

    void selectHelper();
    Eigen::Vector3f getNormal(const Eigen::Vector3i& face);
    void updateMesh(const std::vector<Eigen::Vector3i> &triangles,
                    const std::vector<Eigen::Vector3f> &vertices,
                           std::vector<Eigen::Vector3f>& verts,
                           std::vector<Eigen::Vector3f>& normals,
                           std::vector<Eigen::Vector3f>& colors);

    void updateMesh(const std::vector<Eigen::Vector3i> &triangles,
                    const std::vector<Eigen::Vector3f> &vertices,
                    const std::vector<Eigen::Vector3f> &normals,
                    std::vector<Eigen::Vector3f>& verts,
                    std::vector<Eigen::Vector3f>& norms,
                    std::vector<Eigen::Vector3f>& colors);
};
