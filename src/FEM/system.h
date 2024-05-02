#include "Eigen/Dense"
#include "iostream"
#include "graphics/shape.h"
#include <unordered_set>
#include "amplitude.h"

using namespace Eigen;

class System {
public:
    System();
    System(std::vector<Vector3f>& vertices, std::vector<Vector3i>& faces);

    double getTime() {return t;};
    MatrixXf getState() { return _particleStates; };
    MatrixXf getStateDerivative(MatrixXf& currParticleStates); // derivEval loop from slides

    void setTime(double tNew) { t = tNew;};
    void setState(MatrixXf state) { _particleStates = state; };
    void setParticleMass(double pMass) {particleMass = pMass/numParticles;};

    std::vector<Vector3f> getVertices(); // extracts particle positions from particleStates variable and returns them

    // collisions
    void setFallingShapesList(std::vector<Shape*>& fsArray) {_fallingShapes = fsArray;};
    void setFallingShape(Shape* shape) {fallingShape = shape;};
    void setWaterSurfaceShape(Shape* shape) {waterSurfaceShape = shape;};

    void setAmplitudeFunction(Amplitude* amp) { _amplitude4d = amp;};

private:
    const int particleStateDim = 6; // dim of an individual particle (position3f + velocity3f)
    double t = 0; // current time in system
    const Vector3f _g = Vector3f(0.0, -9.8, 0.0); // const force of gravity

    // these are filled in during initialize()
    int numParticles;
    double particleMass; // assuming all particles share the mass
    MatrixXf _particleStates; // will have size (6 x numparticles) where 6 is state dim
    MatrixXf _particleForces; // will have size (3 x numparticles) -- this matrix hold force accumultor for each particle
    std::vector<Vector3i> _faces;

    Shape* fallingShape; // this is the core object/shape whose state is being simulated by the system
    std::vector<Shape*> _fallingShapes; // global list of all (other) falling shapes in-scene -- useful for inter-shape collisions in the future

    MatrixXf calculateForces(MatrixXf& currParticleStates);

    Shape* waterSurfaceShape; // get access to water surface mesh
    Amplitude* _amplitude4d;

    void applyBuoyancyForces(int particleId, MatrixXf& currParticleStates, MatrixXf& forceAccumulator); // apply BuoyancyForces induced by water waves on the particle
    Vector3f calculateBuoyancyForce(MatrixXf& currParticleStates);

    // utils
    Vector3f getWorldSpacePos(Vector3f pos, Matrix4f modelMat);
    Vector3f getWorldSpaceDir(Vector3f dir, Matrix4f modelMat);
};
