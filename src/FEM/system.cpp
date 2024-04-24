#include "system.h"

System::System() {

}

System::System(vector<Vector3f>& vertices, vector<Vector3i>& faces) {

    numParticles = vertices.size();
    _particleStates = MatrixXf(6, numParticles); // stores matrial space pos and velocity
    _particleForces = MatrixXf(3, numParticles);

    _faces = faces; // not really using this tbh but too lazy to remove it...

    for (int i=0; i<numParticles; ++i) {
        VectorXf state(6);
        state << vertices[i],  // particle initial position
                 Vector3f(0, 0, 0); // particle initial velocity
        _particleStates.col(i) = state;
    }
}

// derivEval loop from slides
MatrixXf System::getStateDerivative(MatrixXf& currParticleStates) {
    MatrixXf particleStateDerivatives(6, numParticles);

    MatrixXf currParticleForces = _particleForces;
    MatrixXf particleAcc = currParticleForces / particleMass;

    // STEP 1: clear forces
    _particleForces.setZero();

    // STEP 2: calculate external and internal forces - sum all forces into accumulators
    _particleForces = calculateForces(currParticleStates);

    // STEP 3: gather - copy v and f/m into destination array
    for (int i=0; i<numParticles; ++i) {
        VectorXf stateDerivative(6);
        Vector3f particleVelocity = currParticleStates.col(i).tail(3); // get the last 3 components from particle state which contains (pos3d, vel3d)
        stateDerivative << particleVelocity, // v
                           _particleForces.col(i) / particleMass; // f/m
        particleStateDerivatives.col(i) = stateDerivative;
    }

    return particleStateDerivatives;
}

MatrixXf System::calculateForces(MatrixXf& currParticleStates) {
    MatrixXf forceAccumulator(3, numParticles);
    forceAccumulator.setZero();

    // apply external forces
    for (int i=0; i<numParticles; ++i) {
        // add external forces
        forceAccumulator.col(i) += particleMass * _g; // force due to gravity

        // add viscous drag too maybe ???
        Vector3f particleVelocity = currParticleStates.col(i).tail(3); // get the last 3 components from particle state which contains (pos3d, vel3d)
        double kd = 0.1; // 1/4e4;
        // Vector3d a = particleAcc.col(i);
        forceAccumulator.col(i) -= kd * particleVelocity; // * particleMass; // force prop to opp direction of particle's velocity it is moving in

        // force due to the particle hitting/go beyond ground (if any)
        // includes ground friction
        // applyGroundCollisionForces(i, currParticleStates, forceAccumulator);
        // applyFixedSphereCollisionForces(i, currParticleStates, forceAccumulator);
        // applyInterShapeCollisionForces(i, currParticleStates, forceAccumulator);
    }

    // apply internal forces to all particles/nodes
    // applyStressStrainForces(currParticleStates, forceAccumulator);

    return forceAccumulator;
}

vector<Vector3f> System::getVertices() {
    vector<Vector3f> vertices;
    for (int i=0; i<numParticles; ++i) {
        Vector3f state = _particleStates.col(i).head(3); // get the first 3 components from particle state which contains (pos3d, vel3d)
        vertices.push_back(state);
    }
    return vertices;
}
