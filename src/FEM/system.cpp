#include "system.h"
#include "config.h"

System::System() {

}

System::System(std::vector<Vector3f>& vertices, std::vector<Vector3i>& faces) {

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
    Vector3f buoyancyForce = calculateBuoyancyForce(currParticleStates);

    for (int i=0; i<numParticles; ++i) {
        // add external forces
        forceAccumulator.col(i) += particleMass * _g; // force due to gravity

        // add viscous drag too maybe ???
        Vector3f particleVelocity = currParticleStates.col(i).tail(3); // get the last 3 components from particle state which contains (pos3d, vel3d)
        double kd = 0.1; // 1/4e4;
        // Vector3d a = particleAcc.col(i);
        forceAccumulator.col(i) -= kd * particleVelocity;// * particleMass; // force prop to opp direction of particle's velocity it is moving in

        // force due to the particle hitting the water surface
        forceAccumulator.col(i) += buoyancyForce;

    }

    // apply internal forces to all particles/nodes
    // applyStressStrainForces(currParticleStates, forceAccumulator);

    return forceAccumulator;
}

Vector3f System::getWorldSpacePos(Vector3f pos, Matrix4f modelMat) {
    Vector4f particlePos4f;
    particlePos4f << pos, 1.0;
    Vector4f particlePosWorld4f = (modelMat * particlePos4f);
    return particlePosWorld4f.head(3);
}

Vector3f System::getWorldSpaceDir(Vector3f dir, Matrix4f modelMat) {
    Vector4f particleDir4f;
    particleDir4f << dir, 0.0;
    Vector4f particleDirWorld4f = (modelMat * particleDir4f);
    return particleDirWorld4f.head(3);
}


Vector3f System::calculateBuoyancyForce(MatrixXf& currParticleStates) {

    Vector3f buoyancyForce(0, 0, 0);
    Vector3f shapeCOMMatSpace = fallingShape->getShapeCentroid(); // center of mass of the falling shape
    Vector3f shapeCOM = getWorldSpacePos(shapeCOMMatSpace, fallingShape->getModelMatrix()); // COM in world space

    Vector3f particlePosMatSpace = currParticleStates.col(0).head(3); // retrieve the relevant particle position from state
    Vector3f particlePos = getWorldSpacePos(particlePosMatSpace, fallingShape->getModelMatrix()); // particle pos in world space
    float shapeRadius = (particlePos - shapeCOM).norm(); // TODO: get from config file or calc it

    // get smallest/lowest y-coord of all the vertices in the falling shape obj
    float smallestShapeY = particlePos.y();
    for (int i=1; i<numParticles; ++i) {
        particlePosMatSpace = currParticleStates.col(i).head(3); // retrieve the relevant particle position from state
        particlePos = getWorldSpacePos(particlePosMatSpace, fallingShape->getModelMatrix()); // particle pos in world space
        smallestShapeY = std::min(smallestShapeY, particlePos.y());
    }

    // get largest/highest y-coord of all the vertices in the water surface mesh
    float largestWaterSurfaceY = waterSurfaceShape->getVertices()[0].y();
    for (const Vector3f& v : waterSurfaceShape->getVertices()) {
        largestWaterSurfaceY = std::max(largestWaterSurfaceY, v.y());
    }

    if (smallestShapeY - largestWaterSurfaceY > -shapeRadius/4) { // this means the entirity of the shape is above the water surface
        return buoyancyForce; // return 0 force
    }

    // formula from paper
//    float h = (largestWaterSurfaceY - smallestShapeY); // TODO: don't fully get what this is
//    float waterHeight = largestWaterSurfaceY;
//    float V = M_PI * shapeRadius * shapeRadius * (waterHeight - shapeRadius - h); // "volume" of submerged portion
//    float rho = 1; config.fluidDensity;
//    buoyancyForce = -1 * rho * V * _g;

    // naive way
     buoyancyForce = Vector3f(0, 1, 0); // * (largestWaterSurfaceY - smallestShapeY);

    return buoyancyForce;
}

std::vector<Vector3f> System::getVertices() {
    std::vector<Vector3f> vertices;
    for (int i=0; i<numParticles; ++i) {
        Vector3f state = _particleStates.col(i).head(3); // get the first 3 components from particle state which contains (pos3d, vel3d)
        vertices.push_back(state);
    }
    return vertices;
}
