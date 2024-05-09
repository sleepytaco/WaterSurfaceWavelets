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

    config.driveForce = Vector3d::Zero();
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
        // force due to gravity
        forceAccumulator.col(i) += particleMass * _g;

        // force due to buoyancy
        forceAccumulator.col(i) += buoyancyForce;

        Vector3f particleVelocity = currParticleStates.col(i).tail(3); // get the last 3 components from particle state which contains (pos3d, vel3d)
        double kd = 0.05; // 1/4e4;
        // Vector3d a = particleAcc.col(i);
        forceAccumulator.col(i) -= kd * particleVelocity;// * particleMass; // force prop to opp direction of particle's velocity it is moving in

        Eigen::Matrix3f mat = AngleAxisf(config.boatRotation, Vector3f::UnitY()).toRotationMatrix();
        forceAccumulator.col(i) += mat*config.driveForce.cast<float>();

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
    float shapeRadius = 5; // TODO: get from config file or calc it

    // get smallest/lowest y-coord of all the vertices in the falling shape obj
    float smallestShapeY = particlePos.y();
    for (int i=1; i<numParticles; ++i) {
        particlePosMatSpace = currParticleStates.col(i).head(3); // retrieve the relevant particle position from state
        particlePos = getWorldSpacePos(particlePosMatSpace, fallingShape->getModelMatrix()); // particle pos in world space
        smallestShapeY = std::min(smallestShapeY, particlePos.y());
    }

    Vector2d particleSurfacePos = Vector2d(particlePos.x() / (config.dimXY*config.meshScale) * (config.xMax - config.xMin), particlePos.z() / (config.dimXY*config.meshScale) * (config.yMax - config.yMin));
    auto [displacement, normal] = _amplitude4d->waterHeight(particleSurfacePos);
    float waterSurfaceY = displacement.y();

    if (smallestShapeY - waterSurfaceY > 0) { // this means the entirity of the shape is above the water surface
        return buoyancyForce; // return 0 force
    }

    // formula from paper
    float h = fmin(waterSurfaceY - smallestShapeY, shapeRadius); // TODO: don't fully get what this is
    float V = M_PI * shapeRadius * shapeRadius * h; // "volume" of submerged portion
    float rho = config.fluidDensity;
    buoyancyForce = -_g * rho * V;

    Vector2d idxSpacePos = _amplitude4d->posToIdxSpace(particleSurfacePos);
    double idxSpaceX = idxSpacePos.x();
    double idxSpaceY = idxSpacePos.y();
    int i = floor(idxSpaceX);
    int j = floor(idxSpaceY);

    Vector3f particleVelocity = currParticleStates.col(0).tail(3);
    double rigidEnergy = config.objMass * particleVelocity.squaredNorm() * 0.5 + config.objMass * config.g * h;
    if (prevRigidEnergy == -1) prevRigidEnergy = rigidEnergy;
    double rigidEnergyDelta = rigidEnergy - prevRigidEnergy;
    rigidEnergyDelta = rigidEnergyDelta;
    prevRigidEnergy = rigidEnergy;

    for (int theta=0; theta<=config.dimTheta; ++theta) { // b
        double currentAmp = _amplitude4d->m_currentAmplitude.get(i, j, theta, 0);
//        double fluidEnergy = 0.5 * rho * config.g * currentAmp * currentAmp;
        double newAmp = 2 / (rho * 1000 * config.g) * (rigidEnergyDelta / config.dimTheta);
        _amplitude4d->m_currentAmplitude.set(i, j, theta, 0, newAmp);
    }

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

void System::rotateBoat(int direction) {
    config.boatRotation += direction * 0.01 * M_PI;
    float totalX = 0;
    float totalZ = 0;
    for (int i=0; i<numParticles; ++i) {
        Vector3f pos = _particleStates.col(i).head(3);
        totalX += pos.x();
        totalZ += pos.z();
    }
    Vector3f center = Vector3f(totalX/numParticles, 0, totalZ/numParticles);

    for (int i=0; i<numParticles; ++i) {
        Vector3f pos = _particleStates.col(i).head(3);
        Eigen::Matrix3f mat = AngleAxisf(direction * 0.01 * M_PI, Vector3f::UnitY()).toRotationMatrix();
        _particleStates.col(i).head(3) = mat * (pos - center) + center;
    }
}
