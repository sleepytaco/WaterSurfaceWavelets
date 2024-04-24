#include "solver.h"

Solver::Solver() { }


void Solver::eulerStep(System* sys, double h) {
    // input h is a particular timestep forward in future
    double t = sys->getTime();
    MatrixXf x0 = sys->getState();
    MatrixXf deltaX = sys->getStateDerivative(x0); // derivEval returns derivative of current state x0

    // update particles state to move a timestep h forward
    sys->setState(x0 + h*deltaX);
    sys->setTime(t+h);
}

void Solver::midPointMethodStep(System* sys, double h) {
    // input h is a particular timestep forward in future
    double t0 = sys->getTime();
    MatrixXf x0 = sys->getState(); // contains (x, v)
    MatrixXf deltaX0 = sys->getStateDerivative(x0); // derivEval returns derivative of current state x0 --- contains (x_dot, v_dot)

    MatrixXf x1 = x0 + 0.5*h*deltaX0; // state after half step
    // double t1 = t0 + 0.5*h; // time after half step
    // the sys is now at time t1
    // calc and save x' and v' at time t1
    MatrixXf deltaX1 = sys->getStateDerivative(x1);

    // update particles state to move a FULL timestep h forward
    sys->setState(x0 + h*deltaX1);
    sys->setTime(t0+h);
}

void Solver::RK4(System* sys, double h) {
    // https://lpsa.swarthmore.edu/NumInt/NumIntFourth.html RK4 algo notes

    // input h is a particular timestep forward in future
    double t0 = sys->getTime();
    MatrixXf x0 = sys->getState(); // contains (x, v)
    MatrixXf k1 = sys->getStateDerivative(x0); // derivEval returns derivative of current state x0 --- contains (x_dot, v_dot)

    MatrixXf temp = x0 + 0.5*h*k1;
    MatrixXf k2 = sys->getStateDerivative(temp);
    temp = x0 + 0.5*h*k2;
    MatrixXf k3 = sys->getStateDerivative(temp);
    temp = x0 + h*k3;
    MatrixXf k4 = sys->getStateDerivative(temp);

    // update particles state to move a timestep h forward
    sys->setState(x0 + (1/(double)6 * (k1+2*k2+2*k3+k4) * h));
    sys->setTime(t0 + h);
}
