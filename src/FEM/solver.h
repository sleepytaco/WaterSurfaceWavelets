#include "system.h"

class Solver {
public:
    Solver();

    void eulerStep(System* sys, double h);
    void midPointMethodStep(System* sys, double h);
    void RK4(System* sys, double h);
};
