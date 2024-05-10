//#include "config.h"
//#include "Eigen/Dense"
//using namespace Eigen;

//enum LevelSetShape { Circle, Square };

//// not exactly a levelset rn but we ball
//// this is reflecting stuff based off implicit equations
//class LevelSet
//{
//public:
//    // must specify levelset shape during creation
//    LevelSet(LevelSetShape shape);

//    int shape = LevelSetShape::Circle;
//    Vector2d shapeCenter = Vector2d(config.dXY/2, config.dXY/2);
//    void setShapeCenter(Vector2d center) {shapeCenter = center;};

//    // default circle params
//    double radius = 1;
//    void setCircleParams(double r) {radius = r;};

//    // default square params
//    // or, width, height, center of square --- can implicitly define the boundaries this way...?
//    double squareWidth = 1;
//    double squareHeight = 1;
//    void setSquareParams(double w, double h) {squareWidth=w; squareHeight=h;};

//    void circleReflection(Vector2d& advPos, int& thetaIdx);
//    void squareReflection(Vector2d& advPos, int& thetaIdx);

//    // aggregates all the possible reflections
//    void boundaryReflection(Vector2d& advPos, int& thetaIdx);
//};
