#include "amplitude.h"

void Amplitude::init_ls(std::vector<Vector3f>& terrain){
    this->terrain = terrain;
}


std::optional<Vector2d> Amplitude::boundary_check(Vector2d pos){

    Vector2d idxPos = posToIdxSpace(pos);

    if(idxPos(0)>=0 && idxPos(1)>=0 && idxPos(0) < config.dimXY -1
        && idxPos(1) < config.dimXY - 1){
        int x1 = floor(idxPos(0));
        int y1 = floor(idxPos(1));
        int x2 = x1 + 1;
        int y2 = y1 + 1;
        double dx = idxPos(0) - x1;
        double dy = idxPos(1) - y1;

        Vector3f q11 = this->terrain[x1 + config.dimXY * y1];
        Vector3f q21 = this->terrain[x2 + config.dimXY * y1];
        Vector3f q12 = this->terrain[x1 + config.dimXY * y2];
        Vector3f q22 = this->terrain[x2 + config.dimXY * y2];

        double interpolatedValue = (q11 * (1 - dx) * (1 - dy) +
                                   q12 * (1 - dx) * dy +
                                   q21 * dx * (1 - dy) +
                                    q22 * dx * dy).y();

        if(interpolatedValue > 0){
            double d1 = (q21 - q11).y();
            double d2 = (q12 - q11).y();
            Vector3f v1 = Vector3f(5, d1, 0);
            Vector3f v2 = Vector3f(0, d2, 5);
            Vector3f normal = v1.cross(v2);
            Vector2d normal2d = Vector2d(normal.x(), normal.z()).normalized();
            return std::optional<Vector2d>{normal2d};
        }

    }
    return std::nullopt;

}
