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

        float r1 = q11(1) * (q21(0) - idxPos(0))/5.0
                      + q21(1) * (idxPos(0) - q11(0))/5.0;

        float r2 = q12(1) * (q21(0) - idxPos(0))/5.0
                      + q22(1) * (idxPos(0) - q11(0))/5.0;

        float p = r1 * (q21(2) - idxPos(1))/5.0
                   + r2 * (idxPos(1)  - q11(2))/5.0;

        if(interpolatedValue > 0){
            Vector3f d1 = q21 - q11;
            Vector3f d2 = q12 - q11;
            Vector3f cross = d1.cross(d2);
            Vector2d normal = Vector2d(cross(0), cross(2)).normalized();
            return std::optional<Vector2d>{normal};
        }

    }
    return std::nullopt;

}
