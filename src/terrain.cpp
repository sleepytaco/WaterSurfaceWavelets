#include "terrain.h"

Terrain::Terrain()
{
    m_lookup = 1024;
    randlookup.reserve(m_lookup);

    // Initialize random number generator
    std::srand(1230);

    // Populate random vector lookup table
    for (int i = 0; i < m_lookup; i++)
    {
        randlookup.push_back(Vector2d(std::rand() * 2.0 / RAND_MAX - 1.0,
                                            std::rand() * 2.0 / RAND_MAX - 1.0));
    }
}

Terrain::~Terrain()
{
    randlookup.clear();
}

void Terrain::init(){
    std::vector<Vector3i> triangles;

    for(int i = 0; i < resolution; ++i){
        for(int j = 0; j < resolution; ++j){
            m_vertices.push_back(getPosition(i, j));
        }
    }

    for(int i = 0; i < resolution - 1; ++i){
        for(int j = 0; j < resolution - 1; ++j){
            int c1 = i + j * resolution;
            int c2 = i + 1 + j * resolution;
            int c3 = i + (j + 1) * resolution;
            int c4 = i + 1 + (j + 1) * resolution;

            triangles.push_back(Vector3i(c3, c4, c1));
            triangles.push_back(Vector3i(c4, c2, c1));
        }

    }

    this->m_shape.init(m_vertices, triangles);
    m_shape.setColor(0.65,0.41,0.33);

}


Vector3f Terrain::getPosition(int x, int y) {
    // Normalizing the planar coordinates to a unit square
    // makes scaling independent of sampling resolution.
    float xn = (1.0 * x) / resolution;
    float yn = (1.0 * y) / resolution;
    float z = getHeight(xn, yn) * resolution * config.meshScale;
    float half = resolution * config.meshScale * 0.5;
    return Vector3f(y * config.meshScale - half,
                    z,
                    x * config.meshScale - half);
}

double Terrain::getHeight(double x, double y){
    double z_2 = computePerlin(x * 2, y * 2) / 2;
    double z_4 = computePerlin(x * 4, y * 4)/4;
    double z_8 = computePerlin(x * 8, y * 8)/16;
    double z_16 = computePerlin(x * 16, y * 16)/16;

    // Task 7: combine multiple different octaves of noise to produce fractal perlin noise

    // Return 0 as placeholder
    return z_2 + z_4 + z_8 + z_16;
}

float interpolate(float A, float B, float alpha) {
    // Task 4: implement your easing/interpolation function below
    float ease_output = 3 * pow(alpha, 2) - 2 * pow(alpha, 3);

    // Return 0 as placeholder
    return A + ease_output * (B - A);
}




double Terrain::computePerlin(double x, double y){
    int x_low = floor(x);
    int x_high = x_low + 1;
    int y_low = floor(y);
    int y_high = y_low + 1;

    Vector2d rand_ll = sample_vector(y_low, x_low);
    Vector2d rand_lh = sample_vector(y_high, x_low);
    Vector2d rand_hl = sample_vector(y_low, x_high);
    Vector2d rand_hh = sample_vector(y_high, x_high);



    // Task 2: compute offset vectors
    Vector2d offset_ll =  Vector2d(x - x_low, y - y_low);
    Vector2d offset_hh =  Vector2d(x - x_high, y - y_high);
    Vector2d offset_lh =  Vector2d(x - x_low, y - y_high);
    Vector2d offset_hl =  Vector2d(x - x_high, y - y_low);

    // Task 3: compute the dot product between the grid point direction vectors and its offset vectors
    // float A = ... // dot product between top-left direction and its offset
    float A = rand_ll.dot(offset_ll);
    float B = rand_hl.dot(offset_hl);
    float C = rand_lh.dot(offset_lh);
    float D = rand_hh.dot(offset_hh);

    return interpolate(interpolate(A, B, x - x_low), interpolate(C, D, x - x_low), y- y_low);
}

Vector2d Terrain::sample_vector(int x, int y)
{
    std::hash<int> intHash;
    int index = intHash(x * 41 + y * 43) % m_lookup;
    return randlookup.at(index);
}
