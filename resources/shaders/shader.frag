#version 330 core
out vec4 fragColor;

in vec3 worldPos;
in vec3 worldNorm;
in vec3 normal_cameraSpace;

uniform vec4 cameraPos; // position of the camera

uniform samplerCube u_skybox;
uniform float u_reflection;
uniform float u_refraction;
uniform float u_materialRefractiveIndex;

uniform int   wire  = 0;
uniform float red   = 1.0;
uniform float green = 1.0;
uniform float blue  = 1.0;
uniform float alpha = 1.0;
uniform int drawmode = 0;

// ---------- SOURCE: https://people.computing.clemson.edu/~jtessen/reports/papers_files/coursenotes2004.pdf ----------
uniform vec4 upwelling = vec4(0.0, 0.2, 0.3, 1.0);
uniform vec4 sky = vec4(0.69, 0.84, 1.0, 1.0);
uniform vec4 air = vec4(0.1, 0.1, 0.1, 1.0);
uniform float nSnell = 1.34;
uniform float Kdiffuse = 0.91;

void main() {
    if (drawmode == 0) {
        vec3 nI = normalize(worldPos - vec3(cameraPos));
        vec3 nN = normalize(worldNorm);

        float costhetai = abs(dot(nI, nN));
        float thetai = acos(costhetai);
        float sinthetat = sin(thetai) / nSnell;
        float thetat = asin(sinthetat);

        float reflectivity;

        if (thetai == 0.0) {
            reflectivity = (nSnell - 1) / (nSnell + 1);
            reflectivity = reflectivity * reflectivity;
        } else {
            float fs = sin(thetat - thetai) / sin(thetat + thetai);
            float ts = tan(thetat - thetai) / tan(thetat + thetai);
            reflectivity = 0.5 * (fs * fs + ts * ts);
        }

        vec3 dPE = worldPos - vec3(cameraPos);
        float dist = /*length(dPE) **/ Kdiffuse;
        dist = exp(-dist);

        vec4 skyColor = sky;

        skyColor = texture(u_skybox, nN);

        fragColor = dist * (reflectivity * skyColor + (1 - reflectivity) * upwelling) + (1 - dist) * air;
    }

    else if (drawmode == 1) {
        vec3 lightDir = normalize(vec3(0, 0.5, 1));
        float c = clamp(dot(normal_cameraSpace, lightDir), 0, 1);
        fragColor = vec4(red * c, green * c, blue * c, 1);
    }
}

//void main() {

//    vec3 normalizedWorldNorm = normalize(worldNorm);

//    vec3 vecToPos = normalize(worldPos - vec3(cameraPos));
//    float reflectance = pow((1.f-u_materialRefractiveIndex) / (1.f+u_materialRefractiveIndex), 2);

//    float fresnel = max(reflectance + (1.f-reflectance) * pow((1.f-cos(dot(normalizedWorldNorm, -vecToPos))), 5.f), 0.1f);

//    vec3 refractVec = normalize(refract(vecToPos, normalizedWorldNorm, 1.f/u_materialRefractiveIndex));
//    vec3 reflectVec = normalize(reflect(vecToPos, normalizedWorldNorm));

//    fragColor = fresnel * u_reflection * texture(u_skybox, reflectVec) + (1.f-fresnel) * u_refraction * texture(u_skybox, refractVec);

////    // Do lighting in camera space
////    vec3 lightDir = normalize(vec3(0, 0.5, 1));
////    float c = clamp(dot(normal_cameraSpace, lightDir), 0, 1);

////    fragColor = vec4(red * c, green * c, blue * c, 1);
//}
