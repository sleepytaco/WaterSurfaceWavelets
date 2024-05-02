#version 330 core
// from https://learnopengl.com/Advanced-OpenGL/Cubemaps
layout (location = 0) in vec3 position;

out vec3 uvwCoord;

uniform mat4 u_projectionMatrix;
uniform mat4 u_viewMatrix;

void main() {

    uvwCoord = vec3(position.x, position.y, position.z);
    gl_Position = u_projectionMatrix * u_viewMatrix * vec4(position, 1.0);
}
