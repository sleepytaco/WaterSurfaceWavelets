#version 330 core
// from https://learnopengl.com/Advanced-OpenGL/Cubemaps
out vec4 fragColor;

in vec3 uvwCoord;

uniform samplerCube u_skybox;

void main() { fragColor = texture(u_skybox, uvwCoord); }
