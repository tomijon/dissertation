#version 460 core

layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 vertNormal;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

out vec4 fragPos;
out vec3 normal;

void main() {
	mat4 mvp = projection * view * model;
	gl_Position = mvp * vec4(aPos, 1.0);
	fragPos = vec4(aPos, 1.0);
	normal = vertNormal;
}