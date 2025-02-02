#version 460 core

in vec4 fragPos;
in vec3 normal;

uniform vec3 lightDirection;
uniform vec4 position;
uniform vec4 drawColor;

out vec4 finalColor;

void main() {
	finalColor = drawColor;
	
	if (drawColor.x == 0) {
		if (pow(fragPos.x - position.x, 2) + pow(fragPos.z - position.z, 2) > 4000) {
			finalColor.a = 0;
		}
		return;
	}

	float diffuse = (max(dot(normalize(-normal), lightDirection), -1) + 1) / 2;
	finalColor.xyz = finalColor.xyz * diffuse;
}