#version 460 core

#define PERLIN 1
#define VALUE 2
#define LAYERED 3
#define pi 3.1415926535897932384626433832795

layout(local_size_x = 32, local_size_y = 1, local_size_z = 32) in;

// Uniforms.
uniform int n_chunks;
uniform int chunkWidth;
uniform int chunkHeight;
uniform int style;

// Chunk coordinate input buffer.
layout(std430, binding = 0) buffer Coordinates {
    int coordinates[];
};

// Data output buffer.
layout(std430, binding = 1) buffer Data {
    float data[];
};


// Convert the chunk number, x and z to an index in the data array.
int to_data_index(int chunk, int x, int z) {
    return (chunk * chunkWidth * chunkWidth) + (z * chunkWidth) + x;
}


// -- NOISE GENERATION FUNCTIONS ---------------------------------------
float lerp(float v1, float v2, float weight) {
    return v1 - ((v1 - v2) * weight);
}

// Random number generator for value noise.
float value_noise(int x, int z, float seed) {
    float hash = pow(x, 2) + (3 * x) + (2 * x * z) + z + pow(z, 2);
    hash /= 2.0f;
    float in_cos = cos(hash + seed);
    float in_sin = sin(hash + seed);
    return (-1 * ((pow(4, in_cos) + (6 * pow(in_sin, 2))) / 4.0f) + 2) * 0.5;
}


// Interpolation weight function.
float weight(float value) {
    int range = int(value / ((chunkWidth - 1) / 2.0f));
    int cos_mult = -int(((range % 2) - 0.5f) * 2);
    int addition = (2 * range) + 1;

    float scalar = (2 * pi) / (chunkWidth - 1);
    return ((cos_mult * cos((scalar * value) + pi)) + addition) / 4;
}


// Generate a value for the corner.
float[4] vertex_noise(int chunkX, int chunkZ, float seed) {
    float corners[4];

    int offset = 0;
    for (int vert = 0; vert < 4; vert++) {
        int x_diff = offset & 1;
        int z_diff = (offset & 2) >> 1;
        corners[vert] = value_noise(chunkX + x_diff, chunkZ + z_diff, seed);
        offset++;
    }
    return corners;
}


// Generate the noise.
float gen_value_noise(int chunkX, int chunkZ, int x, int z, float seed) {
    float corners[4] = vertex_noise(chunkX, chunkZ, seed);

    float v01 = lerp(corners[0], corners[1], weight(x));
    float v23 = lerp(corners[2], corners[3], weight(x));
    return lerp(v01, v23, weight(z));
}



void value(int x, int z, int layers) {
    float seed = 0;
    float seed_interval = 1.618; // Golden Ratio.
    for (int i = 0; i < layers; i++) {
        for (int c = 0; c < n_chunks; c++) {
            float noise = gen_value_noise(coordinates[(c * 2)],
                coordinates[(c * 2) + 1], x, z, seed);
            data[to_data_index(c, x, z)] += (noise - 0.5) / layers;   
        }
        seed += seed_interval;
    }
}


float gradient_value(int x, int z, float seed, float freq) {
    float hash = int(pow(x, 2) + (3 * x) + (2 * x * z) + z + pow(z, 2));
    hash /= 2.0f;
    return sin(freq * (hash + seed));
}


float[8] gen_gradient_vectors(int chunkX, int chunkZ, float seed, float freq) {
    float gradients[8];

    int offset = 0;
    for (int vert = 0; vert < 4; vert++) {
        int x_diff = offset & 1;
        int z_diff = (offset & 2) >> 1;

        // Values of the gradient.
        gradients[(vert * 2) + 0] = gradient_value(chunkX + x_diff, chunkZ + z_diff, seed, freq);
        gradients[(vert * 2) + 1] = gradient_value(chunkX + x_diff, chunkZ + z_diff, seed * 1.618, freq);
        offset++;
    }
    return gradients;
}


float calculate_perlin_noise(float gradients[8], int x, int z) {
    float true_width = chunkWidth - 1;

    // Generate offset vectors.
    vec2 off0 = vec2((0 - x) / true_width, (0 - z) / true_width);
    vec2 off1 = vec2((true_width - x) / true_width, (0 - z) / true_width);
    vec2 off2 = vec2((0 - x) / true_width, (true_width - z) / true_width);
    vec2 off3 = vec2((true_width - x) / true_width, (true_width - z) / true_width);

    // Calculate dot products.
    float dot0 = dot(off0, vec2(gradients[0], gradients[1]));
    float dot1 = dot(off1, vec2(gradients[2], gradients[3]));
    float dot2 = dot(off2, vec2(gradients[4], gradients[5]));
    float dot3 = dot(off3, vec2(gradients[6], gradients[7]));

    // Interpolate between values.
    float v01 = lerp(dot0, dot1, weight(x));
    float v23 = lerp(dot2, dot3, weight(x));
    return lerp(v01, v23, weight(z)); 
}




void perlin(int x, int z, int octaves, float freq_mult, float amp_mult) {
    for (int c = 0; c < n_chunks; c++) {
        int chunkX = coordinates[(c * 2)];
        int chunkZ = coordinates[(c * 2) + 1];

        float total_noise = 0;
        float seed = 0;
        float seed_interval = 1.618; // Golden Ratio.
        float oct_freq = 4;
        float oct_amp = 5;

        float amp_total = (oct_amp * (1 - pow(amp_mult, octaves))) / (1 - amp_mult);
        float amp_weight;

        // Calculate the octaves.
        for (int octave = 0; octave < octaves; octave++) {
            float gradients[8] = gen_gradient_vectors(chunkX, chunkZ, seed, oct_freq);

            amp_weight = oct_amp / amp_total;
            total_noise += calculate_perlin_noise(gradients, x, z) * amp_weight;

            // Update freq and amp and seed.
            oct_freq *= freq_mult;
            oct_amp *= amp_mult;
            seed += seed_interval;
        }

        data[to_data_index(c, x, z)] += total_noise;
    }
  
}


float random(float seed) {
    return sin(fract(seed * 123.43255) * 3.14);
}


// Noise generation starts.
void main() {
    int x = int(gl_GlobalInvocationID.x);
    int z = int(gl_GlobalInvocationID.z);

    for (int c = 0; c < n_chunks; c++) {
        data[to_data_index(c, x, z)] = 0;
    }
    

    if (x >= chunkWidth || z >= chunkWidth) {
        return;
    }

    switch (style) {
        case PERLIN:
            perlin(x, z, 4, 20.0f, 1.5f);
            break;
        case VALUE:
            value(x, z, 2);
            break;

        case LAYERED:
            value(x, z, 2);
            perlin(x, z, 4, 20.0f, 1.5f);

            for (int c = 0; c < n_chunks; c++) {
                data[to_data_index(c, x, z)] /= 2;
            }
            break;
    }
}