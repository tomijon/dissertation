#version 460 core

layout(local_size_x = 31, local_size_y = 1, local_size_z = 31) in;

// Uniforms.
uniform int n_chunks;
uniform int chunkWidth;
uniform int chunkHeight;

// Constant arrays used in surface generation.
const int n_vertices = 12;
const int n_permutations = 15;
const int n_corners = 8;

const float surf_vertices[3 * n_vertices] = float[3 * n_vertices](
    0, 0.5, 0,  0.5, 0, 0,  1, 0.5, 0,  0.5, 1, 0,
    0, 0.5, 1,  0.5, 0, 1,  1, 0.5, 1,  0.5, 1, 1,
    1, 0, 0.5,  1, 1, 0.5,  0, 1, 0.5,  0, 0, 0.5
);

const int permutations[n_permutations * n_vertices] = int[n_permutations * n_vertices](
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    4, 11, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    4, 11, 6, 11, 8, 6, 0, 0, 0, 0, 0, 0,
    4, 11, 5, 7, 6, 9, 0, 0, 0, 0, 0, 0,
    0, 5, 11, 5, 0, 6, 6, 0, 2, 0, 0, 0,
    0, 2, 4, 2, 6, 4, 0, 0, 0, 0, 0, 0,
    0, 5, 11, 5, 0, 6, 6, 0, 2, 4, 7, 10,
    4, 11, 5, 2, 3, 9, 0, 0, 0, 0, 0, 0,
    4, 11, 6, 11, 8, 6, 2, 3, 9, 0, 0, 0,
    4, 7, 10, 2, 3, 9, 5, 8, 6, 0, 0, 0,
    10, 11, 7, 7, 11, 5, 9, 1, 3, 8, 1, 9,
    0, 5, 4, 3, 9, 0, 5, 0, 9, 8, 5, 9,
    0, 10, 3, 4, 11, 5, 1, 2, 8, 7, 6, 9,
    10, 3, 4, 4, 3, 5, 3, 2, 5, 5, 2, 8,
    10, 5, 11, 5, 10, 2, 2, 10, 3, 6, 2, 5
);
const int permMap[n_permutations] = int[n_permutations](
    0, 8, 12, 9, 196, 204,
    198, 24, 28, 22, 90, 216,
    105, 232, 228);

// Rotating the active corners.
const int xRotations[n_corners] = int[n_corners](1, 5, 3, 7, 0, 4, 2, 6);
const int yRotations[n_corners] = int[n_corners](2, 3, 6, 7, 0, 1, 4, 5);
const int zRotations[n_corners] = int[n_corners](1, 3, 0, 2, 5, 7, 4, 6);

// Reserve rotations of the vertices.
const int reverseX[n_vertices] = int[n_vertices](2, 8, 6, 9, 0, 11, 4, 10, 5, 7, 3, 1);
const int reverseY[n_vertices] = int[n_vertices](10, 3, 9, 7, 11, 1, 8, 5, 2, 6, 4, 0);
const int reverseZ[n_vertices] = int[n_vertices](1, 2, 3, 0, 5, 6, 7, 4, 9, 10, 11, 8);


// Chunk Coordinates.
layout(std430, binding = 0) buffer Coordinates {
    int coordinates[];
};

// Data buffer.
layout(std430, binding = 1) buffer Data {
    float data[];
};

// Vertices output.
layout(std430, binding = 2) buffer Vertices {
    float vertices[];
};


// Convert the chunk number, x and z to an index in the data array.
int to_data_index(int chunk, int x, int z) {
    return (chunk * chunkWidth * chunkWidth) + (z * chunkWidth) + x;
}


// Finds if the corner at the x y z position (of vertices) is active.
// Returns 0 if inactive, 1 if active.
int is_active(int chunk, int x, int y, int z) {
    if (y == 0) return 1;
    float coordinate_height = (data[to_data_index(chunk, x, z)] + 0.5) * chunkHeight;
    if (y <= coordinate_height) {
        return 1;
    }
    return 0;
}


// Calculates the active corners of the current voxel.
int find_active(int chunk, int x, int y, int z) {
    int active_corners = 0;

    for (int cur_z = 0; cur_z < 2; cur_z++) {
        for (int cur_y = 0; cur_y < 2; cur_y++) {
            for (int cur_x = 0; cur_x < 2; cur_x++) {
                active_corners = (active_corners << 1) | is_active(
                    chunk, x + cur_x, y + cur_y, z + cur_z);
            }
        }
    }
    return active_corners;
}


// Determines if the active corners should be inverted (have its bits
// flipped).
bool can_invert(int active_corners) {
    int total = 0;
    for (int i = 0; i < 8; i++) {
        total += (active_corners >> i) & 1;
    }
    return total > 4;
}


void finish(int permutation, int c, int x, int y, int z, bool invert, int x_rot, int y_rot, int z_rot, int active_corners) {
    int current_perm[n_vertices];
    for (int i = 0; i < n_vertices; i++) {
        if (invert) {
            current_perm[i] = permutations[(permutation * n_vertices) + i];
        } else {
            current_perm[n_vertices - 1 - i] = permutations[(permutation * n_vertices) + i];
        }
    }

    // Undo Rotations.
    for (int nX = 0; nX < x_rot; nX++) {
        for (int vert = 0; vert < n_vertices; vert++){ 
            current_perm[vert] = reverseX[current_perm[vert]];
        }
    }

    for (int nY = 0; nY < y_rot; nY++) {
        for (int vert = 0; vert < n_vertices; vert++){ 
            current_perm[vert] = reverseY[current_perm[vert]];
        }
    }

    for (int nZ = 0; nZ < z_rot; nZ++) {
        for (int vert = 0; vert < n_vertices; vert++){ 
            current_perm[vert] = reverseZ[current_perm[vert]];
        }
    }
    
    // Vertex offset locations based on chunk width and position.
    int x_offset = ((chunkWidth - 1) * coordinates[c * 2]) + x;
    int y_offset = y;
    int z_offset = ((chunkWidth - 1) * coordinates[(c * 2) + 1]) + z;

    // Get vertices and normals.
    for (int vert = 0; vert < n_vertices; vert += 3) {
        // Collect 3 vertices that make up a surface.
        vec3 vertex1 = vec3(0);
        vec3 vertex2 = vec3(0);
        vec3 vertex3 = vec3(0);

        vertex1.x = surf_vertices[(current_perm[vert] * 3)] + x_offset;
        vertex1.y = surf_vertices[(current_perm[vert] * 3) + 1] + y_offset;
        vertex1.z = surf_vertices[(current_perm[vert] * 3) + 2] + z_offset;

        vertex2.x = surf_vertices[(current_perm[vert + 1] * 3)] + x_offset;
        vertex2.y = surf_vertices[(current_perm[vert + 1] * 3) + 1] + y_offset;
        vertex2.z = surf_vertices[(current_perm[vert + 1] * 3) + 2] + z_offset;

        vertex3.x = surf_vertices[(current_perm[vert + 2] * 3)] + x_offset;
        vertex3.y = surf_vertices[(current_perm[vert + 2] * 3) + 1] + y_offset;
        vertex3.z = surf_vertices[(current_perm[vert + 2] * 3) + 2] + z_offset;

        // Calculate normal of previous surface.
        vec3 u = vec3(0);
        vec3 v = vec3(0);

        u.x = vertex2.x - vertex1.x;
        u.y = vertex2.y - vertex1.y;
        u.z = vertex2.z - vertex1.z;

        v.x = vertex3.x - vertex1.x;
        v.y = vertex3.y - vertex1.y;
        v.z = vertex3.z - vertex1.z;

        vec3 normal = vec3(0);
        normal.x = (u.y * v.z) - (u.z * v.y);
        normal.y = (u.z * v.x) - (u.x * v.z);
        normal.z = (u.x * v.y) - (u.y * v.x);

        int width = chunkWidth - 1;
        int height = chunkHeight - 1;
        int surf = n_vertices * 2 * 3;
        int mesh = int(pow(width, 2)) * height * surf;

        int fetch = ((c * mesh) + (z + (y * width) + (x * width * height)) * surf);

        vertices[fetch + (vert * 6)] = vertex1.x;
        vertices[fetch + (vert * 6) + 1] = vertex1.y;
        vertices[fetch + (vert * 6) + 2] = vertex1.z;
        vertices[fetch + (vert * 6) + 3] = normal.x;
        vertices[fetch + (vert * 6) + 4] = normal.y;
        vertices[fetch + (vert * 6) + 5] = normal.z;

        vertices[fetch + (vert * 6) + 6] = vertex2.x;
        vertices[fetch + (vert * 6) + 7] = vertex2.y;
        vertices[fetch + (vert * 6) + 8] = vertex2.z;
        vertices[fetch + (vert * 6) + 9] = normal.x;
        vertices[fetch + (vert * 6) + 10] = normal.y;
        vertices[fetch + (vert * 6) + 11] = normal.z;

        vertices[fetch + (vert * 6) + 12] = vertex3.x;
        vertices[fetch + (vert * 6) + 13] = vertex3.y;
        vertices[fetch + (vert * 6) + 14] = vertex3.z;
        vertices[fetch + (vert * 6) + 15] = normal.x;
        vertices[fetch + (vert * 6) + 16] = normal.y;
        vertices[fetch + (vert * 6) + 17] = normal.z;
    }
}


void calculate(int c, int x, int y, int z) {
    int active_corners = find_active(c, x, y, z);

    // Inversion.
    bool invert = can_invert(active_corners);
    if (invert) {
        active_corners = (~active_corners) & 255;
    }

    // Rotate until permutation is found.
    for (int z_rot = 0; z_rot < 2; z_rot++) {
        int activeBeforeY = active_corners;
        for (int y_rot = 0; y_rot < 4; y_rot++) {
            int activeBeforeX = active_corners;
            for (int x_rot = 0; x_rot < 4; x_rot++) {
                for (int index = 0; index < n_permutations; index++) {
                    if (active_corners == permMap[index]) {
                        finish(index, c, x, y, z, invert, x_rot, y_rot, z_rot, active_corners);
                        return;
                    }
                }

                // Rotate on x axis.
                int new_active = 0;
                int lastBit = 0;

                // Bitwise magic.
                for (int i = 0; i < 8; i++) {
                    lastBit = (active_corners >> (7 - xRotations[i])) & 1;
                    new_active = (new_active << 1) + lastBit;
                }
                active_corners = new_active;
            }

            // Rotate on y axis.
            int new_active = 0;
            int lastBit = 0;

            // Bitwise magic.
            for (int i = 0; i < 8; i++) {
                lastBit = (activeBeforeX >> (7 - yRotations[i])) & 1;
                new_active = (new_active << 1) + lastBit;
            }
            active_corners = new_active;
        }

        // Rotate on z axis.
        int new_active = 0;
        int lastBit = 0;

        // Bitwise magic.
        for (int i = 0; i < 8; i++) {
            lastBit = (activeBeforeY >> (7 - zRotations[i])) & 1;
            new_active = (new_active << 1) + lastBit;
        }
        active_corners = new_active;
    }
}


// Start marching.
void march() {
    int x = int(gl_GlobalInvocationID.x);
    int z = int(gl_GlobalInvocationID.z);

    // Safeguard.
    if (x >= chunkWidth - 1 || z >= chunkWidth - 1) {
        return;
    }

    // Marching the column.
    for (int c = 0; c < n_chunks; c++) {
        for (int y = 0; y < chunkHeight - 1; y++) {
             calculate(c, x, y, z);
        }
    }
}

void main() {
    march();
}