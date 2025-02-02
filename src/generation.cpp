#include <cmath>
#include <chrono>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

#include <fstream>

#include <glad/glad.h>
#include <glm/glm.hpp>

#include <generation.hpp>
#include <camera.hpp>

using std::cout;

namespace Generation {
    /** The 12 vertex positions that will be used to generate the
     *  surfaces in the marching cubes algorithm.
     */
    const GLfloat *vertices = new float[36] {
        1, 0.5, 0,  0.5, 0, 0,  0, 0.5, 0,  0.5, 1, 0,
        1, 0.5, 1,  0.5, 0, 1,  0, 0.5, 1,  0.5, 1, 1,
        0, 0, 0.5,  0, 1, 0.5,  1, 1, 0.5,  1, 0, 0.5
    };

    const Permutation *permutations = new Permutation[15] {
        {0, new GLuint[0]{}},
        {3, new GLuint[3]{4, 11, 5}},
        {6, new GLuint[6]{4, 11, 6, 11, 8, 6}},
        {6, new GLuint[6]{4, 11, 5, 7, 6, 9}},
        {9, new GLuint[9]{0, 5, 11, 5, 0, 6, 6, 0, 2}},
        {6, new GLuint[6]{0, 2, 4, 2, 6, 4}},
        {12, new GLuint[12]{0, 5, 11, 5, 0, 6, 6, 0, 2, 4, 7, 10}},
        {6, new GLuint[6]{4, 11, 5, 2, 3, 9}},
        {9, new GLuint[9]{4, 11, 6, 11, 8, 6, 2, 3, 9}},
        {9, new GLuint[9]{4, 7, 10, 2, 3, 9, 5, 8, 6}},
        {12, new GLuint[12]{10, 11, 7, 7, 11, 5, 9, 1, 3, 8, 1, 9}},
        {12, new GLuint[12]{0, 5, 4, 3, 9, 0, 5, 0, 9, 8, 5, 9}},
        {12, new GLuint[12]{0, 10, 3, 4, 11, 5, 1, 2, 8, 7, 6, 9}},
        {12, new GLuint[12]{10, 3, 4, 4, 3, 5, 3, 2, 5, 5, 2, 8}},
        {12, new GLuint[12]{10, 5, 11, 5, 10, 2, 2, 10, 3, 6, 2, 5}}
    };

    const int *permMap = new int[15]{
        0, 8, 12, 9, 196,
        204, 198, 24, 28, 22,
        90, 216, 105, 232, 228
    };

    const int *xRotations = new int[8]{1, 5, 3, 7, 0, 4, 2, 6};
    const int *yRotations = new int[8]{2, 3, 6, 7, 0, 1, 4, 5};
    const int *zRotations = new int[8]{1, 3, 0, 2, 5, 7, 4, 6};
    const int *reverseX = new int[12]{2, 8, 6, 9, 0, 11, 4, 10, 5, 7, 3, 1};
    const int *reverseY = new int[12]{10, 3, 9, 7, 11, 1, 8, 5, 2, 6, 4, 0};
    const int *reverseZ = new int[12]{1, 2, 3, 0, 5, 6, 7, 4, 9, 10, 11, 8};


    /** Chunk manager constructor.
     * 
     *  Camera view diameter is used to initialise the 2D array of
     *  chunks that holds chunks that are currently rendered.
     */
    Manager::Manager(NoiseType style, int width, int height, Camera::View *camera,
            GenMode mode, GLuint generate, GLuint march, GLuint renderShader)
        : style(style), chunkWidth(width), chunkHeight(height), camera(camera),
            mode(mode), generate(generate), march(march), renderShader(renderShader) {

        // Generate rendered chunk array.
        chunks = new Chunk**[camera->view_diameter()];
        for (int i = 0; i < camera->view_diameter(); i++) {
            chunks[i] = new Chunk*[camera->view_diameter()];
        }

        // Populate the array.
        std::vector<int> position = camera_position();
        for (int x = 0; x < camera->view_diameter(); x++) {
            for (int z = 0; z < camera->view_diameter(); z++) {
                chunks[z][x] = generate_chunk_CPU(
                    position[0] - camera->view_range() + x,
                    position[1] - camera->view_range() + z
                );
            }
        }
        
        last_position[0] = position[0];
        last_position[1] = position[1];
    }


    /** Update method for updating the chunk array.
     * 
     *  TODO:   This function needs nerfed its too big. There must be a
     *          shorter way. Just gotta figure out how to trim the shift
     *          section as thats the more complicated bit.
     */
    void Manager::update() {
        std::vector<int> new_position = camera_position();
        // x axis shift.
        if (new_position[0] != last_position[0]) {
            int dx = new_position[0] - last_position[0];

            // Moving too fast.
            if (std::abs(dx) > 1) {
                throw std::runtime_error("Moving too fast!\n");
            }

            // Camera moved left so shift everything right.
            if (dx == -1) {
                // Delete all the chunks furthest right because they are
                // about to be sent to the shadow realm.
                for (int z = 0; z < camera->view_diameter(); z++) {
                    delete_chunk(chunks[z][camera->view_diameter() - 1]);
                }
                // Shift everything right.
                for (int z = 0; z < camera->view_diameter(); z++) {
                    for (int x = camera->view_diameter() - 2; x >= 0; x--) {
                        chunks[z][x + 1] = chunks[z][x];
                    }
                }
                // Fill in the gaps.
                if (mode == CPU) {
                    for (int z = 0; z < camera->view_diameter(); z++) {
                        Chunk *old = chunks[z][1];
                        Chunk *new_chunk = generate_chunk_CPU(old->x - 1, old->z);
                        chunks[z][0] = new_chunk;
                    }
                } else if (mode == GPU) {
                    int n_chunks = camera->view_diameter();
                    int coordinates[n_chunks][2];

                    // Fetch coordinates of all chunks.
                    for (int z = 0; z < camera->view_diameter(); z++) {
                        Chunk *old = chunks[z][1];
                        coordinates[z][0] = old->x - 1;
                        coordinates[z][1] = old->z;
                    }

                    // Generate and insert chunks.
                    Chunk **generated = generate_chunks_GPU(n_chunks, coordinates);
                    for (int z = 0; z < camera->view_diameter(); z++) {
                        chunks[z][0] = generated[z];
                    }
                }

            // Otherwise the camera moved right so shift everything left.
            } else {
                // Delete the far left chunks.
                for (int z = 0; z < camera->view_diameter(); z++) {
                    delete_chunk(chunks[z][0]);
                }
                // Shift everything left.
                for (int z = 0; z < camera->view_diameter(); z++) {
                    for (int x = 1; x < camera->view_diameter(); x++) {
                        chunks[z][x - 1] = chunks[z][x];
                    }
                }
                // Fill in the gaps.
                if (mode == CPU) {
                    for (int z = 0; z < camera->view_diameter(); z++) {
                        Chunk *old = chunks[z][camera->view_diameter() - 2];
                        Chunk *new_chunk = generate_chunk_CPU(old->x + 1, old->z);
                        chunks[z][camera->view_diameter() - 1] = new_chunk;
                    }
                } else if (mode == GPU) {
                    int n_chunks = camera->view_diameter();
                    int coordinates[n_chunks][2];

                    // Fetch coordinates of all chunks.
                    for (int z = 0; z < camera->view_diameter(); z++) {
                        Chunk *old = chunks[z][camera->view_diameter() - 2];
                        coordinates[z][0] = old->x + 1;
                        coordinates[z][1] = old->z;
                    }

                    // Generate and insert chunks.
                    Chunk **generated = generate_chunks_GPU(n_chunks, coordinates);
                    for (int z = 0; z < camera->view_diameter(); z++) {
                        chunks[z][camera->view_diameter() - 1] = generated[z];
                    }
                }
            }
            last_position[0] = new_position[0];
            cout << "New Position: " << new_position[0] << ", " << new_position[1] << std::endl;
        }

        // z axis shift.
        if (new_position[1] != last_position[1]) {
            int dz = new_position[1] - last_position[1];

            // Moving too fast.
            if (std::abs(dz) > 1) {
                throw std::runtime_error("Moving too fast!\n");
            }

            // Camera moved up so shift everything down.
            if (dz == -1) {
                // Delete all the chunks at the bottom because they are
                // about to be sent to the shadow realm.
                for (int x = 0; x < camera->view_diameter(); x++) {
                    delete_chunk(chunks[camera->view_diameter() - 1][x]);
                }
                // Shift everything down.
                for (int x = 0; x < camera->view_diameter(); x++) {
                    for (int z = camera->view_diameter() - 2; z >= 0; z--) {
                        chunks[z + 1][x] = chunks[z][x];
                    }
                }
                // Fill in the gaps.
                if (mode == CPU) {
                    for (int x = 0; x < camera->view_diameter(); x++) {
                        Chunk *old = chunks[1][x];
                        Chunk *new_chunk = generate_chunk_CPU(old->x, old->z - 1);
                        chunks[0][x] = new_chunk;
                    }
                } else if (mode == GPU) {
                    int n_chunks = camera->view_diameter();
                    int coordinates[n_chunks][2];

                    // Fetch coordinates of all chunks.
                    for (int x = 0; x < camera->view_diameter(); x++) {
                        Chunk *old = chunks[1][x];
                        coordinates[x][0] = old->x;
                        coordinates[x][1] = old->z - 1;
                    }

                    // Generate and insert chunks.
                    Chunk **generated = generate_chunks_GPU(n_chunks, coordinates);
                    for (int x = 0; x < camera->view_diameter(); x++) {
                        chunks[0][x] = generated[x];
                    }
                }

            // Otherwise the camera moved down so shift everything up.
            } else {
                // Delete the top chunks.
                for (int x = 0; x < camera->view_diameter(); x++) {
                    delete_chunk(chunks[0][x]);
                }
                // Shift everything up.
                for (int x = 0; x < camera->view_diameter(); x++) {
                    for (int z = 1; z < camera->view_diameter(); z++) {
                        chunks[z - 1][x] = chunks[z][x];
                    }
                }
                // Fill in the gaps.
                if (mode == CPU) {
                    for (int x = 0; x < camera->view_diameter(); x++) {
                        Chunk *old = chunks[camera->view_diameter() - 2][x];
                        Chunk *new_chunk = generate_chunk_CPU(old->x, old->z + 1);
                        chunks[camera->view_diameter() - 1][x] = new_chunk;
                    }
                } else if (mode == GPU) {
                    int n_chunks = camera->view_diameter();
                    int coordinates[n_chunks][2];

                    // Fetch coordinates of all chunks.
                    for (int x = 0; x < camera->view_diameter(); x++) {
                        Chunk *old = chunks[camera->view_diameter() - 2][x];
                        coordinates[x][0] = old->x;
                        coordinates[x][1] = old->z + 1;
                    }

                    // Generate and insert chunks.
                    Chunk **generated = generate_chunks_GPU(n_chunks, coordinates);
                    for (int x = 0; x < camera->view_diameter(); x++) {
                        chunks[camera->view_diameter() - 1][x] = generated[x];
                    }
                }
            }
            last_position[1] = new_position[1];
            cout << "New Position: " << new_position[0] << ", " << new_position[1] << std::endl;
        }
    }


    /** Performs the necessary operations to generate a chunk
     *  on the CPU. The chunk is dynamically allocated when marched.
     */
    Chunk * Manager::generate_chunk_CPU(int x, int z) {
        auto start = std::chrono::high_resolution_clock::now();
        

        int **data = style == PERLIN ? perlin(x, z, 4, 10, 1.5) : value(x, z, 4);

        Chunk *chunk = new Chunk;
        chunk->x = x;
        chunk->z = z;

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> duration = end - start;
        std::cout << "CPU Generation Time: " << std::fixed << std::setprecision(3) << duration.count() << "ms\n";

        start = std::chrono::high_resolution_clock::now();
        march_cubes_CPU(chunk, data);

        // Delete the data.
        for (int y = 0; y < chunkHeight; y++) {
            delete[] data[y];
        }
        delete[] data;

        end = std::chrono::high_resolution_clock::now();
        duration = end - start;
        std::cout << "CPU March Time: " << std::fixed << std::setprecision(3) << duration.count() << "ms\n";
        return chunk;
    }


    /** Translates the camera's world position into chunk based
     * coordinates.
     * 
     * Returns a 2D array representing the x and z chunk position of the
     * camera.
     */
    std::vector<int> Manager::camera_position() {
        std::vector<int> position(2);
        glm::vec4 cam_pos = camera->get_position();

        position[0] = std::floor(cam_pos.x / ((float)chunkWidth - 1));
        position[1] = std::floor(cam_pos.z / ((float)chunkWidth - 1));
        return position;
    }


    /** Frees memory allocated by the chunk, including the vertex object
     *  buffers.
     * 
     *  WARNING: ASSUMES THE BUFFERS ARE ALREADY UNBOUND.
     */
    void Manager::delete_chunk(Chunk *chunk) {
        // Freeing the GPU buffers first.
        glDeleteBuffers(1, &chunk->vertices);
        glDeleteVertexArrays(1, &chunk->vertexArrayObject);
        delete chunk;
    }


    /** Calulates the weight that should be used inside the lerp
     *  function.
     * 
     *  Specially made for extra terrain smoothness :ok_hand:
     */
    float Manager::weight(float value) {
        int range = (int)(value / ((chunkWidth - 1) / 2.0f));
        int cos_mult = -(int)(((range % 2) - 0.5f) * 2);
        int addition = (2 * range) + 1;

        float scalar = (2 * pi) / (chunkWidth - 1);
        return ((cos_mult * std::cos((scalar * value) + pi)) + addition) / 4;
    }


    /** Generate a value for the value noise generation method.
     */
    float Manager::value_noise(int x, int z, float seed) {
        float hash = std::pow(x, 2) + (3 * x) + (2 * x * z) + z + std::pow(z, 2);
        hash /= 2.0f;
        float cos_section = std::cos(hash + seed);
        float sin_section = std::sin(hash + seed);
        return (-1 * ((std::pow(4, cos_section) + (6 * std::pow(sin_section, 2))) / 4.0f) + 2) * 0.5;
    }


    /** Generate the random values for each of the 4 corners of the
     *  chunk.
     */
    void Manager::gen_vertex_noise(float corners[4], int chunkX, int chunkZ, float seed) {
        int offset = 0;
        for (int vert = 0; vert < 4; vert++) {
            int x_diff = offset & 1;
            int z_diff = (offset & 2) >> 1;
            corners[vert] = value_noise(chunkX + x_diff, chunkZ + z_diff, seed);
            offset++;
        }
    }


    /** Calculates the noise value at the current x, z coordinate in the
     *  chunk.
     * 
     *  Performs linear interpolation between the 4 corners of the
     *  chunk.
     */
    float Manager::gen_value_noise(int chunkX, int chunkZ, int x, int z, float seed) {
        float corners[4];
        gen_vertex_noise(corners, chunkX, chunkZ, seed);

        // Lerp
        float v01 = lerp(corners[0], corners[1], weight(x));
        float v23 = lerp(corners[2], corners[3], weight(x));
        return lerp(v01, v23, weight(z));
    }


    /** Creates the heightmap data for the chunk using value noise.
     * 
     *  This function dynamically allocates memory for the returned
     *  value.
     */
    int ** Manager::value(int chunkX, int chunkZ, int layers) {
        int **data = new int*[chunkHeight]{};
        for (int y = 0; y < chunkHeight; y++) {
            data[y] = new int[chunkWidth]{};
        }

        // Loop through the 32x32 chunk.
        for (int x = 0; x < chunkWidth; x++) {
            for (int z = 0; z < chunkWidth; z++) {
                float total_noise = 0;
                float seed = 0;
                float seed_interval = 1.618;

                // Calculate the layers.
                for (int layer = 0; layer < layers; layer++) {
                    total_noise += gen_value_noise(chunkX, chunkZ, x, z, seed);
                    seed += seed_interval;
                }
                total_noise /= layers;

                // Insert Noise.
                for (int y = 0; y < (int) ((chunkHeight / 2) * total_noise) + (chunkHeight / 2); y++) {
                    data[y][z] += 1 << x;
                }
            }
        }
        return data;
    }


    /** Generates a value for one of the ordinates in the gradient
     *  vector.
     */
    float Manager::gradient_value(int x, int z, float seed, float freq) {
        float hash = std::pow(x, 2) + (3 * x) + (2 * x * z) + z + std::pow(z, 2);
        hash /= 2.0f;
        return sin(freq * (hash + seed));
    }


    /** Generates all the gradient vectors (2D gradient vectors for each
     *  of the 4 corners).
     * 
     *  Performs this on an already created array because dynamically
     *  allocating it 32 x 32 times is gonna hurt.
     */
    void Manager::gen_gradient_vectors(float gradients[4][2], int chunkX, int chunkZ, float seed, float freq) {
        int offset = 0;
        for (int vert = 0; vert < 4; vert++) {
            int x_diff = offset & 1;
            int z_diff = (offset & 2) >> 1;

            // Values of the gradient.
            gradients[vert][0] = gradient_value(chunkX + x_diff, chunkZ + z_diff, seed, freq);
            gradients[vert][1] = gradient_value(chunkX + x_diff, chunkZ + z_diff, seed * 1.618, freq);
            offset++;
        }
    }


    /** Calculates the current height value at the x, z coordinate.
     * 
     *  Performs linear interpolation between the dot products of the
     *  4 corners.
     */
    float Manager::calculate_perlin_noise(float gradients[4][2], int x, int z) {
        float true_width = chunkWidth - 1;

        // Generate offset vectors.
        float off0[2] = {(0 - x) / true_width, (0 - z) / true_width};
        float off1[2] = {(true_width - x) / true_width, (0 - z) / true_width};
        float off2[2] = {(0 - x) / true_width, (true_width - z) / true_width};
        float off3[2] = {(true_width - x) / true_width, (true_width - z) / true_width};

        // Calculate dot products.
        float dot0 = dot(off0, gradients[0]);
        float dot1 = dot(off1, gradients[1]);
        float dot2 = dot(off2, gradients[2]);
        float dot3 = dot(off3, gradients[3]);

        // Interpolate between values.
        float v01 = lerp(dot0, dot1, weight(x));
        float v23 = lerp(dot2, dot3, weight(x));
        return lerp(v01, v23, weight(z)); 
    }


    /** Generates the heightmap for the chunk using perlin noise.
     * 
     *  Octaves - How many times to apply the noise.
     *  freq - How much to multiply the initial frequency by each octave
     *  amp - How much to multiple the initial amplitude by each octave.
     * 
     *  This function dynamically allocates memory for the returned
     *  value.
     */ 
    int ** Manager::perlin(int chunkX, int chunkZ, int octaves, float freq_mult, float amp_mult) {
        int **data = new int*[chunkHeight]{};
        for (int y = 0; y < chunkHeight; y++) {
            data[y] = new int[chunkWidth]{};
        }

        // Each section in the chunk.
        for (int x = 0; x < chunkWidth; x++) {
            for (int z = 0; z < chunkWidth; z++) {
                float total_noise = 0;
                float seed = 0;
                float seed_interval = 1.618; // Golden Ratio.
                float oct_freq = 4;
                float oct_amp = 5;

                float amp_total = (oct_amp * (1 - std::pow(amp_mult, octaves))) / (1 - amp_mult);
                float amp_weight;

                // Calculate the octaves.
                for (int octave = 0; octave < octaves; octave++) {
                    float gradients[4][2];
                    gen_gradient_vectors(gradients, chunkX, chunkZ, seed, oct_freq);

                    amp_weight = oct_amp / amp_total;
                    total_noise += calculate_perlin_noise(gradients, x, z) * amp_weight;

                    // Update freq and amp and seed.
                    oct_freq *= freq_mult;
                    oct_amp *= amp_mult;
                    seed += seed_interval;
                }

                // Insert Noise.
                for (int y = 0; y < (int) ((chunkHeight / 2) * total_noise) + (chunkHeight / 2); y++) {
                    data[y][z] += 1 << x;
                }
            }
        }
        return data;
    }


    /** Perform the marching cubes algorithm on the CPU.
     */
    void Manager::march_cubes_CPU(Chunk *chunk, int **data) {
        // Vertices 0, 1, 2 and Normals 3, 4, 5
        std::vector<Rotations> *permsToReverse = new std::vector<Rotations>();
        uint8_t active = 0; // Gets overwritten each iteration.

        // Loop through data.
        for (int y = 0; y < chunkHeight - 1; y++) {
            for (int z = 0; z < chunkWidth - 1; z++) {
                for (int x = 0; x < chunkWidth - 1; x++) {
                    // Generate number that will map to permutation.
                    active = (active << 2) + (mask & (data[y][z] >> x));
                    active = (active << 2) + (mask & (data[y + 1][z] >> x));
                    active = (active << 2) + (mask & (data[y][z + 1] >> x));
                    active = (active << 2) + (mask & (data[y + 1][z + 1] >> x));
                    
                    permsToReverse->push_back(getMatch((int) active, x, y, z));
                }
            }
        }

        // Reverse the permutations.
        int size = permsToReverse->size();
        Rotations *reversed = new Rotations[size];
        for (int n_perm = 0; n_perm < size; n_perm++) {
            reversed[n_perm] = reverse((*permsToReverse)[n_perm]);
        }
        delete permsToReverse;
        buildMesh(chunk, size, reversed);
    }


    /** Finds a matching permutation.
     */
    Rotations Manager::getMatch(int active, int x, int y, int z) {
        // Count and flip active corners more than 4.
        int total = 0;
        bool invert = false;
        for (int i = 0; i < 8; i++) {
            total += (active >> i) & 1;
        }
        if (total > 4) {
            active = (~active) & 255;
            invert = true;
        } 

        Permutation *permutation;

        for (int zRot = 0; zRot < 2; zRot++) {
            int activeBeforeY = active;
            for (int yRot = 0; yRot < 4; yRot++) {
                int activeBeforeX = active;
                for (int xRot = 0; xRot < 4; xRot++) {
                    // Perform test.
                    for (int index = 0; index < 15; index++) {
                        if (active == permMap[index]) {
                            Rotations rotation = {
                                &(permutations[index]),
                                xRot, yRot, zRot,
                                x, y, z,
                                invert
                            };
                            return rotation;
                        }
                    }
                    active = rotateActive(active, xRotations);
                }
                active = rotateActive(activeBeforeX, yRotations);
            }
            active = rotateActive(activeBeforeY, zRotations);
        }
        return Rotations{};
    }


    /** Rotates the active vertices around the axis determined by
     *  rotations.
     */
    int Manager::rotateActive(int active, const int *rotations) {
        int newActive = 0;
        int lastBit = 0;

        // Bitwise magic.
        for (int i = 0; i < 8; i++) {
            lastBit = (active >> (7 - rotations[i])) & 1;
            newActive = (newActive << 1) + lastBit;
        }
        return newActive;
    }


    /** Reverses the rotations performed to achieve the target
     *  permutation, on a copy of the permutation.
     */
    Rotations Manager::reverse(Rotations rotation) {
        // Generate copy.
        Permutation *copy = new Permutation{
            rotation.permutation->n_vertices,
            new GLuint[rotation.permutation->n_vertices]{},
        };
        // std::memcpy(copy->vertices, rotation.permutation->vertices, sizeof(GLuint) * copy->n_vertices);

        // Reverse if necessary.
        if (rotation.invert) {
            for (int vert = 0; vert < copy->n_vertices; vert++) {
                copy->vertices[vert] = rotation.permutation->vertices[copy->n_vertices - 1 - vert];
            }
        } else {
            std::memcpy(copy->vertices, rotation.permutation->vertices, sizeof(GLuint) * copy->n_vertices);
        }

        // Reverse rotations in reverse order.
        for (int nX = 0; nX < rotation.x_rotations; nX++) {
            for (int vert = 0; vert < copy->n_vertices; vert++) {
                copy->vertices[vert] = reverseX[copy->vertices[vert]];
            }
        }

        for (int nY = 0; nY < rotation.y_rotations; nY++) {
            for (int vert = 0; vert < copy->n_vertices; vert++) {
                copy->vertices[vert] = reverseY[copy->vertices[vert]];
            }
        }                        

        for (int nZ = 0; nZ < rotation.z_rotations; nZ++) {
            for (int vert = 0; vert < copy->n_vertices; vert++) {
                copy->vertices[vert] = reverseZ[copy->vertices[vert]];
            }
        }
        rotation.permutation = copy;
        return rotation;
    }


    /** Builds a mesh from an array of permutations and stores the
     *  vertex information in the chunk.
     */
    void Manager::buildMesh(Chunk *chunk, int totalPerms, Rotations *rotations) {
        std::vector<GLfloat> *mesh = new std::vector<GLfloat>();

        for (int n_perm = 0; n_perm < totalPerms; n_perm++) {
            const Permutation *permutation = rotations[n_perm].permutation;

            int xOffset = ((chunkWidth - 1) * chunk->x) + rotations[n_perm].x;
            int yOffset = rotations[n_perm].y;
            int zOffset = ((chunkWidth - 1) * chunk->z) + rotations[n_perm].z;

            float surfaceVertices[9];
            float *vertex1 = &surfaceVertices[0];
            float *vertex2 = &surfaceVertices[3];
            float *vertex3 = &surfaceVertices[6];

            for (int vert = 0; vert < permutation->n_vertices; vert += 3) {
                surfaceVertices[0] = vertices[(permutation->vertices[vert] * 3)] + xOffset;
                surfaceVertices[1] = vertices[(permutation->vertices[vert] * 3) + 1] + yOffset;
                surfaceVertices[2] = vertices[(permutation->vertices[vert] * 3) + 2] + zOffset;
                surfaceVertices[3] = vertices[(permutation->vertices[vert + 1] * 3)] + xOffset;
                surfaceVertices[4] = vertices[(permutation->vertices[vert + 1] * 3) + 1] + yOffset;
                surfaceVertices[5] = vertices[(permutation->vertices[vert + 1] * 3) + 2] + zOffset;
                surfaceVertices[6] = vertices[(permutation->vertices[vert + 2] * 3)] + xOffset;
                surfaceVertices[7] = vertices[(permutation->vertices[vert + 2] * 3) + 1] + yOffset;
                surfaceVertices[8] = vertices[(permutation->vertices[vert + 2] * 3) + 2] + zOffset;

                // Calculate normal.
                float u[3];
                float v[3];

                u[0] = vertex2[0] - vertex1[0];
                u[1] = vertex2[1] - vertex1[1];
                u[2] = vertex2[2] - vertex1[2];

                v[0] = vertex3[0] - vertex1[0];
                v[1] = vertex3[1] - vertex1[1];
                v[2] = vertex3[2] - vertex1[2];

                float normal[3];
                normal[0] = (u[1] * v[2]) - (u[2] * v[1]);
                normal[1] = (u[2] * v[0]) - (u[0] * v[2]);
                normal[2] = (u[0] * v[1]) - (u[1] * v[0]);

                // Add to mesh.
                for (int surfVert = 0; surfVert < 9; surfVert += 3) {
                    for (int ord = 0; ord < 3; ord++) {
                        mesh->push_back(surfaceVertices[surfVert + ord]);
                    }

                    for (int ord = 0; ord < 3; ord++) {
                        mesh->push_back(normal[ord]);
                    }
                }
            }
            delete[] permutation->vertices;
            delete permutation;
        }

        // Generate buffers.
        glGenVertexArrays(1, &(chunk->vertexArrayObject));
        glCreateBuffers(1, &(chunk->vertices));
        glBindVertexArray(chunk->vertexArrayObject);
        glBindBuffer(GL_ARRAY_BUFFER, chunk->vertices);
        glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * mesh->size(), mesh->data(), GL_DYNAMIC_DRAW);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), NULL);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_TRUE, 6 * sizeof(GLfloat), (void*)(3 * sizeof(GLfloat)));
        glEnableVertexAttribArray(1);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindVertexArray(0);

        chunk->n_faces = mesh->size() / 2;
        delete mesh;
        delete rotations;
    }


    /** Generates a number of chunks on the GPU.
     * 
     *  n_chunks - will be sent as a uniform to the generate shader.
     *  coordinates - the x, z location of each chunk that is to be
     *                  generated.
     * 
     *  This function calls glUseProgram. It does not reassign the
     *  previous program on cleanup.
     */
    Chunk** Manager::generate_chunks_GPU(int n_chunks, int coordinates[][2]) {
        auto start = std::chrono::high_resolution_clock::now();
        
        glUseProgram(generate);
        glUniform1i(glGetUniformLocation(generate, "n_chunks"), n_chunks);
        glUniform1i(glGetUniformLocation(generate, "chunkWidth"), chunkWidth);
        glUniform1i(glGetUniformLocation(generate, "chunkHeight"), chunkHeight);
        glUniform1i(glGetUniformLocation(generate, "style"), (int) style);

        // Shader Storage Buffers.
        GLuint chunk_positions;
        GLuint data;

        int converted[n_chunks * 2];
        for (int i = 0; i < n_chunks; i++) {
            converted[(i * 2)] = coordinates[i][0];
            converted[(i * 2) + 1] = coordinates[i][1];
        }

        // Copying coordinates into the shader buffer.
        glGenBuffers(1, &chunk_positions);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, chunk_positions);
        glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(int) * 2 * n_chunks,
            converted, GL_STREAM_DRAW);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, chunk_positions);

        // Creating the data buffer to store height values.
        glGenBuffers(1, &data);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, data);
        glBufferData(GL_SHADER_STORAGE_BUFFER,
            sizeof(float) * n_chunks * chunkWidth * chunkWidth,
            nullptr, GL_STREAM_DRAW);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, data);

        // Initilise chunk objects.
        Chunk **chunks = new Chunk*[n_chunks]{};
        for (int c = 0; c < n_chunks; c++) {
            chunks[c] = new Chunk;
            chunks[c]->x = coordinates[c][0];
            chunks[c]->z = coordinates[c][1];
            chunks[c]->n_faces = (12 * 3 * 2 * (chunkWidth - 1) * (chunkWidth - 1) * (chunkHeight - 1));
        }

        // Generate.
        GLuint fence;
        GLsync sync = glFenceSync(GL_SYNC_GPU_COMMANDS_COMPLETE, 0);
        glDispatchCompute(chunkWidth, 1, chunkWidth);
        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

        // Timing for generation.
        GLenum waitReturn = GL_UNSIGNALED;
        while (waitReturn != GL_ALREADY_SIGNALED && waitReturn != GL_CONDITION_SATISFIED) {
            waitReturn = glClientWaitSync(sync, GL_SYNC_FLUSH_COMMANDS_BIT, 1);
        }
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> duration = end - start;
        std::cout << "Total Chunks: " << n_chunks << std::endl;
        std::cout << "GPU Generation Time: " << std::fixed << std::setprecision(3) << duration.count() << "ms\n";

        // Marching cubes with timing.
        start = std::chrono::high_resolution_clock::now();
        march_cubes_GPU(chunks, n_chunks, chunk_positions, data);

        waitReturn = GL_UNSIGNALED;
        while (waitReturn != GL_ALREADY_SIGNALED && waitReturn != GL_CONDITION_SATISFIED) {
            waitReturn = glClientWaitSync(sync, GL_SYNC_FLUSH_COMMANDS_BIT, 1);
        }
        end = std::chrono::high_resolution_clock::now();
        duration = end - start;
        std::cout << "GPU March Time: " << std::fixed << std::setprecision(3) << duration.count() << "ms\n";
        

        glDeleteSync(sync);
        glDeleteBuffers(1, &chunk_positions);
        glDeleteBuffers(1, &data);        
        return chunks;
    }

    void Manager::march_cubes_GPU(Chunk **chunks, int n_chunks,
        GLuint chunk_positions, GLuint data)
    {
        glUseProgram(march);
        glUniform1i(glGetUniformLocation(march, "n_chunks"), n_chunks);
        glUniform1i(glGetUniformLocation(march, "chunkWidth"), chunkWidth);
        glUniform1i(glGetUniformLocation(march, "chunkHeight"), chunkHeight);
        int mesh_size = sizeof(float) * 12 * 2 * 3 * std::pow(chunkWidth - 1, 2) * (chunkHeight - 1);

        // Generating vertex array.
        GLuint vertices;
        glGenBuffers(1, &vertices);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, vertices);
        glBufferData(GL_SHADER_STORAGE_BUFFER, mesh_size * n_chunks,
            nullptr, GL_STREAM_DRAW);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, vertices);

        // Dispatch
        glDispatchCompute(chunkWidth - 1, 1, chunkWidth - 1);
        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

        //glUseProgram(renderShader);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

        // Create VBO and copy vertices across for each chunk.
        for (int c = 0; c < n_chunks; c++) {
            // std::vector<GLfloat> *mesh = new std::vector<GLfloat>();
            Chunk *current = chunks[c];

            // Initialise vertex array object.
            glGenVertexArrays(1, &(current->vertexArrayObject));
            glBindVertexArray(current->vertexArrayObject);

            // Initialise vertices buffer.
            glGenBuffers(1, &(current->vertices));
            glBindBuffer(GL_ARRAY_BUFFER, current->vertices);
            glBufferData(GL_ARRAY_BUFFER, mesh_size, nullptr, GL_DYNAMIC_DRAW);

            // Copy over the data.
            int read_offset = mesh_size * c;
            glBindBuffer(GL_SHADER_STORAGE_BUFFER, vertices);
            glCopyBufferSubData(GL_SHADER_STORAGE_BUFFER, GL_ARRAY_BUFFER,
                read_offset, 0, mesh_size);

            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), NULL);
            glEnableVertexAttribArray(0);
            glVertexAttribPointer(1, 3, GL_FLOAT, GL_TRUE, 6 * sizeof(GLfloat), (void*)(3 * sizeof(GLfloat)));
            glEnableVertexAttribArray(1);

            glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
            glBindBuffer(GL_ARRAY_BUFFER, 0);
            glBindVertexArray(0);

            // Output vertices.
            glBindBuffer(GL_ARRAY_BUFFER, current->vertices);
            // delete mesh;
        }
        glDeleteBuffers(1, &vertices);
    }
}