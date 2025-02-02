#ifndef GENERATION_H
#define GENERATION_H

#include <vector>
#include <camera.hpp>

namespace Generation {
    // Forward Type Declarations
    typedef float GLfloat;
    typedef unsigned int GLuint;

    /** A structure that stores information about 1 of the 15 different
     *  base marching cubes permutations.
     */
    typedef struct Permutation {
        int n_vertices;
        GLuint *vertices;
    } Permutation;


    /** A structure for storing rotation information about a voxel.
     *  Stores the x, y and z rotations that were required to get to the
     *  permutation.
     * 
     *  Also stores the x, y, and z coordinate and a boolean that
     *  determines if the order of the vertices needs to be flipped in
     *  order to invert the winding order.
     */
    typedef struct Rotations {
        const Permutation *permutation;
        int x_rotations, y_rotations, z_rotations;
        int x, y, z;
        bool invert;
    } Rotations;


    /** Struct that contains the vertex objects for chunk rendering.
     *  Also contains the x and z position of the chunk in the world.
     */
    typedef struct Chunk {
        GLuint vertexArrayObject;
        GLuint vertices;
        int n_faces;
        int x, z;
    } Chunk;


    /** The different noise types supported. */
    enum NoiseType {
        PERLIN = 1,
        VALUE = 2,
        LAYERED = 3
    };


    /** Choose between CPU or GPU generation */
    enum GenMode {
        CPU = 1,
        GPU = 2
    };


    /** Class that manages generation of the terrain and which chunks
     *  to generate.
     */
    class Manager {
    public:
        const NoiseType style;
        const GenMode mode;
        const int chunkWidth, chunkHeight;
        const GLuint generate, march, renderShader;
        Camera::View *camera;

        Manager(NoiseType style, int width, int height, Camera::View *camera,
            GenMode mode, GLuint generate, GLuint march, GLuint renderShader);
        ~Manager();
        Chunk * get_chunk(int x, int z) { return chunks[z][x]; }
        std::vector<int> camera_position();
        void update();

    private:
        Chunk ***chunks;
        int last_position[2];

        // -- Methods --------------------------------------------------
        void delete_chunk(Chunk *chunk);

        /** Responsible for making the data that is passed to the
         * marching cubes methods for conversion. */
        Chunk* generate_chunk_CPU(int x, int z);
        Chunk** generate_chunks_GPU(int n_chunks, int coordinates[][2]);

        /** Responsible for actually making the chunk instance. */
        void march_cubes_CPU(Chunk* chunk, int **data);
        void march_cubes_GPU(Chunk **chunks, int n_chunks,
            GLuint chunk_positions, GLuint data);

        /** CPU Marching Cubes Functions */
        Rotations getMatch(int active, int x, int y, int z);
        int rotateActive(int active, const int *rotations);
        Rotations reverse(Rotations rotations);
        void buildMesh(Chunk *chunk, int totalPerms, Rotations* rotations);

        // -- CPU Noise Generation Functions ---------------------------
        // Noise Utility Functions.
        inline float lerp(float v1, float v2, float weight) {
            return v1 - ((v1 - v2) * weight);
        }
        float weight(float x);

        // Value noise functions.
        float value_noise(int x, int z, float seed);
        void gen_vertex_noise(float corners[4], int chunkX, int chunkZ, float seed);
        float gen_value_noise(int chunkX, int chunkZ, int x, int z, float seed);
        int ** value(int chunkX, int chunkZ, int layers);

        // Perlin noise functions.
        inline float dot(float vec1[2], float vec2[2]) {
            return (vec1[0] * vec2[0]) + (vec1[1] * vec2[1]);
        }
        float gradient_value(int x, int z, float seed, float freq);
        void gen_gradient_vectors(float gradients[4][2], int chunkX, int chunkZ, float seed, float freq);
        float calculate_perlin_noise(float gradients[4][2], int x, int z);
        int ** perlin(int chunkX, int chunkZ, int octaves, float freq_mult, float amp_mult);
    };


    // Constants.
    const float pi = 3.14159265359;
    extern const GLfloat *vertices;
    extern const Permutation *permutations;
    extern const int *permMap;
    const uint8_t mask = 3;

    // Rotation Maps.
    extern const int *xRotations;
    extern const int *yRotations;
    extern const int *zRotations;
    extern const int *reverseX;
    extern const int *reverseY;
    extern const int *reverseZ;
}
#endif