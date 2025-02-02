#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <camera.hpp>
#include <generation.hpp>

using std::cout;
using std::endl;
using std::string;


// Some (one) globals.
float view_angle = 140;


/** Reads the shader source file.
 * 
 *  Parameters:
 *  string fileName - filename of the shader source.
 * 
 *  Returns:
 *  string containing the whole source file.
 */
string readFile(string fileName) {
    std::ifstream file (fileName);
    string source = "";

    if (file.is_open()) {
        std::stringstream sourceBuffer;
        sourceBuffer << file.rdbuf();
        source = sourceBuffer.str();
        file.close();
        return source;
    } else throw std::runtime_error("Shader file not found: " + fileName);
}


/** Compiles a shader and returns the GLuint it points to.
 * 
 *  Parameters:
 *  const char* source - One long string containing all the shader
 *                       source code.
 * 
 *  Returns:
 *  GL_FALSE if failed. GL_TRUE if success.
 */
 GLboolean compileShader(GLuint shader, const char* source) {
    if (!source) throw std::invalid_argument("Source is null.");
    glShaderSource(shader, 1, &source, NULL);
    glCompileShader(shader);
    GLboolean success;
    glGetShaderiv(shader, GL_COMPILE_STATUS, (int*)&success);
    return success;
}


/** Sets which version of OpenGL core should be used for glfw.
 * 
 *  Example:
 *  major = 4, minor = 6 results in version 4.6
 *  
 *  Parameters:
 *  int major - the major version number
 *  int minor - the minor version number 
 */
void setGLVersion(int major, int minor) {
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, major);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, minor);
}


/** Load glad and enable some gl options. */
void initOpenGL() {
    gladLoadGL();
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); 
}


/** Creates a window using glfw.
 *  
 *  Sets the GL version to 4.6
 *  Enables opengl and sets window as the current context.
 *  Makes window the input mode.
 *  
 *  Parameters:
 *  int width - the width of the window in pixels.
 *  int height - the height of the window in pixels.
 *  const char* title - the window title text. 
 * 
 *  Returns:
 *  A GLFWwindow pointer to the window.
 */
GLFWwindow* createWindow(int width, int height, const char* title) {
    if (!glfwInit()) throw std::runtime_error("Failed to initialize GLFW.");
    setGLVersion(4, 6);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_DEPTH_BITS, 24);
    
    GLFWwindow* window = glfwCreateWindow(width, height, title, NULL, NULL);
    if (!window) throw std::runtime_error("Failed to create window.");
    
    glfwMakeContextCurrent(window);
    initOpenGL();
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_HIDDEN);
    return window;
}


/** Create a shader.
 * 
 *  This function technically allocates memory. Remember to delete your
 *  shader when finished with.
 * 
 *  Parameters:
 *  string fileName - filename of the shader to compile.
 * 
 *  Returns:
 *  GLuint of the shader.
 */
GLuint createShader(string fileName, GLenum type) {
    GLuint shader = glCreateShader(type);
    string source = readFile(fileName);
    GLboolean result = compileShader(shader, source.c_str());

    // Output compilation errors.
    if (result == GL_FALSE) {
        int length, BUF_SIZE = 1024;
        char buffer[BUF_SIZE];

        glGetShaderInfoLog(shader, BUF_SIZE, &length, buffer);
        cout << "Error code: " << glGetError() << endl;
        cout << buffer << endl;
        throw std::runtime_error("Shader Compilation Failed.");
    }
    return shader;
}


/** Intermediate linking function that will display error output if
 *  there is any.
 */
void linkProgram(GLuint program) {
    glLinkProgram(program);

    // Check for linking errors.
    int success, length, BUF_SIZE = 1024;
    char buffer[BUF_SIZE];
    glGetProgramiv(program, GL_LINK_STATUS, &success);

    if (!success) {
        glGetProgramInfoLog(program, BUF_SIZE, &length, buffer);
        cout << buffer << endl;
        throw std::runtime_error("Shader Linking Failed.");
    }
}


/** Links a vertex shader and a fragment shader into a single program.
 * 
 *  Parameters:
 *  GLuint vertex - the vertex shader;
 *  GLuint fragment - the fragment shader;
 * 
 *  Returns:
 *  GLuint representing the shader.
 */
GLuint linkShaderProgram(GLuint vertex, GLuint fragment) {
    GLuint shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertex);
    glAttachShader(shaderProgram, fragment);
    linkProgram(shaderProgram);

    glDeleteShader(vertex);
    glDeleteShader(fragment);
    return shaderProgram;
}


/** Create a shader program for a compute shader.
 */
GLuint computeProgram(string name) {
    GLuint compute = createShader(name, GL_COMPUTE_SHADER);
    GLuint program = glCreateProgram();
    glAttachShader(program, compute);
    linkProgram(program);
    glDeleteShader(compute);
    return program;
}


/** Initilise the vertex and fragment shaders.
 *  
 *  Looks for shaders called vertex.glsl and fragment.glsl
 *  
 *  Parameters:
 *  string location - the location of the vertex and fragment shaders.
 * 
 */
GLuint initRenderShader(string location) {
    GLuint vertexShader = createShader(location + "/vertex.glsl", GL_VERTEX_SHADER);
    GLuint fragmentShader = createShader(location + "/fragment.glsl", GL_FRAGMENT_SHADER);
    return linkShaderProgram(vertexShader, fragmentShader);
}


/** Scroll wheel callback function.
 * 
 *  Returns the change in y.
 */
void scroll_wheel_callback(GLFWwindow *window, double dx, double dy) {
    view_angle = std::max(0.0, std::min(180.0, view_angle + (5 * dy)));
}

 
/** Main method. Pretty self explanatory. */
int main(int n_args, char *args[]) {
    // Window and OpenGL initialisation.
    GLFWwindow* window = createWindow(1920, 1080, "Procedural Terrain");
    GLuint renderShader = initRenderShader("shaders");

    glfwSetScrollCallback(window, scroll_wheel_callback);
    
    // Compute shader creation.
    GLuint generate = computeProgram("shaders/generate.glsl");
    cout << "generate compiled and linked\n";
    GLuint march = computeProgram("shaders/march.glsl");
    cout << "march compiled and linked\n";

    glUseProgram(generate);
    glDispatchCompute(32, 1, 32);
    glUseProgram(march);
    glDispatchCompute(31, 1, 31);
    glUseProgram(renderShader);

    // Colors and uniforms and stuff.
    glm::vec4 triangleColor = glm::vec4(0.4, 0.66, 0.26, 1);
    glUniform4fv(glGetUniformLocation(renderShader, "drawColor"),
            1, glm::value_ptr(triangleColor));
    glClearColor(0, 0, 0, 1);

    glm::vec3 lightDirection = glm::normalize(glm::vec3(-60, -40, -60));
    glUniform3fv(glGetUniformLocation(renderShader, "lightDirection"),
        1, glm::value_ptr(lightDirection));

    glm::mat4 model = glm::mat4(1.0f);
    glUniformMatrix4fv(glGetUniformLocation(renderShader, "model"),
        1, GL_FALSE, glm::value_ptr(model));

    // Timing variables.
    auto last = std::chrono::high_resolution_clock::now();
    auto now = std::chrono::high_resolution_clock::now();
    std::chrono::duration<float> duration = last - now;
    float dt = duration.count();

    int range = 10;
    int width = 16;
    int height = 64;
    Generation::GenMode mode = Generation::CPU;
    Generation::NoiseType type = Generation::PERLIN;

    // Terrain generation specific variables.
    Camera::View *camera = new Camera::View(window, renderShader, 1920, 1080, range);
    Generation::Manager *manager = new Generation::Manager(
        type, width, height, camera, mode, generate, march, renderShader);

    // Program loop.
    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();

        // Leave the program.
        if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
            glfwSetWindowShouldClose(window, GLFW_TRUE);
            continue;
        }

        // Timings.
        now = std::chrono::high_resolution_clock::now();
        duration = now - last;
        dt = duration.count();
        last = now;

        // cout << std::fixed << std::setprecision(1) << 1 / dt << endl;

        // Update stuff.
        camera->update(dt);
        manager->update();

        // Begin rendering.
        glUseProgram(renderShader);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glUniform4fv(glGetUniformLocation(renderShader, "drawColor"),
            1, glm::value_ptr(triangleColor));



        float max_angle = view_angle * 3.14 / 180.0;
        float max_value = std::cos(max_angle);


        int toggle = 0;

        // Triangle drawing.
        for (int x = 0; x < camera->view_diameter(); x++) {
            for (int z = 0; z < camera->view_diameter(); z++) {
                Generation::Chunk *chunk = manager->get_chunk(x, z);

                // Determine angle between camera look vector and the vector
                // between the camera position and chunk position.
                std::vector<int> cam_pos = manager->camera_position();
                glm::vec3 direction = glm::vec3(0);
                direction.x = cam_pos[0] - chunk->x;
                direction.z = cam_pos[1] - chunk->z;
                direction = glm::normalize(direction);

                glm::vec3 look_vec = glm::vec3(0);
                look_vec.x = camera->forward().x;
                look_vec.z = camera->forward().z;
                look_vec = glm::normalize(look_vec);

                float cosangle = glm::dot(direction, look_vec);

                if (cosangle < max_value) {
                    glBindVertexArray(chunk->vertexArrayObject);
                    glDrawArrays(GL_TRIANGLES, 0, chunk->n_faces);
                    glBindBuffer(GL_ARRAY_BUFFER, 0);
                    continue;
                }
            }
        }
        glfwSwapBuffers(window);
        glBindVertexArray(0);
    }
    return 0;
}