#ifndef CAMERA_HPP
#define CAMERA_HPP

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>


namespace Camera {

    // Controls used for camera operation.
    enum controls_t {
        NONE = 0,
        FORWARD = 1,
        BACK = 2,
        LEFT = 3,
        RIGHT = 4,
        UP = 5,
        DOWN = 6,
        SPEED = 7
    };
    const int n_controls = 8;


    // Keys related to the generalised control names.
    const int keys[n_controls] = {
        GLFW_KEY_UNKNOWN,       // NONE
        GLFW_KEY_W,             // FORWARD
        GLFW_KEY_S,             // BACK
        GLFW_KEY_A,             // LEFT
        GLFW_KEY_D,             // RIGHT
        GLFW_KEY_SPACE,         // UP
        GLFW_KEY_LEFT_SHIFT,    // DOWN
        GLFW_KEY_LEFT_CONTROL   // SPEED
    };


    // Used to determine which bit the movement applies to.
    const unsigned int move_flag_values[n_controls + 1] {
        0,              // NONE
        1,              // FORWARD
        1 << 1,         // BACK
        1 << 2,         // LEFT
        1 << 3,         // RIGHT
        1 << 4,         // UP
        1 << 5,         // DOWN
        1 << 6          // SPEED
    };


    /** View class for handling camera matrices and movement. */
    class View {
    public:
        double sensitivity = 20.0;
        double speed = 10.0;

        View(GLFWwindow *window, GLuint renderShader, int width, int height,
            int range);

        int view_range() { return range; }
        int view_diameter() { return (range * 2) + 1; }
        float get_pitch() { return pitch; }
        float get_yaw() { return yaw; }
        glm::vec4 get_position() { return position; }
        glm::vec3 forward() { return forwardVector; }
        glm::vec3 right() { return rightVector; }
        glm::vec3 up() {return upVector; }
        glm::mat4 view_matrix() { return viewMatrix; }

        void update(double dt) {
            updateMovementInputs();
            updatePitchYaw(dt);
            updateDirectionVectors();
            moveCamera(dt);
            updateViewMatrix();
        }

    private:
        GLFWwindow *window;
        GLuint renderShader;

        glm::vec4 position = glm::vec4(0);
        glm::vec3 forwardVector;
        glm::vec3 rightVector;
        glm::vec3 upVector;

        glm::mat4 viewMatrix;
        glm::mat4 projectionMatrix;

        double pitch = 90;
        double yaw = 0;
        double mouseX = 400;
        double mouseY = 400;
        unsigned int movement_flags = 0;
        int range;

        void updateMovementInputs();
        void moveCamera(double dt);
        void updatePitchYaw(double dt);
        void updateDirectionVectors();
        void updateViewMatrix();
    };
}
#endif