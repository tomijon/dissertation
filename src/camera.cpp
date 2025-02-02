#include <iostream>
#include <glm/gtc/type_ptr.hpp>

#include <camera.hpp>

using std::cout;

namespace Camera {
    /** Constructor for the View class.
     * 
     * Creates the projection matrix and updates the shader projection
     * uniform.
     * 
     * Range refers to the radius.
     */
    View::View(GLFWwindow *window, GLuint renderShader, int width, int height,
        int range) : window(window), renderShader(renderShader), range(range) {

        projectionMatrix = glm::perspective(
            glm::radians(45.0f), (float)width / (float)height, 0.1f, 2000.0f);

        glUniformMatrix4fv(glGetUniformLocation(renderShader, "projection"),
            1, GL_FALSE, glm::value_ptr(projectionMatrix));
    }


    /** Update the movement flags for active keys.
     * 
     *  Loops through every control key, activating the corresponding bit
     *  if the key is pressed, setting the bit to zero if released.
     */
    void View::updateMovementInputs() {
        for (int k = 0; k < n_controls; k++) {
            int state = glfwGetKey(window, keys[k]);

            // On Key Pressed.
            if (state == GLFW_PRESS) {
                movement_flags |= move_flag_values[k];
            } else if (state == GLFW_RELEASE) {
                movement_flags &= ~move_flag_values[k];
            }
        }
    }


    /** Move the position of the camera.
     *  
     *  Updates the position uniform in the render shader program.
     */
    void View::moveCamera(double dt) {
        if (movement_flags == move_flag_values[NONE]) return;

        int current_speed = speed;
        glm::vec3 resultant = glm::vec3(0);

        // Horrizontal plane movement.
        if (movement_flags & move_flag_values[FORWARD]) {
            resultant += forwardVector;
        }
        if (movement_flags & move_flag_values[LEFT]) {
            resultant += -rightVector;
        }
        if (movement_flags & move_flag_values[BACK]) {
            resultant += -forwardVector;
        }
        if (movement_flags & move_flag_values[RIGHT]) {
            resultant += rightVector;
        }
        resultant.y = 0; // Remove vertical movement.

        // Vertical plane movement.
        if (movement_flags & move_flag_values[UP]) {
            resultant += glm::vec3(0, 1, 0);
        }
        if (movement_flags & move_flag_values[DOWN]) {
            resultant += glm::vec3(0, -1, 0);
        }

        // Other movement related controls.
        if (movement_flags & move_flag_values[SPEED]) {
            current_speed *= 3;
        }

        // Prevent divide by zero.
        if (glm::length(resultant) == 0) return;

        resultant = glm::normalize(resultant);
        position.x += resultant.x * current_speed * dt;
        position.y += resultant.y * current_speed * dt;
        position.z += resultant.z * current_speed * dt;

        glUniform4fv(glGetUniformLocation(renderShader, "position"),
            1, glm::value_ptr(position));
    }


    /** Update the pitch and yaw angles of the camera.
     * 
     *  Pitch and Yaw change based on change in the x and y position of
     *  the mouse.
     */
    void View::updatePitchYaw(double dt) {
        double x, y;
        
        glfwGetCursorPos(window, &x, &y);
        glfwSetCursorPos(window, mouseX, mouseY);

        double dx = x - mouseX;
        double dy = y - mouseY;

        if (dx == 0 && dy == 0) { return; }

        // Pitch and Yaw calculations.
        yaw = std::max(std::min(yaw + (sensitivity* dy * dt), 90.0), -90.0);
        pitch += sensitivity * dx * dt;
        if (pitch < 0) pitch += 360;
        else if (pitch > 360) pitch -= (360 * (pitch / 360));
    }


    /** Updates the direction vectors forward, right and up. */
    void View::updateDirectionVectors() {
        float cosYaw = std::cos(glm::radians(yaw));
        float sinYaw = std::sin(glm::radians(yaw));
        float cosPitch = std::cos(glm::radians(pitch));
        float sinPitch = -std::sin(glm::radians(pitch));

        glm::mat4 yawRotation = glm::mat4(
            1, 0, 0, 0,
            0, cosYaw, sinYaw, 0,
            0, -sinYaw, cosYaw, 0,
            0, 0, 0, 1
        );

        glm::mat4 pitchRotation = glm::mat4(
            cosPitch, 0, -sinPitch, 0,
            0, 1, 0, 0,
            sinPitch, 0, cosPitch, 0,
            0, 0, 0, 1
        );

        forwardVector = glm::normalize(glm::vec3(
            pitchRotation * yawRotation * glm::vec4(0.0f, 0.0f, 1.0f, 0.0f)));

        rightVector = glm::normalize(glm::cross(
            forwardVector, glm::vec3(0.0f, 1.0f, 0.0f)));

        upVector = glm::normalize(glm::cross(rightVector, forwardVector));
    }


    /** Updates the view matrix.
     * 
     * Updates the view uniform in the render shader.
     */
    void View::updateViewMatrix() {
        viewMatrix = glm::lookAt(glm::vec3(position),
            glm::vec3(position) + forwardVector, upVector);

        glUniformMatrix4fv(glGetUniformLocation(renderShader, "view"),
            1, GL_FALSE, glm::value_ptr(viewMatrix));
    }
}