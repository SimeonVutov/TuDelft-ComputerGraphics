#pragma once
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
DISABLE_WARNINGS_POP()
#include <array>
#include <span>
#include <vector>

static glm::mat4 rotationMatrix(float angle, const glm::vec3& axis)
{
    return glm::rotate(glm::identity<glm::mat4>(), angle, axis);
}
static glm::mat4 translationMatrix(const glm::vec3& translation)
{
    return glm::translate(glm::identity<glm::mat4>(), translation);
}

// !!! DO NOT MODIFY !!!
struct TopSpinState {
    float armAngle { 0 };
    float seatsAngle { 0 };
    float harnessAngle { 0 };
};
// !!! DO NOT MODIFY !!!
struct TopSpinMatrices {
    glm::mat4 supports;
    glm::mat4 arms;
    glm::mat4 seats;
    std::array<glm::mat4, 44> harnesses;
};

// ==========================
// =====    EXERCISE    =====
// ==========================
inline TopSpinMatrices computeTopSpinTransformations(const TopSpinState& topSpinState)
{
    static constexpr int seatsPerRow = 22;
    const float armAngleInRadians = glm::radians(topSpinState.armAngle);
    const float seatsAngleInRadians = glm::radians(topSpinState.seatsAngle);
    const float harnessAngleInRadians = glm::radians(topSpinState.harnessAngle);

    TopSpinMatrices matrices;
    // Place the supports at the origin.
    matrices.supports = glm::identity<glm::mat4>(); // DO NOT CHANGE THIS LINE
    matrices.arms = glm::translate(matrices.supports, glm::vec3(0, 11, 0));
    matrices.arms = glm::rotate(matrices.arms, armAngleInRadians, glm::vec3(1, 0, 0));
    matrices.seats = glm::translate(matrices.arms, glm::vec3(0, -7.8f, 0));
    matrices.seats = glm::rotate(matrices.seats, seatsAngleInRadians, glm::vec3(1, 0, 0));
    matrices.seats = glm::translate(matrices.seats, glm::vec3(0, -3.0f, 0));    

    const float rowOffset[2] = {0.95f, -0.35f};
    const float seatWidth = 0.7f;
    const float startSeatX = -8.0f + (seatWidth / 2.0f);
    const float endSeatX = 8.0f - (seatWidth / 2.0f);
    const float gapBetweenSeats = (16.0f - (seatsPerRow * seatWidth)) / (seatsPerRow - 1);
    for(int row = 0; row < 2; row++) {
        for(int col = 0; col < seatsPerRow; col++) {
            glm::mat4 harness = matrices.seats;
            float x = startSeatX + col * (gapBetweenSeats + seatWidth);
            float y = 1.5f;
            float z = rowOffset[row];

            harness = glm::translate(harness, glm::vec3(x, y, z));
            harness = glm::rotate(harness, harnessAngleInRadians, glm::vec3(1, 0, 0));

            matrices.harnesses[(row * 22) + col] = harness;
        }
    }

    return matrices;
}
