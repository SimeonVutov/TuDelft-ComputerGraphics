#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <span>
#include <tuple>
#include <utility>
#include <vector>

// ==================================
// ========    Exercise 1    ========
// ==================================
std::pair<float, float> statistics(std::span<const float> values)
{
    if (values.empty()) {
        return { 0.0f, 0.0f };
    }

    float mean = 0.0f;
    for(float v : values) {
        mean += v;
    }

    mean /= static_cast<float>(values.size());
    
    float std = 0.0f;
    for(float v : values) {
        std += (v - mean) * (v - mean);
    }
    std /= static_cast<float>(values.size());
    std = std::sqrt(std);

    return std::pair<float, float> {mean, std};
}

// ==================================
// ========    Exercise 2    ========
// ==================================
struct Tree {
    float value;
    std::vector<Tree> children;
};

float countRecursive(const Tree &node, int level, bool countOnlyEvenLevels) {
    float sum = 0.0f;

    if(!countOnlyEvenLevels || (level % 2 ==0)) {
        sum += node.value;
    }

    for(const auto &child : node.children) {
        sum += countRecursive(child, level + 1, countOnlyEvenLevels);
    }

    return sum;
}
// Your implementation goes here. Feel free to define helper functions/structs/classes in this
// file if you need to. Make sure to put them above this function to prevent linker warnings.
float countTree(const Tree& tree, bool countOnlyEvenLevels)
{
    return countRecursive(tree, 0, countOnlyEvenLevels);
}
