//
// Created by jallmenroeder on 07/12/2020.
//

#pragma once

#include <glm/glm.hpp>

/**
 * Interface class for BRDFs. Here, a BRDF has a fixed view direction and roughness and can be sampled and evaluated
 * for light directions.
 */
class BRDF {
public:
    virtual float eval(const glm::vec3& V) const = 0;
    virtual float pdf(const glm::vec3& V) const = 0;
    virtual glm::vec3 sample(const glm::vec2& uv) const = 0;
};