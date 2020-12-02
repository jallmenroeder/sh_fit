//
// Created by jallmenroeder on 02/12/2020.
//

#pragma once

#include <glm/glm.hpp>

class SphericalFunction {
public:
    virtual float eval(const glm::vec3& v) const = 0;
};
