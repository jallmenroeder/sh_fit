//
// Created by jallmenroeder on 02/12/2020.
//

#pragma once

#include <glm/glm.hpp>
#include <vector>
#include <memory>

#include "../Util.h"

class SphericalFunction {
public:
    virtual ~SphericalFunction() = default;
    virtual uptr<SphericalFunction> copy() const = 0;
    virtual float eval(const glm::vec3& V) const = 0;
    virtual int numCoefficients() const = 0;
    virtual void setCoefficients(std::shared_ptr<std::vector<float>> coefficients) = 0;
    virtual std::shared_ptr<std::vector<float>> getCoefficients() = 0;
    virtual void evalBasisArray(const glm::vec3& V, std::vector<float>& array) const = 0;
    virtual std::string getName() const = 0;
};
