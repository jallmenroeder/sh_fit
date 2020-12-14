//
// Created by jallmenroeder on 02/12/2020.
//

#pragma once

#include <glm/glm.hpp>
#include <gsl/gsl_matrix.h>
#include <vector>
#include <memory>

class SphericalFunction {
public:
    virtual float eval(const glm::vec3& V) const = 0;
    virtual int numCoefficients() const = 0;
    virtual void setCoefficients(std::unique_ptr<std::vector<float>> coefficients) = 0;
    virtual float eval_basis(const glm::vec3& V, int idx) const = 0;
};
