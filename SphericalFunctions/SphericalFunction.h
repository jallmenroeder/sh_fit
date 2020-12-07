//
// Created by jallmenroeder on 02/12/2020.
//

#pragma once

#include <glm/glm.hpp>
#include <gsl/gsl_matrix.h>
#include <vector>

class SphericalFunction {
public:
    virtual float eval(const glm::vec3& V) const = 0;
    virtual gsl_matrix* create_lls_matrix(const std::vector<glm::vec3>& lin_trans_samples,
                                          const std::vector<float>& weights) const = 0;
};
