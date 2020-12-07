//
// Created by jallmenroeder on 02/12/2020.
//

#pragma once

#include <glm/glm.hpp>

#include "SphericalFunction.h"

class ClampedCosine : public SphericalFunction {
public:
    explicit ClampedCosine(float amplitude) : m_amplitude(amplitude) {}

    void setAmplitude(float amplitude) { m_amplitude = amplitude; }

    float getAmplitude() const { return m_amplitude; }

    float eval(const glm::vec3& V) const override;
    gsl_matrix* create_lls_matrix(const std::vector<glm::vec3>& samples,
                                  const std::vector<float>& weights) const override;

private:
    float m_amplitude;

};
