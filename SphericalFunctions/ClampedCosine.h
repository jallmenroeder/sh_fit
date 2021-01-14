//
// Created by jallmenroeder on 02/12/2020.
//

#pragma once

#include <glm/glm.hpp>

#include "SphericalFunction.h"

class ClampedCosine : public SphericalFunction {
public:
    explicit ClampedCosine(float amplitude) : m_amplitude(amplitude) {}

    std::unique_ptr<SphericalFunction> copy() const override;

    void setCoefficients(std::shared_ptr<std::vector<float>> coeffs) override;
    std::shared_ptr<std::vector<float>> getCoefficients() override;

    int numCoefficients() const override { return 1; }

    float eval(const glm::vec3& V) const override;

    float eval_basis(const glm::vec3& V, int idx) const override;

private:
    float m_amplitude;

};
