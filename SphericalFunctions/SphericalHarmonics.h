//
// Created by jallmenroeder on 24/11/2020.
//

#pragma once

#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <glm/glm.hpp>
#include <vector>
#include <memory>

#include "SphericalFunction.h"

class SphericalHarmonics : public SphericalFunction {
public:
    explicit SphericalHarmonics(int order);
    std::unique_ptr<SphericalFunction> copy() const override;
    float eval(const glm::vec3& V) const override;
    float eval_basis(const glm::vec3& V, int idx) const override;

    void setCoefficients(std::shared_ptr<std::vector<float>> coeffs) override;
    std::shared_ptr<std::vector<float>> getCoefficients() override;

    int numCoefficients() const override { return m_NUM_COEFFS; }

private:
    const int m_ORDER;
    const int m_NUM_COEFFS;
    std::shared_ptr<std::vector<float>> m_coeffs;
};


