//
// Created by jallmenroeder on 24/11/2020.
//

#pragma once

#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <glm/glm.hpp>
#include <vector>
#include <memory>

class SphericalHarmonics {
public:
    explicit SphericalHarmonics(int order);
    float eval(const glm::vec3& V);
    std::shared_ptr<std::vector<float>> getCoeffs() { return m_coeffs; }
    void setCoeffs(std::shared_ptr<std::vector<float>> coeffs);

private:
    const int m_ORDER;
    const int m_NUM_COEFFS;
    std::shared_ptr<std::vector<float>> m_coeffs;

    static float eval_SH(const glm::vec3& V, const std::vector<float>& coefficients, int order);
};


