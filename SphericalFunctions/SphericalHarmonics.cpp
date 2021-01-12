//
// Created by jallmenroeder on 30/11/2020.
//

#include "SphericalHarmonics.h"

#include "../Util.h"

SphericalHarmonics::SphericalHarmonics(int order)
        : m_ORDER(order),
          m_NUM_COEFFS((order + 1) * (order + 1)),
          m_coeffs(new std::vector<float>(m_NUM_COEFFS, 0.f)) {}


float SphericalHarmonics::eval(const glm::vec3& V) const {
    float sum = 0.f;
    for (int coeff_idx = 0; coeff_idx < m_NUM_COEFFS; coeff_idx++) {
        sum += (*m_coeffs)[coeff_idx] * eval_basis(V, coeff_idx);
    }
    return sum;
}


float SphericalHarmonics::eval_basis(const glm::vec3& V, int idx) const {
    int l = floor(sqrt(idx) + FLT_EPSILON);
    int m = idx - l * (l + 1);
    Spherical V_ = cartesianToSpherical(V);
    if (m > 0) {
        return sqrtf(2) * boost::math::spherical_harmonic_r(l, m, V_.theta, V_.phi);
    } else if (m < 0) {
        return sqrtf(2) * boost::math::spherical_harmonic_i(l, m, V_.theta, V_.phi);
    } else {
        return boost::math::spherical_harmonic_r(l, m, V_.theta, V_.phi);
    }
}


void SphericalHarmonics::setCoefficients(std::shared_ptr<std::vector<float>> coeffs) {
    if (coeffs->size() != m_NUM_COEFFS) {
        printf("Error, SH order was %d, needing %d coefficients but got %ld coefficients, defaulting to 0 instead",
               m_ORDER, m_NUM_COEFFS, coeffs->size());
        return;
    }
    m_coeffs = std::move(coeffs);
}


std::shared_ptr<std::vector<float>> SphericalHarmonics::getCoefficients() {
    return m_coeffs;
}
