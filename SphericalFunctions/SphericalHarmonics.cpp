//
// Created by jallmenroeder on 30/11/2020.
//

#include "SphericalHarmonics.h"

#include "../Util.h"

SphericalHarmonics::SphericalHarmonics(int order)
: m_ORDER(order),
  m_NUM_COEFFS((order + 1) * (order + 1)),
  m_coeffs(std::make_shared<std::vector<float>>(m_NUM_COEFFS, 0.f)) {}


float SphericalHarmonics::eval(const glm::vec3& V) const {
    return eval_SH(V, *m_coeffs, m_ORDER);
}


void SphericalHarmonics::setCoeffs(std::shared_ptr<std::vector<float>> coeffs) {
    if (coeffs->size() != m_NUM_COEFFS) {
        printf("Error, SH order was %d, needing %d coefficients but got %d coefficients, defaulting to 0 instead",
               m_ORDER, m_NUM_COEFFS, coeffs->size());
        return;
    }
    m_coeffs = coeffs;

}


float SphericalHarmonics::eval_SH(const glm::vec3& V, const std::vector<float>& coefficients, int order) {
    Spherical V_ = cartesianToSpherical(V);
    float sum = 0.f;
    int idx = 0;
    for (int l = 0; l <= order; l++) {
        for (int m = -l; m <= l; m++) {
            sum += coefficients[idx] * boost::math::spherical_harmonic_r(l, m, V_.theta, V_.phi);
            idx++;
        }
    }
    return sum;
}