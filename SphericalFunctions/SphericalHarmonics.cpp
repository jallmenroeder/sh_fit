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


gsl_matrix* SphericalHarmonics::create_lls_matrix(const std::vector<glm::vec3>& lin_trans_samples,
                                                  const std::vector<float>& weights) const {
    gsl_matrix* mat = gsl_matrix_alloc(lin_trans_samples.size(), m_NUM_COEFFS);

    std::vector<Spherical> spherical_samples(lin_trans_samples.size());
    for (int i = 0; i < lin_trans_samples.size(); i++) {
        spherical_samples[i] = cartesianToSpherical(lin_trans_samples[i]);
    }

    for (int row = 0; row < lin_trans_samples.size(); row++) {
        int column = 0;
        for (int l = 0; l <= m_ORDER; l++) {
            for (int m = -l; m <= l; m++) {
                double sh_basis_function = boost::math::spherical_harmonic_r(l, m, spherical_samples[row].theta,
                                                                             spherical_samples[row].phi);
                gsl_matrix_set(mat, row, column, sh_basis_function * weights[row]);
                column++;
            }
        }
    }
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