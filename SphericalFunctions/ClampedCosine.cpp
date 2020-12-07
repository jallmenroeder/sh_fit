//
// Created by jallmenroeder on 07/12/2020.
//

#include "ClampedCosine.h"


float ClampedCosine::eval(const glm::vec3& V) const {
    return m_amplitude / M_PIf32 * fmaxf(V.z, 0.f);
}


gsl_matrix* ClampedCosine::create_lls_matrix(const std::vector<glm::vec3>& samples,
                                             const std::vector<float>& weights) const {
    gsl_matrix* mat = gsl_matrix_alloc(samples.size(), 1);
    for (int row = 0; row < samples.size(); row++) {
        gsl_matrix_set(mat, row, 1, eval(samples[row]) * weights[row]);
    }
    return mat;
}

