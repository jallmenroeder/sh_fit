//
// Created by jallmenroeder on 07/12/2020.
//

#include "ClampedCosine.h"


void ClampedCosine::setCoefficients(std::unique_ptr<std::vector<float>> coeffs) {
    if (coeffs->size() != 1) {
        printf("Error, clamped cosine can only handle 1 coefficient but got %d", coeffs->size());
        return;
    }
    m_amplitude = (*coeffs)[0];
}


float ClampedCosine::eval(const glm::vec3& V) const {
    return m_amplitude / M_PIf32 * fmaxf(V.z, 0.f);
}


float ClampedCosine::eval_basis(const glm::vec3& V, int idx) const {
    if (idx != 0) {
        printf("error, clamped cosine can only be evaluated for index 0");
        return 0.f;
    }
    return eval(V);
}

