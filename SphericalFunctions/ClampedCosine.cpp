//
// Created by jallmenroeder on 07/12/2020.
//

#include "ClampedCosine.h"


void ClampedCosine::setCoefficients(std::shared_ptr<std::vector<float>> coeffs) {
    if (coeffs->size() != 1) {
        printf("Error, clamped cosine can only handle 1 coefficient but got %ld", coeffs->size());
        return;
    }
    m_amplitude = (*coeffs)[0];
}


std::unique_ptr<SphericalFunction> ClampedCosine::copy() const {
    return std::make_unique<ClampedCosine>(1.f);
}


std::shared_ptr<std::vector<float>> ClampedCosine::getCoefficients() {
    std::vector<float> coeffs = {m_amplitude};
    return std::make_shared<std::vector<float>>(coeffs);
}


float ClampedCosine::eval(const glm::vec3& V) const {
    return m_amplitude * eval_basis(V, 0);
}


float ClampedCosine::eval_basis(const glm::vec3& V, int idx) const {
    if (idx != 0) {
        printf("error, clamped cosine can only be evaluated for index 0");
        return 0.f;
    }
    return fmaxf(V.z, 0.f) / M_PIf32;
}

