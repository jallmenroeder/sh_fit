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
    return std::make_unique<ClampedCosine>();
}


std::shared_ptr<std::vector<float>> ClampedCosine::getCoefficients() {
    std::vector<float> coeffs = {m_amplitude};
    return std::make_shared<std::vector<float>>(coeffs);
}


float ClampedCosine::eval(const glm::vec3& V) const {
    auto eval_array = std::vector<float>(1);
    evalBasisArray(V, eval_array);
    return m_amplitude * eval_array[0];
}


void ClampedCosine::evalBasisArray(const glm::vec3& V, std::vector<float> &array) const {
	assert(!array.empty());
	array[0] = fmaxf(V.z, 0.f) / M_PIf32;
}

