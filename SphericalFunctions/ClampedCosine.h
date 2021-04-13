//
// Created by jallmenroeder on 02/12/2020.
//

#pragma once

#include <glm/glm.hpp>

#include "SphericalFunction.h"

class ClampedCosine : public SphericalFunction {
public:
    explicit ClampedCosine() : m_amplitude(1.f) {}


    void setCoefficients(std::shared_ptr<std::vector<float>> coeffs) override {
        if (coeffs->size() != 1) {
            printf("Error, clamped cosine can only handle 1 coefficient but got %ld", coeffs->size());
            return;
        }
        m_amplitude = (*coeffs)[0];
    }


    uptr<SphericalFunction> copy() const override {
        return std::make_unique<ClampedCosine>();
    }


    std::shared_ptr<std::vector<float>> getCoefficients() override {
        std::vector<float> coeffs = {m_amplitude};
        return std::make_shared<std::vector<float>>(coeffs);
    }


    int numCoefficients() const override { return 1; }
	std::string getName() const override {return "cos"; }


    float eval(const glm::vec3& V) const override {
        auto eval_array = std::vector<float>(1);
        evalBasisArray(V, eval_array);
        return m_amplitude * eval_array[0];
    }


    void evalBasisArray(const glm::vec3& V, std::vector<float> &array) const override {
        assert(!array.empty());
        array[0] = fmaxf(V.z, 0.f) / M_PIf32;
    }

private:
    float m_amplitude;
};
