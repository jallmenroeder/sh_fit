//
// Created by jallmenroeder on 02/12/2020.
//

#pragma once

#include <glm/glm.hpp>

#include "SphericalFunction.h"

class ClampedCosine : public SphericalFunction {
public:
    explicit ClampedCosine() : m_amplitude(1.f) {}

    std::unique_ptr<SphericalFunction> copy() const override;

    void setCoefficients(std::shared_ptr<std::vector<float>> coeffs) override;
    std::shared_ptr<std::vector<float>> getCoefficients() override;

    int numCoefficients() const override { return 1; }
	std::string getName() const override {return "cos"; }

    float eval(const glm::vec3& V) const override;

    void evalBasisArray(const glm::vec3& V, std::vector<float>& array) const override;

private:
    float m_amplitude;

};
