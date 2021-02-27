//
// Created by jallmenroeder on 24/11/2020.
//

#pragma once

#include <glm/glm.hpp>
#include <vector>
#include <memory>

#include "SphericalFunction.h"

class SphericalHarmonics : public SphericalFunction {
public:
    explicit SphericalHarmonics(int order);
    std::unique_ptr<SphericalFunction> copy() const override;
    float eval(const glm::vec3& V) const override;
    void evalBasisArray(const glm::vec3& V, std::vector<float>& array) const override;

    void setCoefficients(std::shared_ptr<std::vector<float>> coeffs) override;
    std::shared_ptr<std::vector<float>> getCoefficients() override;

    int numCoefficients() const override { return m_NUM_COEFFS; }
    std::string getName() const override { return "sh_n" + std::to_string(m_ORDER); }

private:
	static float K_l_m(int l, int m);

    const int m_ORDER;
    const int m_NUM_COEFFS;
    const int m_LEGENDRE_SIZE;
    std::shared_ptr<std::vector<float>> m_coeffs;
    std::vector<float> m_normalization;
};


