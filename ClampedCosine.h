//
// Created by jallmenroeder on 02/12/2020.
//

#pragma once

#include <glm/glm.hpp>

#include "SphericalFunction.h"

class ClampedCosine : public SphericalFunction {
public:
    explicit ClampedCosine(float amplitude) : m_amplitude(amplitude) {}

    void setAmplitude(float amplitude) { m_amplitude = amplitude; }
    float getAmplitude() const { return m_amplitude; }

    float eval(const glm::vec3& v) const override {
        return m_amplitude * fmaxf(v.z, 0.f);
    }

private:
    float m_amplitude;

};
