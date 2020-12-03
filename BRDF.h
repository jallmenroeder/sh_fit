//
// Created by jallmenroeder on 24/11/2020.
//

#pragma once

#include <glm/glm.hpp>

#include "Util.h"
#include "SphericalFunction.h"

static const float F0 = .4f;

class BRDF : public SphericalFunction {
public:
    BRDF(const glm::vec3 view_dir, float roughness) : m_view_dir(view_dir), m_roughness(roughness) {}
    float eval(const glm::vec3& v) const override;
    glm::vec3 sample(const glm::vec2& uv);

private:
    static float eval_ggx(const glm::vec3& V, const glm::vec3& L, float roughness, float& pdf);

    glm::vec3 m_view_dir;
    float m_roughness;
};