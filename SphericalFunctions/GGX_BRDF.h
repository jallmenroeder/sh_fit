//
// Created by jallmenroeder on 24/11/2020.
//

#pragma once

#include <glm/glm.hpp>

#include "../Util.h"
#include "BRDF.h"

static const float F0 = .4f;

class GGX_BRDF : public BRDF {
public:
    GGX_BRDF(const glm::vec3 view_dir, float roughness) : m_view_dir(view_dir), m_roughness(roughness) {}

    float eval(const glm::vec3& V) const override;
    float pdf(const glm::vec3& V) const override;
    glm::vec3 sample(const glm::vec2& uv) const override;

private:
    static float eval_ggx(const glm::vec3& V, const glm::vec3& L, float roughness);

    glm::vec3 m_view_dir;
    float m_roughness;
};