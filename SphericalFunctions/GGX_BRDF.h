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
    GGX_BRDF(const glm::vec3& view_dir, float roughness);

    /**
     * Returns the !cosine weighted! GGX BRDF
     */
    float eval(const glm::vec3& L) const override;
    float pdf(const glm::vec3& L) const override;
    glm::vec3 sample(const glm::vec2& uv) const override;

private:
    float eval_ggx(const glm::vec3& V, const glm::vec3& L) const;

    glm::vec3 m_V;
    float m_roughness;
    float m_roughness2;

    // precomputed values for sampling
    glm::vec3 m_Vh;
	glm::vec3 m_T1;
	glm::vec3 m_T2;
};