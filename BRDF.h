//
// Created by jallmenroeder on 24/11/2020.
//

#pragma once

#include <glm/glm.hpp>

#include "Util.h"

static const float f0 = .4f;

static float eval_ggx(const glm::vec3& V, const glm::vec3& L, float roughness, float& pdf) {

    if (V.z <= 0.f || L.z <= 0.f) {
        pdf = 0.f;
        return 0.f;
    }

    // Halfvector
    glm::vec3 H = glm::normalize(V + L);

    // Normal Distribution D
    float roughness2 = roughness * roughness;
    float denom = (H.z * roughness2 - H.z) * H.z + 1.f;
    float D = roughness2 / denom / denom;

    // Fresnel
    float F = f0 + (1.f - f0) * powf(1.f - cdot(L, H), 5);

    // Visibility
    float G_V = L.z * sqrtf((-V.z * roughness2 + V.z) * V.z + roughness2);
    float G_L = V.z * sqrtf((-L.z * roughness2 + L.z) * L.z + roughness2);
    float vis = 0.5f / (G_V + G_L);
    pdf = fabsf(D * H.z / 4.0f / glm::dot(V,H));
    return D * vis * F / (float)M_PI;
}


static glm::vec3 sample(const glm::vec3& L, float roughness, const glm::vec2& uv) {
    const float phi = 2.0f * (float)M_PI * uv.x;
    const float r = roughness * sqrtf(uv.y / (1.0f - uv.y));
    const glm::vec3 N = normalize(glm::vec3(r*cosf(phi), r*sinf(phi), 1.0f));
    return -L + 2.0f * N * dot(N, L);
}
