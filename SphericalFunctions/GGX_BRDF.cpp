//
// Created by jallmenroeder on 02/12/2020.
//

#include "GGX_BRDF.h"


float GGX_BRDF::eval_ggx(const glm::vec3& V, const glm::vec3& L) const {

    if (V.z <= 0.f || L.z <= 0.f) {
        return 0.f;
    }

    // Halfvector
    glm::vec3 H = glm::normalize(V + L);

    // Normal Distribution D
    float denom = (H.z * m_roughness2 - H.z) * H.z + 1.f;
    float D = m_roughness2 / denom / denom;

    // Fresnel
    float F = F0 + (1.f - F0) * powf(1.f - cdot(L, H), 5);

    // Visibility
    float G_V = L.z * sqrtf((-V.z * m_roughness2 + V.z) * V.z + m_roughness2);
    float G_L = V.z * sqrtf((-L.z * m_roughness2 + L.z) * L.z + m_roughness2);
    float vis = 0.5f / (G_V + G_L);
    return D * vis * F / M_PIf32 * V.z;
}


float GGX_BRDF::pdf(const glm::vec3& V) const {
    if (V.z <= 0) return 0.f;

    const auto H = glm::normalize(V + m_view_dir);
    float theta = acosf(H.z);
    float roughness_2 = m_roughness * m_roughness;
    float denom = (roughness_2 - 1.f) * H.z * H.z + 1.f;
    denom *= denom * M_PIf32;
    return roughness_2 * H.z * sinf(theta) / denom;
}


glm::vec3 GGX_BRDF::sample(const glm::vec2& uv) const {
    const float theta = acosf(sqrtf((1 - uv.x) / (uv.x * (m_roughness * m_roughness - 1.f) + 1.f)));
    const float phi = 2.0f * M_PIf32 * uv.y;
    const glm::vec3 N = sphericalToCartesian({theta, phi});
    return glm::normalize(-m_view_dir + 2.0f * N * dot(N, m_view_dir));
}


float GGX_BRDF::eval(const glm::vec3& v) const {
    return eval_ggx(m_view_dir, v);
}