//
// Created by jallmenroeder on 02/12/2020.
//

#include "GGX_BRDF.h"


float GGX_BRDF::eval_ggx(const glm::vec3& V, const glm::vec3& L, float roughness, float& pdf) {

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
    float F = F0 + (1.f - F0) * powf(1.f - cdot(L, H), 5);

    // Visibility
    float G_V = L.z * sqrtf((-V.z * roughness2 + V.z) * V.z + roughness2);
    float G_L = V.z * sqrtf((-L.z * roughness2 + L.z) * L.z + roughness2);
    float vis = 0.5f / (G_V + G_L);
    pdf = fabsf(D * H.z / 4.0f / glm::dot(V, H));
    return D * vis * F / M_PIf32;
}


glm::vec3 GGX_BRDF::sample(const glm::vec2& uv) const {
    const float phi = 2.0f * M_PIf32 * uv.x;
    const float r = m_roughness * sqrtf(uv.y / (1.0f - uv.y));
    const glm::vec3 N = normalize(glm::vec3(r * cosf(phi), r * sinf(phi), 1.0f));
    return glm::normalize(-m_view_dir + 2.0f * N * dot(N, m_view_dir));
}


float GGX_BRDF::eval(const glm::vec3& v, float& pdf) const {
    return eval_ggx(m_view_dir, v, m_roughness, pdf);
}