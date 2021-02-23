//
// Created by jallmenroeder on 02/12/2020.
//

#include "GGX_BRDF.h"


GGX_BRDF::GGX_BRDF(const glm::vec3& view_dir, float roughness)
	: m_V(view_dir),
	  m_roughness(roughness),
	  m_roughness2(roughness * roughness) {

	// precomputation for sampling, never change for one BRDF configuration
	// Sampling the GGX Distribution of visible normals, see http://jcgt.org/published/0007/04/01/ for more detail
	m_Vh = normalize(glm::vec3(m_roughness * m_V.x, m_roughness * m_V.y, m_V.z));
	float lensq = m_Vh.x * m_Vh.x + m_Vh.y * m_Vh.y;
	m_T1 = lensq > 0 ? glm::vec3(-m_Vh.y, m_Vh.x, 0) * (1.f / sqrtf(lensq)) : glm::vec3(1, 0, 0);
	m_T2 = cross(m_Vh, m_T1);
	m_G1_V = 1.f / (1.f + ((-1.f + sqrtf(1.f + (m_V.x * m_V.x * m_roughness2 + m_V.y * m_V.y * m_roughness2) / (m_V.z * m_V.z))) * .5f));
}

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
    // add L.z because we evaluate the cosine weighted BRDF
    return D * vis * F / M_PIf32 * L.z;
}


float GGX_BRDF::pdf(const glm::vec3 &L) const {
	// Halfvector
	glm::vec3 H = glm::normalize(L + m_V);

	// Normal Distribution D
	float temp = (H.x * H.x / m_roughness2 + H.y * H.y / m_roughness2 + H.z * H.z);
	float D = 1.f / (M_PIf32 * m_roughness2 * temp * temp);

	return (m_G1_V * fmaxf(0, glm::dot(m_V, H)) * D) / (4.f * glm::dot(m_V, H) * m_V.z);
}


// Sampling the GGX Distribution of visible normals
// http://jcgt.org/published/0007/04/01/
glm::vec3 GGX_BRDF::sample(const glm::vec2& uv) const {
	float r = sqrtf(uv.x);
	float phi = 2.f * M_PIf32 * uv.y;
	float t1 = r * cosf(phi);
	float t2 = r * sinf(phi);
	float s = 0.5f * (1.f + m_Vh.z);
	t2 = (1.f - s)*sqrtf(1.f - t1*t1) + s*t2;
	glm::vec3 Nh = t1*m_T1 + t2*m_T2 + sqrtf(fmaxf(0.f, 1.f - t1*t1 - t2*t2))*m_Vh;
	glm::vec3 Ne = normalize(glm::vec3(m_roughness * Nh.x, m_roughness * Nh.y, fmaxf(0.0, Nh.z)));
	return glm::reflect(-m_V, Ne);
}


float GGX_BRDF::eval(const glm::vec3& L) const {
    return eval_ggx(m_V, L);
}