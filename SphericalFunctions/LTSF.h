//
// Created by jallmenroeder on 07/12/2020.
//

#pragma once

#include <memory>

#include "BRDF.h"
#include "SphericalFunction.h"


class LTSF : public BRDF {
public:
    LTSF(std::unique_ptr<SphericalFunction> spherical_function, const glm::mat3& M, const glm::vec3& view_dir,
         float roughness);
    float eval(const glm::vec3& V, float& pdf) const override;
    float eval_ltsf_basis(const glm::vec3& V, int idx) const;
    glm::vec3 sample(const glm::vec2& uv) const override;
    void findFit();

private:
    glm::vec3 linearly_transform_vec(const glm::vec3& V, float& jacobian) const;

    std::unique_ptr<SphericalFunction> m_spherical_function;
    std::unique_ptr<BRDF> m_target_function;
    glm::mat3 m_M;
    glm::mat3 m_M_inv;
    float m_det_M_inv;
    glm::vec3 m_view_dir;
    float m_roughness;

};
