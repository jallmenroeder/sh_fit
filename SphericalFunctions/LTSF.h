//
// Created by jallmenroeder on 07/12/2020.
//

#pragma once

#include <memory>

#include "BRDF.h"
#include "SphericalFunction.h"
#include "gsl/gsl_vector.h"


class LTSF : public BRDF {
public:
    LTSF(std::unique_ptr<SphericalFunction> spherical_function, const glm::mat3& M, const glm::vec3& view_dir,
         float roughness);

    void update(const glm::mat3& M);
    void setErrorResolution(int resolution);

    float eval(const glm::vec3& V) const override;
    float pdf(const glm::vec3& V) const override;
    float evalLtsfBasis(const glm::vec3& V, int idx) const;

    glm::vec3 sample(const glm::vec2& uv) const override;
    void findSphericalExpansion();
    double calculateError() const;

    /**
     * This function maps a linear transformation (defined by 4 parameters in vector x) to the difference between LTSF
     * and target function. It follows the GSL specifications for multidimensional functions for minimization
     * (https://www.gnu.org/software/gsl/doc/html/multimin.html#providing-a-function-to-minimize).
     * The minimum of this function is thus the closest approximation of the target function with LTSF.
     * @param x         4-component vector containing the linear transformation
     *                                  / a 0 b \
     *                  (a, b, c, d) -> | 0 c 0 |
     *                                  \ d 0 1 /
     * @param params    possible parameters, is needed by GSL specification but not used here
     * @return          difference between LTSF and target function
     */
    double minimizeFunc(const gsl_vector* x, void* params);

private:
    glm::vec3 linearly_transform_vec(const glm::vec3& V, float& jacobian) const;

    std::unique_ptr<SphericalFunction> m_spherical_function;
    std::unique_ptr<BRDF> m_target_function;
    glm::mat3 m_M;
    glm::mat3 m_M_inv;
    float m_det_M_inv;
    glm::vec3 m_view_dir;
    float m_roughness;
    int m_error_resolution = 64;

};
