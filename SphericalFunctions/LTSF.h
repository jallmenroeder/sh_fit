//
// Created by jallmenroeder on 07/12/2020.
//

#pragma once

#include <memory>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multimin.h>

#include "BRDF.h"
#include "SphericalFunction.h"

struct Idx {
    int view;
    int roughness;
};


class LTSF : public BRDF {
public:
    LTSF(std::unique_ptr<SphericalFunction> spherical_function, const glm::mat3& M, Idx idx, int LUT_dimension);
    ~LTSF();

    void update(const glm::mat3& M);

    Idx getIdx() const { return m_idx; };
    std::shared_ptr<glm::mat3> getLinearTransformation() const { return m_M; };
    std::shared_ptr<glm::mat3> getInvLinearTransformation() const { return m_M_inv; };
    std::shared_ptr<std::vector<float>> getCoefficients() const { return m_spherical_function->getCoefficients(); }
    std::unique_ptr<SphericalFunction> getSphericalFunctionCopy() const;
    float getResidual() const { return m_residual; }

    float eval(const glm::vec3& V) const override;
    float pdf(const glm::vec3& V) const override;
    void evalLtsfBasisArray(const glm::vec3& V, std::vector<float>& array, float& jacobian) const;

    glm::vec3 sample(const glm::vec2& uv) const override;
    double findSphericalExpansion();

    /**
     * This function maps a linear transformation (defined by 4 parameters in vector m_multimin_x) to the difference between LTSF
     * and target function. It follows the GSL specifications for multidimensional functions for minimization
     * (https://www.gnu.org/software/gsl/doc/html/multimin.html#providing-a-function-to-minimize).
     * The minimum of this function is thus the closest approximation of the target function with LTSF.
     * @param x         4-component vector containing the linear transformation
     *                                  / a 0 b \
     *                  (a, b, c, d) -> | 0 c 0 |
     *                                  \ d 0 1 /
     * @param params    pointer to the LTSF instance since this function must be static
     * @return          difference between LTSF and target function
     */
    static double minimizeFunc(const gsl_vector* x, void* params);

    void findFit();
private:
    glm::vec3 linearly_transform_vec(const glm::vec3& V, float& jacobian) const;

    std::unique_ptr<SphericalFunction> m_spherical_function;
    std::unique_ptr<BRDF> m_target_function;
    std::shared_ptr<glm::mat3> m_M;
    std::shared_ptr<glm::mat3> m_M_inv;
    float m_det_M_inv;
    Idx m_idx;
    glm::vec3 m_view_dir;
    float m_roughness;
    std::vector<float> m_basis_eval;
    float m_residual;

    // gsl multifit
	gsl_matrix* m_gsl_mat;
	gsl_matrix* m_gsl_cov;
	gsl_vector* m_gsl_target_vector;
	gsl_vector* m_gsl_coefficients;
	gsl_vector* m_gsl_weights;
	gsl_multifit_linear_workspace* m_workspace;

	// gsl multimin
	const gsl_multimin_fminimizer_type *m_multimin_type;
	gsl_multimin_fminimizer *m_multimin_workspace;
	gsl_vector *m_multimin_step_size, *m_multimin_x;
};
