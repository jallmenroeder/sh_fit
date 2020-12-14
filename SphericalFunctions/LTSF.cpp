//
// Created by jallmenroeder on 07/12/2020.
//

#include "LTSF.h"

#include <gsl/gsl_multifit.h>
#include <iostream>

#include "../Util.h"
#include "GGX_BRDF.h"


LTSF::LTSF(std::unique_ptr<SphericalFunction> spherical_function, const glm::mat3& M, const glm::vec3& view_dir,
           float roughness) :
        m_spherical_function(std::move(spherical_function)),
        m_target_function(new GGX_BRDF(view_dir, roughness)),
        m_M(M),
        m_M_inv(glm::inverse(M)),
        m_det_M_inv(glm::determinant(m_M_inv)),
        m_view_dir(view_dir),
        m_roughness(roughness) {}


float LTSF::eval(const glm::vec3& V, float& pdf) const {
    // constant pdf due to uniform sampling
    pdf = 1.f / 2.f * M_PIf32;

    // evaluate original spherical function
    float jacobian;
    glm::vec3 trans_V = linearly_transform_vec(V, jacobian);
    return m_spherical_function->eval(trans_V) * jacobian;
}


float LTSF::eval_ltsf_basis(const glm::vec3& V, int idx) const {
    float jacobian;
    glm::vec3 trans_V = linearly_transform_vec(V, jacobian);
    return m_spherical_function->eval_basis(trans_V, idx) * jacobian;
}

// uniform sampling on the hemisphere
glm::vec3 LTSF::sample(const glm::vec2& uv) const {
    const float pi_2_v = 2.f * M_PIf32 * uv.y;
    const float sqrt_1_u = sqrtf(1.f - uv.x * uv.x);
    return glm::normalize(glm::vec3(cosf(pi_2_v) * sqrt_1_u, sinf(pi_2_v) * sqrt_1_u, uv.y));
}


void LTSF::findFit() {
//    int rows = NUM_SAMPLES * NUM_SAMPLES * 2;
    int rows = NUM_SAMPLES * NUM_SAMPLES * 4;
    int columns = m_spherical_function->numCoefficients();
    auto gsl_mat = gsl_matrix_alloc(rows, columns);
    auto gsl_cov = gsl_matrix_alloc(columns, columns);
    auto gsl_target_vector = gsl_vector_alloc(rows);
    auto gsl_coefficients = gsl_vector_alloc(columns);

    float pdf_ltsf = 1.f / 2.f * M_PIf32;
//    int row_idx = 0;

    std::vector<float> jacobian(rows);
//    for (int i = 0; i < NUM_SAMPLES; i++) {
//        for (int j = 0; j < NUM_SAMPLES; j++) {
    for (int i = 0; i < NUM_SAMPLES; i++) {
        for (int j = 0; j < 4 * NUM_SAMPLES; j++) {
            int row_idx = i * 4 * NUM_SAMPLES + j;
            const float theta = ((float)i + 0.5f) / (float)NUM_SAMPLES * M_PIf32 / 2.f;
            const float phi = ((float)j + 0.5f) / (float)NUM_SAMPLES * M_PIf32 * 2.f;

            auto V = sphericalToCartesian({theta, phi});
            float _;
            float eval_target = m_target_function->eval(V, _);
            float weight = sinf(theta);
            for (int column_idx = 0; column_idx < columns; column_idx++) {
                float eval_ltsf = eval_ltsf_basis(V, column_idx);
                gsl_matrix_set(gsl_mat, row_idx, column_idx, eval_ltsf * weight);
            }
            // we seek a fit for the cosine weighted BRDF, therefore cos(theta)
            gsl_vector_set(gsl_target_vector, row_idx, eval_target * cosf(theta) * weight);

//            const float u = ((float)j + 0.5f) / (float)NUM_SAMPLES;
//            const float v = ((float)i + 0.5f) / (float)NUM_SAMPLES;
//            glm::vec2 uv(u, v);
//
//            // sample LTSF
//            {
//                glm::vec3 V = sample(uv);
//                float pdf_target;
//                float eval_target = m_target_function->eval(V, pdf_target);
//                float weight = pdf_ltsf / (pdf_target + pdf_ltsf);
//
//                for (int column_idx = 0; column_idx < columns; column_idx++) {
//                    float eval_ltsf = eval_ltsf_basis(V, column_idx);
//                    gsl_matrix_set(gsl_mat, row_idx, column_idx, eval_ltsf * weight);
//                }
//                // we seek a fit for the cosine weighted BRDF, therefore V.z
//                gsl_vector_set(gsl_target_vector, row_idx, eval_target * fmax(0.f, V.z) * weight);
//            }
//
//            row_idx++;
//
//            // sample BRDF
//            {
//                glm::vec3 V = m_target_function->sample(uv);
//                float pdf_target;
//                float eval_brdf = m_target_function->eval(V, pdf_target);
//                float weight = pdf_target / (pdf_target + pdf_ltsf);
//
//                for (int column_idx = 0; column_idx < columns; column_idx++) {
//                    float eval_ltsf = eval_ltsf_basis(V, column_idx);
//                    gsl_matrix_set(gsl_mat, row_idx, column_idx, eval_ltsf * weight);
//                }
//                // we seek a fit for the cosine weighted BRDF, therefore V.z
//                gsl_vector_set(gsl_target_vector, row_idx, eval_brdf * fmax(0.f, V.z) * weight);
//            }
//            row_idx++;
        }
    }

    // initialise gsl linear least squares solver and solve
    auto workspace = gsl_multifit_linear_alloc(rows, columns);
    double chi_squared;
    gsl_multifit_linear(gsl_mat, gsl_target_vector, gsl_coefficients, gsl_cov, &chi_squared, workspace);
    gsl_multifit_linear_free(workspace);

    // set fitted coefficients to spherical function
    auto coefficients = std::unique_ptr<std::vector<float>>(new std::vector<float>(columns));
    std::cout << "Chi Squared: " << chi_squared << std::endl;
    std::cout << "Coefficients: ";
    for (int column = 0; column < columns; column++) {
        (*coefficients)[column] = gsl_vector_get(gsl_coefficients, column);
        std::cout << (*coefficients)[column] << " ";
    }
    std::cout << std::endl;
    m_spherical_function->setCoefficients(std::move(coefficients));

    // free resource from gsl variables
    gsl_matrix_free(gsl_mat);
    gsl_vector_free(gsl_target_vector);
    gsl_vector_free(gsl_coefficients);
}


// inversely transform vector, calculate jacobian
glm::vec3 LTSF::linearly_transform_vec(const glm::vec3& V, float& jacobian) const {
    glm::vec3 trans_V = m_M_inv * V;
    const float norm = glm::length(trans_V);
    jacobian = m_det_M_inv / (norm * norm * norm);
    return glm::normalize(trans_V);
}
