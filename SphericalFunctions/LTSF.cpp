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


float LTSF::eval(const glm::vec3& V) const {
    // evaluate original spherical function
    float jacobian;
    glm::vec3 trans_V = linearly_transform_vec(V, jacobian);
    return m_spherical_function->eval(trans_V) * jacobian;
}


float LTSF::pdf(const glm::vec3& V) const {
    return V.z / M_PIf32;
}


float LTSF::evalLtsfBasis(const glm::vec3& V, int idx) const {
    float jacobian;
    glm::vec3 trans_V = linearly_transform_vec(V, jacobian);
    return m_spherical_function->eval_basis(trans_V, idx) * jacobian;
}

// cosine weighted sampling on the hemisphere
// see http://www.pbr-book.org/3ed-2018/Monte_Carlo_Integration/2D_Sampling_with_Multidimensional_Transformations.html#Cosine-WeightedHemisphereSampling
glm::vec3 LTSF::sample(const glm::vec2& uv) const {
    // map uv to [-1,1]
    glm::vec2 mapped_uv = 2.f * uv - glm::vec2(1.f, 1.f);
    glm::vec2 disc_sample;
    // handle degeneracy at origin
    if (mapped_uv.x == 0.f && mapped_uv.y == 0.f) {
        disc_sample = glm::vec2(0.f, 0.f);
    } else {
        float theta, r;
        if (abs(mapped_uv.x) > abs(mapped_uv.y)) {
            r = mapped_uv.x;
            theta = M_PIf32 * (mapped_uv.y / mapped_uv.x) / 4.f;
        } else {
            r = mapped_uv.y;
            theta = M_PIf32 / 2.f - M_PIf32 / 4.f * (mapped_uv.x / mapped_uv.y);
        }
        disc_sample = r * glm::vec2(std::cos(theta), std::sin(theta));
    }
    float z = std::sqrt(std::max(0.f, 1.f - disc_sample.x * disc_sample.x - disc_sample.y * disc_sample.y));
    return {disc_sample.x, disc_sample.y, z};
}


void LTSF::findFit() {
    int rows = NUM_SAMPLES * NUM_SAMPLES * 2;
    int columns = m_spherical_function->numCoefficients();
    auto gsl_mat = gsl_matrix_alloc(rows, columns);
    auto gsl_cov = gsl_matrix_alloc(columns, columns);
    auto gsl_target_vector = gsl_vector_alloc(rows);
    auto gsl_coefficients = gsl_vector_alloc(columns);
    int row_idx = 0;

    std::vector<float> jacobian(rows);
    for (int i = 0; i < NUM_SAMPLES; i++) {
        for (int j = 0; j < NUM_SAMPLES; j++) {
            const float u = ((float)j + 0.5f) / (float)NUM_SAMPLES;
            const float v = ((float)i + 0.5f) / (float)NUM_SAMPLES;
            glm::vec2 uv(u, v);

            // sample LTSF
            {
                glm::vec3 V = sample(uv);
                float pdf_target = m_target_function->pdf(V);
                float eval_target = m_target_function->eval(V);
                float pdf_ltsf = pdf(V);
                float weight = pdf_ltsf / (pdf_target + pdf_ltsf);

                for (int column_idx = 0; column_idx < columns; column_idx++) {
                    float eval_ltsf = evalLtsfBasis(V, column_idx);
                    gsl_matrix_set(gsl_mat, row_idx, column_idx, eval_ltsf * weight);
                }
                // we seek a fit for the cosine weighted BRDF, therefore V.z
                gsl_vector_set(gsl_target_vector, row_idx, eval_target * fmax(0.f, V.z) * weight);
            }

            row_idx++;

            // sample BRDF
            {
                glm::vec3 V = m_target_function->sample(uv);
                float pdf_target = m_target_function->pdf(V);
                float eval_brdf = m_target_function->eval(V);
                float pdf_ltsf = pdf(V);
                float weight = pdf_target / (pdf_target + pdf_ltsf);

                for (int column_idx = 0; column_idx < columns; column_idx++) {
                    float eval_ltsf = evalLtsfBasis(V, column_idx);
                    gsl_matrix_set(gsl_mat, row_idx, column_idx, eval_ltsf * weight);
                }
                // we seek a fit for the cosine weighted BRDF, therefore V.z
                gsl_vector_set(gsl_target_vector, row_idx, eval_brdf * fmax(0.f, V.z) * weight);
            }
            row_idx++;
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


// calculates mean squared error over hemisphere
double LTSF::calculateError(int resolution) const {
    double error = 0.;
    for (int i = 0; i < resolution; i++) {
        for (int j = 0; j < 4 * resolution; j++) {
            float theta = ((float)i + .5f) / (float)resolution / 2.f * M_PIf32;
            float phi = ((float)j + .5f) / (float)resolution / 4.f / 2.f * M_PIf32;

            auto V = sphericalToCartesian({theta, phi});
            double current_error = eval(V) - m_target_function->eval(V);
            current_error *= current_error;
            error += current_error * sin(theta);
        }
    }
    return sqrt(error / (4. * (double)resolution * (double)resolution));
}


// inversely transform vector, calculate jacobian
glm::vec3 LTSF::linearly_transform_vec(const glm::vec3& V, float& jacobian) const {
    glm::vec3 trans_V = m_M_inv * V;
    const float norm = glm::length(trans_V);
    jacobian = m_det_M_inv / (norm * norm * norm);
    return glm::normalize(trans_V);
}
