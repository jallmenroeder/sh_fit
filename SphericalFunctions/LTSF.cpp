//
// Created by jallmenroeder on 07/12/2020.
//

#include "LTSF.h"

#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multimin.h>
#include <iostream>

#include "../Util.h"
#include "GGX_BRDF.h"


LTSF::LTSF(std::unique_ptr<SphericalFunction> spherical_function, const glm::mat3& M, Idx idx, int LUT_dimension) :
        m_spherical_function(std::move(spherical_function)),
        m_M(std::make_shared<glm::mat3>(M)),
        m_M_inv(std::make_shared<glm::mat3>(glm::inverse(M))),
        m_det_M_inv(glm::determinant(*m_M_inv)),
        m_idx(idx),
        m_error_resolution(2 * NUM_SAMPLES) {
    m_view_dir = sphericalToCartesian({((float)idx.view + .5f) / (float)LUT_dimension / 2.f * M_PIf32, 0.f});
    m_roughness = ((float)idx.roughness + 5.f ) / (float)LUT_dimension;
    m_roughness *= m_roughness;
    m_target_function = std::make_unique<GGX_BRDF>(m_view_dir, m_roughness);
}


void LTSF::update(const glm::mat3 &M) {
    if (abs(glm::determinant(M)) < FLT_EPSILON) {
        *m_M = glm::mat3(1.f);
    } else {
        *m_M = M;
    }
    *m_M_inv = glm::inverse(*m_M);
    m_det_M_inv = glm::determinant(*m_M_inv);
}


void LTSF::setErrorResolution(int resolution) {
    m_error_resolution = resolution;
}


std::unique_ptr<SphericalFunction> LTSF::getSphericalFunctionCopy() const {
    return m_spherical_function->copy();
}


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


double LTSF::findSphericalExpansion() {
    int rows = NUM_SAMPLES * NUM_SAMPLES * 2;
    int columns = m_spherical_function->numCoefficients();
    auto gsl_mat = gsl_matrix_alloc(rows, columns);
    auto gsl_cov = gsl_matrix_alloc(columns, columns);
    auto gsl_target_vector = gsl_vector_alloc(rows);
    auto gsl_coefficients = gsl_vector_alloc(columns);
    int row_idx = 0;

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
                gsl_vector_set(gsl_target_vector, row_idx, eval_target * weight);
            }

            row_idx++;

            // sample BRDF
            {
                glm::vec3 V = m_target_function->sample(uv);
                float pdf_target = m_target_function->pdf(V);
                float eval_target = m_target_function->eval(V);
                float pdf_ltsf = pdf(V);
                float weight = pdf_target / (pdf_target + pdf_ltsf);

                for (int column_idx = 0; column_idx < columns; column_idx++) {
                    float eval_ltsf = evalLtsfBasis(V, column_idx);
                    gsl_matrix_set(gsl_mat, row_idx, column_idx, eval_ltsf * weight);
                }
                // we seek a fit for the cosine weighted BRDF, therefore V.z
                gsl_vector_set(gsl_target_vector, row_idx, eval_target * weight);
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
    auto coefficients = std::make_unique<std::vector<float>>(columns);
    for (int column = 0; column < columns; column++) {
        (*coefficients)[column] = gsl_vector_get(gsl_coefficients, column);
    }
    m_spherical_function->setCoefficients(std::move(coefficients));

    // free resource from gsl variables
    gsl_matrix_free(gsl_mat);
    gsl_vector_free(gsl_target_vector);
    gsl_vector_free(gsl_coefficients);

    return chi_squared;
}


// calculates mean squared error using multiple importance sampling
double LTSF::calculateError() const {
    double error = 0.;
    for (int i = 0; i < m_error_resolution; i++) {
        for (int j = 0; j < m_error_resolution; j++) {
            const float u = ((float)j + 0.5f) / (float)m_error_resolution;
            const float v = ((float)i + 0.5f) / (float)m_error_resolution;
            glm::vec2 uv(u, v);

            // sample LTSF
            {
                glm::vec3 V = sample(uv);
                float pdf_target = m_target_function->pdf(V);
                float eval_target = m_target_function->eval(V);
                float pdf_ltsf = pdf(V);
                float weight = pdf_ltsf / (pdf_target + pdf_ltsf);
                double e = abs(eval_target - eval(V));
                e *= e * e;
                error += e * weight;
            }

            // sample BRDF
            {
                glm::vec3 V = m_target_function->sample(uv);
                float pdf_target = m_target_function->pdf(V);
                float eval_target = m_target_function->eval(V);
                float pdf_ltsf = pdf(V);
                float weight = pdf_target / (pdf_target + pdf_ltsf);
                double e = abs(eval_target - eval(V));
                e *= e * e;
                error += e * weight;
            }
        }
    }
    error /= m_error_resolution * m_error_resolution;
    return error;
}


double LTSF::minimizeFunc(const gsl_vector* x, void* params) {
    auto ltsf = static_cast<LTSF*>(params);
    glm::mat3 M(gsl_vector_get(x, 0), 0.f, gsl_vector_get(x, 1),
                0.f, gsl_vector_get(x, 2), 0.f,
                gsl_vector_get(x, 3), 0.f, 1.f);
    ltsf->update(M);
    // TODO: decide if chi-squared is sufficient as cost function or if error needs to be computed
//    ltsf->findSphericalExpansion();
    auto error = ltsf->findSphericalExpansion();
//    auto error = ltsf->calculateError();
    // check for NaN, x != x only true for NaNs
    return error != error ? std::numeric_limits<double>::max() : error;
}


void LTSF::findFit() {
    const gsl_multimin_fminimizer_type *T =
            gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *minimizer_workspace = nullptr;
    gsl_vector *step_size, *x;
    gsl_multimin_function minex_func;

    size_t iter = 0;
    int status;
    double size;

    // Starting point
    x = gsl_vector_alloc(4);
    gsl_vector_set(x, 0, (*m_M)[0][0]);
    gsl_vector_set(x, 1, (*m_M)[0][2]);
    gsl_vector_set(x, 2, (*m_M)[1][1]);
    gsl_vector_set(x, 3, (*m_M)[2][0]);

    // Set initial step sizes to 1
    step_size = gsl_vector_alloc(4);
    gsl_vector_set_all(step_size, 0.05f);

    // Initialize method and iterate
    minex_func.n = 4;
    minex_func.f = minimizeFunc;
    minex_func.params = static_cast<void*>(this);

    minimizer_workspace = gsl_multimin_fminimizer_alloc(T, 4);
    gsl_multimin_fminimizer_set(minimizer_workspace, &minex_func, x, step_size);

    do
    {
        iter++;
        status = gsl_multimin_fminimizer_iterate(minimizer_workspace);

        if (status)
            break;

        size = gsl_multimin_fminimizer_size(minimizer_workspace);
        status = gsl_multimin_test_size(size, 1e-2);

        if (status == GSL_SUCCESS)
        {
//            printf ("converged to minimum at\n");
        }

//        printf ("%5d %10.3e %10.3e %10.3e %10.3e f() = %7.3f size = %.3f\n",
//                iter,
//                gsl_vector_get(minimizer_workspace->x, 0),
//                gsl_vector_get(minimizer_workspace->x, 1),
//                gsl_vector_get(minimizer_workspace->x, 2),
//                gsl_vector_get(minimizer_workspace->x, 3),
//                minimizer_workspace->fval, size);
    }
    while (status == GSL_CONTINUE && iter < 100);

    printf("Found fit for idx: v_%d, r_%d\n", m_idx.view, m_idx.roughness);
    gsl_vector_free(x);
    gsl_vector_free(step_size);
    gsl_multimin_fminimizer_free(minimizer_workspace);
}


// inversely transform vector, calculate jacobian
glm::vec3 LTSF::linearly_transform_vec(const glm::vec3& V, float& jacobian) const {
    glm::vec3 trans_V = *m_M_inv * V;
    const float norm = glm::length(trans_V);
    jacobian = m_det_M_inv / (norm * norm * norm);
    return glm::normalize(trans_V);
}
