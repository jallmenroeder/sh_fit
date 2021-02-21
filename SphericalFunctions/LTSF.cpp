//
// Created by jallmenroeder on 07/12/2020.
//

#include "LTSF.h"

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
    m_view_dir = sphericalToCartesian({((float)idx.view) / (float)(LUT_dimension - 1) / 2.f * M_PIf32, 0.f});
    m_roughness = (float)idx.roughness / (float)(LUT_dimension - 1);
    m_roughness *= m_roughness;
    m_target_function = std::make_unique<GGX_BRDF>(m_view_dir, m_roughness);

    // init multifit
	int rows = NUM_SAMPLES * NUM_SAMPLES * 2;
	int columns = m_spherical_function->numCoefficients();

	m_gsl_mat = gsl_matrix_alloc(rows, columns);
	m_gsl_cov = gsl_matrix_alloc(columns, columns);
	m_gsl_target_vector = gsl_vector_alloc(rows);
	m_gsl_coefficients = gsl_vector_alloc(columns);
	m_workspace = gsl_multifit_linear_alloc(rows, columns);

	// init multimin
	m_multimin_type = gsl_multimin_fminimizer_nmsimplex2;
	m_multimin_x = gsl_vector_alloc(4);
	m_multimin_step_size = gsl_vector_alloc(4);
	m_multimin_workspace = gsl_multimin_fminimizer_alloc(m_multimin_type, 4);
}


LTSF::~LTSF() {
	// free multifit
	gsl_matrix_free(m_gsl_mat);
	gsl_matrix_free(m_gsl_cov);
	gsl_vector_free(m_gsl_target_vector);
	gsl_vector_free(m_gsl_coefficients);
	gsl_multifit_linear_free(m_workspace);

	// free multimin
	gsl_vector_free(m_multimin_x);
	gsl_vector_free(m_multimin_step_size);
	gsl_multimin_fminimizer_free(m_multimin_workspace);
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
        if (fabs(mapped_uv.x) > fabs(mapped_uv.y)) {
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
    int columns = m_spherical_function->numCoefficients();
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
                    gsl_matrix_set(m_gsl_mat, row_idx, column_idx, eval_ltsf * weight);
                }
                gsl_vector_set(m_gsl_target_vector, row_idx, eval_target * weight);
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
                    gsl_matrix_set(m_gsl_mat, row_idx, column_idx, eval_ltsf * weight);
                }
                gsl_vector_set(m_gsl_target_vector, row_idx, eval_target * weight);
            }
            row_idx++;
        }
    }

    // initialise gsl linear least squares solver and solve
    double chi_squared;
    gsl_multifit_linear(m_gsl_mat, m_gsl_target_vector, m_gsl_coefficients, m_gsl_cov, &chi_squared, m_workspace);

    // set fitted coefficients to spherical function
    auto coefficients = std::make_unique<std::vector<float>>(columns);
    for (int column = 0; column < columns; column++) {
        (*coefficients)[column] = gsl_vector_get(m_gsl_coefficients, column);
    }
    m_spherical_function->setCoefficients(std::move(coefficients));

    return chi_squared;
}


double LTSF::minimizeFunc(const gsl_vector* x, void* params) {
    auto ltsf = static_cast<LTSF*>(params);
    glm::mat3 M(gsl_vector_get(x, 0), 0.f, gsl_vector_get(x, 1),
                0.f, gsl_vector_get(x, 2), 0.f,
                gsl_vector_get(x, 3), 0.f, 1.f);
    ltsf->update(M);
    auto error = ltsf->findSphericalExpansion();
    return error != error ? std::numeric_limits<double>::max() : error;
}


void LTSF::findFit() {
    gsl_multimin_function minex_func;

    size_t iter = 0;
    int status;
    double size;

    // Starting point
    gsl_vector_set(m_multimin_x, 0, (*m_M)[0][0]);
    gsl_vector_set(m_multimin_x, 1, (*m_M)[0][2]);
    gsl_vector_set(m_multimin_x, 2, (*m_M)[1][1]);
    gsl_vector_set(m_multimin_x, 3, (*m_M)[2][0]);

    // Set initial step sizes to 1
    gsl_vector_set_all(m_multimin_step_size, 0.05f);

    // Initialize method and iterate
    minex_func.n = 4;
    minex_func.f = minimizeFunc;
    minex_func.params = static_cast<void*>(this);

    gsl_multimin_fminimizer_set(m_multimin_workspace, &minex_func, m_multimin_x, m_multimin_step_size);

    do
    {
        iter++;
        status = gsl_multimin_fminimizer_iterate(m_multimin_workspace);

        if (status)
            break;

        size = gsl_multimin_fminimizer_size(m_multimin_workspace);
        status = gsl_multimin_test_size(size, 1e-7);
    }
    while (status == GSL_CONTINUE && iter < 200);

    if (status == GSL_SUCCESS) {
		printf("Found fit for idx: v_%d, r_%d\n", m_idx.view, m_idx.roughness);
    } else {
    	printf("Error, multimin did not converge for idx: v_%d, r_%d\n", m_idx.view, m_idx.roughness);
    }
}


// inversely transform vector, calculate jacobian
glm::vec3 LTSF::linearly_transform_vec(const glm::vec3& V, float& jacobian) const {
    glm::vec3 trans_V = *m_M_inv * V;
    const float norm = glm::length(trans_V);
    jacobian = m_det_M_inv / (norm * norm * norm);
    return glm::normalize(trans_V);
}
