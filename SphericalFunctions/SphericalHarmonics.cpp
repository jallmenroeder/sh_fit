//
// Created by jallmenroeder on 30/11/2020.
//

#include "SphericalHarmonics.h"

#include "../Util.h"

#include "gsl/gsl_sf_legendre.h"

SphericalHarmonics::SphericalHarmonics(int order)
        : m_ORDER(order),
          m_NUM_COEFFS((order + 1) * (order + 1)),
          m_LEGENDRE_SIZE((order + 1) * (order + 2) / 2),
          m_coeffs(new std::vector<float>(m_NUM_COEFFS, 0.f)) {

	// normalization for SH basis functions
	m_normalization.resize(m_NUM_COEFFS);
	for (int i = 0; i < m_NUM_COEFFS; i++) {
		int l = floor(sqrt(i) + FLT_EPSILON);
		int m = i - l * (l + 1);
		if (m == 0) {
			m_normalization[i] = K_l_m(l, m);
		} else {
			m_normalization[i] = sqrtf(2.f) * K_l_m(l, m);
		}
	}
}


std::unique_ptr<SphericalFunction> SphericalHarmonics::copy() const {
    return std::make_unique<SphericalHarmonics>(m_ORDER);
}


float SphericalHarmonics::eval(const glm::vec3& V) const {
    float sum = 0.f;
    auto eval_array = std::vector<float>(m_NUM_COEFFS);
    evalBasisArray(V, eval_array);
    for (int coeff_idx = 0; coeff_idx < m_NUM_COEFFS; coeff_idx++) {
        sum += (*m_coeffs)[coeff_idx] * eval_array[coeff_idx];
    }
    return sum;
}


void SphericalHarmonics::evalBasisArray(const glm::vec3& V, std::vector<float> &array) const {
	assert(array.size() >= m_NUM_COEFFS);
	double legendre[m_LEGENDRE_SIZE];
	int status = gsl_sf_legendre_array(GSL_SF_LEGENDRE_NONE, m_ORDER, fmin(V.z, 1.f), legendre);
	assert(status == 0);

	for (int i = 0; i < m_NUM_COEFFS; i++) {
		int l = floor(sqrt(i) + FLT_EPSILON);
		int m = i - l * (l + 1);
		int legendre_idx = gsl_sf_legendre_array_index(l, abs(m));
		// needed due to different SH convention
		float sign = i % 2 ? -1.f : 1.f;
		if (m > 0) {
			array[i] = sign * m_normalization[i] * cosf((float)m * atan2f(V.y, V.x)) * legendre[legendre_idx];
		} else if (m < 0) {
			array[i] = sign * m_normalization[i] * sinf((float)(-m) * atan2f(V.y, V.x)) * legendre[legendre_idx];
		} else {
			array[i] = sign * m_normalization[i] * legendre[legendre_idx];
		}
	}
}


void SphericalHarmonics::setCoefficients(std::shared_ptr<std::vector<float>> coeffs) {
    if (coeffs->size() != m_NUM_COEFFS) {
        printf("Error, SH order was %d, needing %d coefficients but got %ld coefficients, defaulting to 0 instead",
               m_ORDER, m_NUM_COEFFS, coeffs->size());
        return;
    }
    m_coeffs = std::move(coeffs);
}


std::shared_ptr<std::vector<float>> SphericalHarmonics::getCoefficients() {
    return m_coeffs;
}

// calculates the normalization constant K^{l}_{m} for SH
float SphericalHarmonics::K_l_m(int l, int m) {
	m = abs(m);
	return sqrtf((float)((2 * l + 1) * factorial(l - m)) / (4.f * M_PIf32 * (float)factorial(l + m)));
}
