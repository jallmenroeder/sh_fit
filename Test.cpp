//
// Created by jallmenroeder on 15.02.21.
//

#include "Test.h"

#include <iostream>
#include <random>

#include "Util.h"

const int NUM_TEST_SAMPLES = 1000000;

bool Test::testBrdfSampling(BRDF* brdf) {
	double uniform_result = 0.;
	double importance_result = 0.;

	// init RNG
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<float> dist(0.0, 1.0);

	for (int i = 0; i < NUM_TEST_SAMPLES; i++) {
		// uniform sampling
		{
			glm::vec2 uv(dist(mt), dist(mt));
			Spherical _v{};
			_v.theta = acosf(1.f - uv.x);
			_v.phi = uv.y * 2.f * M_PIf32;
			double weight = 2. * M_PIf32;
			auto v = sphericalToCartesian({_v.theta, _v.phi});
			auto eval = brdf->eval(v);
			uniform_result += eval * weight;
		}
		// importance sampling
		{
			glm::vec2 uv(dist(mt), dist(mt));
			auto v = brdf->sample(uv);
			float density = brdf->pdf(v);
			double weight = density > FLT_EPSILON ? 1. / density : 0.;
			auto eval = brdf->eval(v);
			importance_result += eval * weight;
		}
	}
	uniform_result /= (double)NUM_TEST_SAMPLES;
	importance_result /= (double)NUM_TEST_SAMPLES;
	double error = fabs(uniform_result - importance_result);
	printf("uniform: %f\n", uniform_result);
	printf("IS: %f\n", importance_result);
	printf("error: %f\n", error);
	return error < 0.001;
}