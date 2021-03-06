//
// Created by jallmenroeder on 24/11/2020.
//

#pragma once

#include <glm/glm.hpp>
#include <memory>

template<class T> using uptr = std::unique_ptr<T>;

struct Spherical {
    float theta, phi;
};

static glm::vec3 sphericalToCartesian(const Spherical& in) {
    return {
            sinf(in.theta) * cosf(in.phi),
            sinf(in.theta) * sinf(in.phi),
            cosf(in.theta)
    };
}


static Spherical cartesianToSpherical(const glm::vec3& cartesian) {
    return {acosf(cartesian.z), atan2f(cartesian.y, cartesian.x)};
}


static float cdot(const glm::vec3& a, const glm::vec3& b) {
    return fmaxf(0.f, glm::dot(a, b));
}


static int factorial(int i) {
	float res = 1;
	for (; i > 0; i--) {
		res *= i;
	}
	return res;
}
