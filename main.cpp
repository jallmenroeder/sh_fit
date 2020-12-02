#include <iostream>
#include <glm/glm.hpp>

#include "Util.h"
#include "BRDF.h"
#include "SphericalHarmonics.h"

int main() {
    std::cout << "Hello, World!" << std::endl;
    glm::vec3 cartesian = {1.f, 1.f, 1.f};
    cartesian = glm::normalize(cartesian);
    Spherical spherical = cartesianToSpherical(cartesian);
    std::cout << spherical.theta << ", " << spherical.phi << std::endl;
    cartesian = sphericalToCartesian(spherical);
    std::cout << cartesian.x << ", " << cartesian.y << ", " << cartesian.z << std::endl;
    float pdf = 0.f;
    glm::vec3 a{3.f, 0.5f, 0.5f};
    glm::vec3 b{0.f, -0.5f, 0.5f};
    glm::vec3 c{0.f, -0.5f, 0.8f};
    a = glm::normalize(a);
    b = glm::normalize(b);
    c = glm::normalize(c);
    std::cout << "GGX: " << eval_ggx(a, b, .2, pdf) << std::endl;
    std::cout << "GGX: " << eval_ggx(a, c, .2, pdf) << std::endl;
    return 0;
}
