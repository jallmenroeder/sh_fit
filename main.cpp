#include <iostream>
#include <glm/glm.hpp>

#include "Util.h"
#include "SphericalFunctions/BRDF.h"

int main() {
    std::cout << "Hello, World!" << std::endl;
    glm::vec3 cartesian = {1.f, 1.f, 1.f};
    cartesian = glm::normalize(cartesian);
    Spherical spherical = cartesianToSpherical(cartesian);
    std::cout << spherical.theta << ", " << spherical.phi << std::endl;
    cartesian = sphericalToCartesian(spherical);
    std::cout << cartesian.x << ", " << cartesian.y << ", " << cartesian.z << std::endl;
    glm::vec3 a{3.f, 0.5f, 0.5f};
    glm::vec3 b{0.f, -0.5f, 0.5f};
    glm::vec3 c{0.f, -0.5f, 0.8f};
    a = glm::normalize(a);
    b = glm::normalize(b);
    c = glm::normalize(c);
    BRDF ggx(a, .2);
    float pdf = 0.f;
    std::cout << "GGX: " << ggx.eval(b, pdf) << std::endl;
    std::cout << "GGX: " << ggx.eval(c, pdf) << std::endl;
    return 0;
}
