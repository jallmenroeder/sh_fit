#include <iostream>
#include <memory>

#include "SphericalFunctions/LTSF.h"
#include "SphericalFunctions/SphericalHarmonics.h"
#include "SphericalFunctions/ClampedCosine.h"
#include "Util.h"
#include "chrono"

int main() {
    auto start = std::chrono::high_resolution_clock::now();
    SphericalFunction* sf = new SphericalHarmonics(4);
//    SphericalFunction* sf = new ClampedCosine(1.f);
    int view_idx = 10;
    int roughness_idx = 10;
    glm::vec3 view_dir = sphericalToCartesian({(float)view_idx / 63.f / 2.f * M_PIf32, 0.f});
    float roughness = ((float)roughness_idx / 63.f);
    roughness *= roughness;
    glm::mat3 M(0.68841448, 0.f, 0.37526232,
                0.f, 0.56166123, 0.f,
                -0.96595106, 0.f, 1.f);
//    glm::mat3 M(0.60742463, 0., 0.27948604,
//                0., 0.49071796, 0.,
//                -0.90142154, 0., 1.);
    LTSF ltsf = LTSF(std::unique_ptr<SphericalFunction>(sf), M, view_dir, roughness);
    ltsf.findFit();
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    printf("Duration: %d ms", (int)duration.count());
    return 0;
}
