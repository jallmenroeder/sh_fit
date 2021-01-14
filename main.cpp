#include <iostream>
#include <memory>

#include "ThreadPool.h"
#include "SphericalFunctions/SphericalHarmonics.h"
#include "SphericalFunctions/ClampedCosine.h"
#include "SphericalFunctions/SphericalFunction.h"
#include "chrono"

int main() {
    auto start = std::chrono::high_resolution_clock::now();

    ThreadPool pool(64);
//    SphericalFunction* sf = new SphericalHarmonics(2);
    SphericalFunction* sf = new ClampedCosine(1.f);
    pool.execute(*sf);
    delete sf;
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    printf("Duration: %d ms", (int)duration.count());
    return 0;
}
