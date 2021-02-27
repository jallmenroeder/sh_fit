#include <iostream>
#include <memory>

#include "ThreadPool.h"
#include "SphericalFunctions/SphericalHarmonics.h"
#include "SphericalFunctions/ClampedCosine.h"
#include "SphericalFunctions/SphericalFunction.h"
#include "chrono"

int main() {
    auto start = std::chrono::high_resolution_clock::now();

    // this does NOT init 64 threads, it just says our LUT is 64x64
    ThreadPool pool(64);
//    SphericalFunction* sf = new SphericalHarmonics(2);
    SphericalFunction* sf = new ClampedCosine(1.f);
    pool.execute(*sf);
    delete sf;
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    printf("Duration: %f min", duration.count() / 60.f);
    return 0;
}
