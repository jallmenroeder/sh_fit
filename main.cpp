#include <iostream>
#include <memory>

#include "ThreadPool.h"
#include "SphericalFunctions/SphericalHarmonics.h"
#include "SphericalFunctions/ClampedCosine.h"
#include "SphericalFunctions/SphericalFunction.h"
#include "chrono"

int main() {
    auto start = std::chrono::high_resolution_clock::now();

    ThreadPool pool;
//    SphericalFunction* sf = new SphericalHarmonics(2);
    SphericalFunction* sf = new ClampedCosine();
    pool.execute(*sf);
    delete sf;
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    printf("Duration: %f min", duration.count() / 60.f);
    return 0;
}
