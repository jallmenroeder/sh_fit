cmake_minimum_required(VERSION 3.17)
project(sh_fit_cpp)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS_RELEASE "-O3")

FIND_PACKAGE(GSL REQUIRED)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

add_executable(sh_fit_cpp
        main.cpp
        Util.h
        SphericalFunctions/GGX_BRDF.h
        SphericalFunctions/GGX_BRDF.cpp
        SphericalFunctions/SphericalHarmonics.h
        SphericalFunctions/SphericalHarmonics.cpp
        SphericalFunctions/SphericalFunction.h
        SphericalFunctions/ClampedCosine.h
        SphericalFunctions/ClampedCosine.cpp
        SphericalFunctions/BRDF.h
        SphericalFunctions/LTSF.cpp
        SphericalFunctions/LTSF.h
        ThreadPool.cpp
        ThreadPool.h
        Test.cpp
        Test.h
        Numpy.h)

TARGET_LINK_LIBRARIES( sh_fit_cpp LINK_PUBLIC GSL::gsl Threads::Threads)