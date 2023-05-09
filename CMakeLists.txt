cmake_minimum_required(VERSION 3.17)
project(sh_fit_cpp)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS_RELEASE "-O3")

include_directories(${PROJECT_SOURCE_DIR})

FIND_PACKAGE(GSL REQUIRED)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

add_subdirectory(matplotplusplus)

add_executable(sh_fit_cpp
        main.cpp
        Util.h
        SphericalFunctions/GGX_BRDF.h
        SphericalFunctions/SphericalHarmonics.h
        SphericalFunctions/SphericalFunction.h
        SphericalFunctions/ClampedCosine.h
        SphericalFunctions/BRDF.h
        SphericalFunctions/LTSF.cpp
        SphericalFunctions/LTSF.h
        plotting/Plotter.cpp
        plotting/Plotter.h
        ThreadPool.cpp
        ThreadPool.h
        Test.h
        Numpy.h)

TARGET_LINK_LIBRARIES(sh_fit_cpp LINK_PUBLIC GSL::gsl Threads::Threads matplot)