cmake_minimum_required(VERSION 3.17)
project(sh_fit_cpp)

set(CMAKE_CXX_STANDARD 11)

FIND_PACKAGE( Boost REQUIRED )
INCLUDE_DIRECTORIES( "/usr/include/boost" )
FIND_PACKAGE(GSL REQUIRED)    # See below (2)

add_executable(sh_fit_cpp
        main.cpp
        Util.h
        SphericalFunctions/GGX_BRDF.h
        SphericalFunctions/GGX_BRDF.cpp
        SphericalFunctions/SphericalHarmonics.h
        SphericalFunctions/SphericalHarmonics.cpp
        SphericalFunctions/SphericalFunction.h
        SphericalFunctions/ClampedCosine.h
        SphericalFunctions/ClampedCosine.cpp)

TARGET_LINK_LIBRARIES( sh_fit_cpp LINK_PUBLIC ${Boost_LIBRARIES} GSL::gsl)