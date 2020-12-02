cmake_minimum_required(VERSION 3.17)
project(sh_fit_cpp)

set(CMAKE_CXX_STANDARD 11)

FIND_PACKAGE( Boost REQUIRED )
INCLUDE_DIRECTORIES( "/usr/include/boost" )

add_executable(sh_fit_cpp
        main.cpp
        BRDF.h
        Util.h
        SphericalHarmonics.h
        SphericalHarmonics.cpp
        SphericalFunction.h
        ClampedCosine.h)

TARGET_LINK_LIBRARIES( sh_fit_cpp LINK_PUBLIC ${Boost_LIBRARIES} )