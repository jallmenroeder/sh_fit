cmake_minimum_required(VERSION 3.17)
project(sh_fit_cpp)

set(CMAKE_CXX_STANDARD 17)

FIND_PACKAGE( Boost REQUIRED )
INCLUDE_DIRECTORIES( "/usr/include/boost" )

add_executable(sh_fit_cpp main.cpp BRDF.h SphericalHarmonics.h Util.h SphericalHarmonics.cpp)

TARGET_LINK_LIBRARIES( sh_fit_cpp LINK_PUBLIC ${Boost_LIBRARIES} )