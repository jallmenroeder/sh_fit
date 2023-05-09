//
// Created by jallmenroeder on 09.05.23.
//

#pragma once

#include <vector>
#include <glm/glm.hpp>

#include "SphericalFunctions/SphericalFunction.h"

class Plotter {
public:
    static void plotResidual(const float* residual, int lutDimension, float max = -1.);
    static void plotSurfacePoints(const std::vector<glm::vec3>& points);
    static void plotSphericalFunction(const SphericalFunction& sphericalFunction);
    static void imagePlot(const std::vector<std::vector<float>>& data);
};
