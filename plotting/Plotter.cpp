//
// Created by jallmenroeder on 09.05.23.
//

#include "Plotter.h"

#include <matplot/matplot.h>

void Plotter::plotResidual(const float* residual, int lutDimension, float max) {
    std::vector<std::vector<double>> data(lutDimension);
    for (int i = 0; i < lutDimension; i++) {
        data[i] = std::vector<double>(lutDimension);

        for (int j = 0; j < lutDimension; j++) {
            float dataPoint = residual[i * lutDimension + j];
            if (max > 0.)
                dataPoint = std::fmin(max, dataPoint);
            data[i][j] = dataPoint;
        }
    }
    matplot::heatmap(data);
    matplot::show();
}

void Plotter::plotSurfacePoints(const std::vector<glm::vec3>& points) {
    std::vector<double> x, y, z;
    for (const auto& point: points) {
        x.emplace_back(point.x);
        y.emplace_back(point.y);
        z.emplace_back(point.z);
    }
    matplot::scatter3(x, y, z);
    matplot::show();
}

void Plotter::plotSphericalFunction(const SphericalFunction& sphericalFunction) {

}

void Plotter::imagePlot(const std::vector<std::vector<float>>& data) {
    matplot::image(data);
    matplot::show();
}
