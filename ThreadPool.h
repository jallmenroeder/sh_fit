//
// Created by jallmenroeder on 13.01.21.
//

#pragma once

#include "SphericalFunctions/SphericalFunction.h"
#include "SphericalFunctions/LTSF.h"

#include <queue>
#include <mutex>
#include <memory>
#include <condition_variable>
#include <glm/mat3x3.hpp>

class LTSF;

/**
 * Thread pool class based on this explanation: https://stackoverflow.com/questions/15752659/thread-pooling-in-c11
 */
class ThreadPool {
public:
    ThreadPool();
    explicit ThreadPool(int num_threads);
    void execute(const SphericalFunction& spherical_function);
    void plotResidual(float max);

private:
    void threadLoop();
    void persistData(const LTSF& ltsf);

    std::queue<uptr<LTSF>> m_queue;
    std::mutex m_queue_mutex;

    std::vector<glm::mat3> m_ltc_params;

    uptr<float[]> m_matrices;
    uptr<float[]> m_inv_matrices;
    uptr<float[]> m_coefficients;
    uptr<float[]> m_residual;
    std::mutex m_data_mutex;

    const int m_LUT_DIMENSION;
    const unsigned int m_NUM_THREADS;
};



