//
// Created by jallmenroeder on 13.01.21.
//

#include "ThreadPool.h"

#include "SphericalFunctions/LTSF.h"

#include <thread>


ThreadPool::ThreadPool(int LUT_dimension)
: m_LUT_DIMENSION(LUT_dimension),
  m_NUM_THREADS(std::thread::hardware_concurrency()) {
    // initialize fitting data
    m_matrices.resize(LUT_dimension);
    m_coefficients.resize(LUT_dimension);
    for (int i = 0; i < LUT_dimension; i++) {
        m_matrices[i].resize(LUT_dimension);
        m_coefficients[i].resize(LUT_dimension);
    }
    printf("number of threads: %d\n", m_NUM_THREADS);
}


void ThreadPool::execute(const SphericalFunction& spherical_function) {
    // TODO: add LTC parameters as initial guess
    for (int i = 0; i < m_LUT_DIMENSION; i++) {
        Idx idx = {i, m_LUT_DIMENSION / 2 - 1};
        addLTSF(std::make_unique<LTSF>(spherical_function.copy(), glm::mat3(1.f), idx, m_LUT_DIMENSION));
    }

    std::vector<std::thread> pool;
    for (unsigned int thread_idx = 0; thread_idx < m_NUM_THREADS; thread_idx++) {
        pool.emplace_back(std::thread(&ThreadPool::infiniteLoopFunction, this));
    }

    for (auto& thread: pool) {
        thread.join();
    }
    printf("Finished execution");
}


void ThreadPool::infiniteLoopFunction() {
    while(true) {
        // acquire next available LTSF from queue
        std::unique_ptr<LTSF> ltsf;
        {
            std::unique_lock<std::mutex> lock(m_queue_mutex);
            if (m_queue.empty()) break;
            ltsf = std::move(m_queue.front());
            m_queue.pop();
        }
        // find linear transformation and coefficients for this LTSF configuration, this is the actual work
        ltsf->findFit();

        // persist data
        auto idx = ltsf->getIdx();
        {
            std::unique_lock<std::mutex> lock(m_data_mutex);
            m_matrices[idx.view][idx.roughness] = ltsf->getInvLinearTransformation();
            m_coefficients[idx.view][idx.roughness] = ltsf->getCoefficients();
        }

        Idx current_idx = ltsf->getIdx();
        if ((current_idx.roughness < m_LUT_DIMENSION / 2) && (current_idx.roughness > 0)) {
            current_idx.roughness--;
        } else if ((current_idx.roughness >= m_LUT_DIMENSION / 2) && (current_idx.roughness < m_LUT_DIMENSION - 1)) {
            current_idx.roughness++;
        } else if (current_idx.roughness == 0) {
            current_idx.roughness = m_LUT_DIMENSION / 2;
        } else {
            printf("finished view dir with index: %d\n", current_idx.view);
            continue;
        }
        auto next = std::make_unique<LTSF>(ltsf->getSphericalFunctionCopy(),
                                           *ltsf->getLinearTransformation(),
                                           current_idx,
                                           m_LUT_DIMENSION);
        {
            std::unique_lock<std::mutex> lock(m_queue_mutex);
            m_queue.push(std::move(next));
        }
    }
}


void ThreadPool::addLTSF(std::unique_ptr<LTSF> ltsf) {
    std::unique_lock<std::mutex> lock(m_queue_mutex);
    m_queue.push(std::move(ltsf));
}
