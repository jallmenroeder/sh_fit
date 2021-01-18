//
// Created by jallmenroeder on 13.01.21.
//

#include "ThreadPool.h"

#include "SphericalFunctions/LTSF.h"

#include "Numpy.h"

#include <thread>


ThreadPool::ThreadPool(int LUT_dimension)
: m_LUT_DIMENSION(LUT_dimension),
  m_NUM_THREADS(std::thread::hardware_concurrency()) {
//  m_NUM_THREADS(1) {
    printf("number of threads: %d\n", m_NUM_THREADS);
}


void ThreadPool::execute(const SphericalFunction& spherical_function) {
    // initialize fitting data
    m_matrices = std::make_unique<float[]>(m_LUT_DIMENSION * m_LUT_DIMENSION * 3 * 3);
    m_inv_matrices = std::make_unique<float[]>(m_LUT_DIMENSION * m_LUT_DIMENSION * 4);
    m_coefficients = std::make_unique<float[]>(m_LUT_DIMENSION * m_LUT_DIMENSION * spherical_function.numCoefficients());

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

    aoba::SaveArrayAsNumpy("cos_mat.npy", m_LUT_DIMENSION, m_LUT_DIMENSION, 3, 3, m_matrices.get());
    aoba::SaveArrayAsNumpy("inv_cos_mat.npy", m_LUT_DIMENSION, m_LUT_DIMENSION, 4, m_inv_matrices.get());
    aoba::SaveArrayAsNumpy("cos_coeff.npy", m_LUT_DIMENSION, m_LUT_DIMENSION, spherical_function.numCoefficients(), m_coefficients.get());

    printf("Finished execution\n");
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

            auto mat = ltsf->getLinearTransformation();
            for (int column = 0; column < 3; column++) {
                for (int row = 0; row < 3; row++) {
                    m_matrices[
                            m_LUT_DIMENSION * 3 * 3 * idx.view
                            + 3 * 3 * idx.roughness
                            + 3 * column
                            + row] = (*mat)[column][row];
//                              + row] = (*mat)[row][column];
                }
            }

            auto inv_mat = ltsf->getInvLinearTransformation();
            // rescale the matrix so the first entry is always 1.
            (*inv_mat) = (*inv_mat) / (*inv_mat)[0][0];
            int column = 2; int row = 0;
            for (int i = 0; i < 4; i++) {
                m_inv_matrices[
                        m_LUT_DIMENSION * 4 * idx.view
                        + 4 * idx.roughness
                        + i] = (*inv_mat)[column][row];
                column += 2;
                if (column > 2) {
                    column = column % 3;
                    row++;
                }
            }

            auto coefficients = ltsf->getCoefficients();
            for (int i = 0; i < coefficients->size(); i++) {
                m_coefficients[
                        m_LUT_DIMENSION * coefficients->size() * idx.view
                        + coefficients->size() * idx.roughness
                        + i] = (*coefficients)[i];
            }
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
