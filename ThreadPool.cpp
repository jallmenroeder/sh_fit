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
    m_residual = std::make_unique<float[]>(m_LUT_DIMENSION * m_LUT_DIMENSION);

    // load initial parameters
    auto ltc_matrices = std::vector<double>();
    aoba::LoadArrayFromNumpy("../cos_mat.npy", ltc_matrices);

    // convert loaded data to right format
    m_ltc_params.resize(m_LUT_DIMENSION);
    for (int i = 0; i < m_LUT_DIMENSION; i++) {
        Idx idx = {i, m_LUT_DIMENSION - 1};
        for (int j = 0; j < 3; j++) {
        	for (int k = 0; k < 3; k++) {
        		m_ltc_params[i][j][k] = ltc_matrices[idx.view * m_LUT_DIMENSION * 3 * 3 + idx.roughness * 3 * 3 + j * 3 + k];
        	}
        }
        // add fitting tasks to execution queue
        m_queue.emplace(new LTSF(spherical_function.copy(), m_ltc_params[i], idx, m_LUT_DIMENSION));
    }

    // start thread execution
    std::vector<std::thread> pool;
    for (unsigned int thread_idx = 0; thread_idx < m_NUM_THREADS; thread_idx++) {
        pool.emplace_back(std::thread(&ThreadPool::threadLoop, this));
    }

    for (auto& thread: pool) {
        thread.join();
    }

    // save fitting data as numpy arrays
    auto func_name = spherical_function.getName();
    aoba::SaveArrayAsNumpy(func_name + "_mat.npy", m_LUT_DIMENSION, m_LUT_DIMENSION, 3, 3, m_matrices.get());
    aoba::SaveArrayAsNumpy("inv_" + func_name + "_mat.npy", m_LUT_DIMENSION, m_LUT_DIMENSION, 4, m_inv_matrices.get());
    aoba::SaveArrayAsNumpy(func_name + "_coeff.npy", m_LUT_DIMENSION, m_LUT_DIMENSION, spherical_function.numCoefficients(), m_coefficients.get());
    aoba::SaveArrayAsNumpy(func_name + "_residual.npy", m_LUT_DIMENSION, m_LUT_DIMENSION, m_residual.get());

    printf("Finished execution\n");
}


void ThreadPool::threadLoop() {
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

        persistData(*ltsf);

        // enqueue next fit with lower roughness
        Idx current_idx = ltsf->getIdx();
    	glm::mat3 guess;
    	if (current_idx.roughness <= 0) {
			printf("finished view dir with index: %d\n", current_idx.view);
			continue;
    	} else {
			current_idx.roughness--;
			guess = *ltsf->getLinearTransformation();
    	}
        auto next = std::make_unique<LTSF>(ltsf->getSphericalFunctionCopy(),
                                           guess,
                                           current_idx,
                                           m_LUT_DIMENSION);
        {
            std::unique_lock<std::mutex> lock(m_queue_mutex);
            m_queue.push(std::move(next));
        }
    }
}


void ThreadPool::persistData(const LTSF& ltsf) {
    auto idx = ltsf.getIdx();
    std::unique_lock<std::mutex> lock(m_data_mutex);

    // persist linear transformation M
    auto mat = ltsf.getLinearTransformation();
    for (int column = 0; column < 3; column++) {
        for (int row = 0; row < 3; row++) {
            m_matrices[
                    m_LUT_DIMENSION * 3 * 3 * idx.view
                    + 3 * 3 * idx.roughness
                    + 3 * column
                    + row] = (*mat)[column][row];
        }
    }

    // persist inverse linear transformation M⁻¹
    auto inv_mat = ltsf.getInvLinearTransformation();
    // rescale the matrix so the first entry is always 1.
    (*inv_mat) = (*inv_mat) / (*inv_mat)[0][0];
    int column = 2; int row = 0;
    for (int i = 0; i < 4; i++) {
        m_inv_matrices[
                m_LUT_DIMENSION * 4 * idx.roughness
                + 4 * idx.view
                + i] = (*inv_mat)[row][column];
        column += 2;
        if (column > 2) {
            column = column % 3;
            row++;
        }
    }

    // persist coefficient(s)
    auto coefficients = ltsf.getCoefficients();
    for (int i = 0; i < coefficients->size(); i++) {
        m_coefficients[
                m_LUT_DIMENSION * coefficients->size() * idx.view
                + coefficients->size() * idx.roughness
                + i] = (*coefficients)[i];
    }

    // persist residual
    m_residual[m_LUT_DIMENSION * idx.view + idx.roughness] = ltsf.getResidual();
}
