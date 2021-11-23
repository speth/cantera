//! @file AdaptivePreconditioner.cpp

// This file is part of Cantera. See License.txt in the top-level
// directory or at https://cantera.org/license.txt for license and
// copyright information.

#include "cantera/numerics/AdaptivePreconditioner.h"
#include "cantera/zeroD/MoleReactor.h"
#include "cantera/zeroD/IdealGasConstPressureMoleReactor.h"
#include "cantera/zeroD/IdealGasMoleReactor.h"
#include "cantera/base/utilities.h"
#include "cantera/base/global.h"
#include <iostream>

namespace Cantera
{

bool AdaptivePreconditioner::operator== (const AdaptivePreconditioner &externalPrecon)
{
    return m_precon_matrix.isApprox(externalPrecon.m_precon_matrix);
}

void AdaptivePreconditioner::operator= (const AdaptivePreconditioner &externalPrecon)
{
    // copy all variables
    m_values = externalPrecon.m_values;
    m_inner = externalPrecon.m_inner;
    m_outer = externalPrecon.m_outer;
    m_sizes = externalPrecon.m_sizes;
    m_starts = externalPrecon.m_starts;
    m_thermo_indices = externalPrecon.m_thermo_indices;
    m_state_indices = externalPrecon.m_state_indices;
    m_nnz = externalPrecon.m_nnz;
    m_identity = externalPrecon.m_identity;
    m_precon_matrix = externalPrecon.m_precon_matrix;
    m_threshold = externalPrecon.m_threshold;
    m_zero = externalPrecon.m_zero;
    m_perturb = externalPrecon.m_perturb;
    m_dimensions = externalPrecon.m_dimensions;
    m_atol = externalPrecon.m_atol;
    m_gamma = externalPrecon.m_gamma;
    m_init = externalPrecon.m_init;
    m_ctr = 0;
    return;
}

double& AdaptivePreconditioner::operator() (size_t row, size_t col)
{
    size_t index = globalIndex(row, col);
    if (m_state_indices.find(index) != m_state_indices.end())
    {
        size_t idx = m_state_indices.at(index);
        return m_values[idx];
    }
    else
    {
        m_zero = 0; // reset value to zero in the event it was changed by assignment
        return m_zero;
    }
}

void AdaptivePreconditioner::initialize(ReactorNet& network)
{
    //! don't use legacy rate constants
    use_legacy_rate_constants(false);
    // reset arrays in case of re-initialization
    m_dimensions.clear();
    m_outer.clear();
    m_inner.clear();
    m_thermo_indices.clear();
    m_values.clear();
    m_sizes.clear();
    m_state_indices.clear();
    m_nnz = 0;
    // set dimensions of preconditioner from network
    size_t totalColLen = network.neq();
    m_dimensions.push_back(totalColLen);
    m_dimensions.push_back(totalColLen);
    if (m_dimensions[0] != m_dimensions[1])
    {
        throw CanteraError("initialize", "specified matrix dimensions are not square");
    }
    // reserve maximum space for vectors making up SparseMatrix
    m_inner.reserve(totalColLen * totalColLen);
    m_outer.reserve(totalColLen * totalColLen);
    m_outer.push_back(0); // first outer starts at zero
    m_thermo_indices.reserve(3 * network.nreactors() * totalColLen);
    // allocate spaces for sizes and set first start to zero
    m_sizes.reserve(network.nreactors() + 1);
    m_sizes.push_back(0);
    // set up sparse patterns
    for (size_t i = 0; i < network.nreactors(); i++)
    {
        Reactor& currReactor = network.reactor(i);
        ReactionDerivativeManager& reaction_derv_mgr = currReactor.reaction_derivative_manager();
        // this is for reinitialization - an error is caused if two
        // preconditioners are initialized with the same network
        if (reaction_derv_mgr.isInit())
        {
            reaction_derv_mgr.reset();
            reaction_derv_mgr.initialize(currReactor);
        }
        size_t currStart = network.reactor_start(i);
        m_starts.push_back(currStart);
        // checking manager
        int column_size = currReactor.neq();
        int species_start = currReactor.species_start();
        // create temporary indices vector and reserve space for it
        vector_int indices;
        std::unordered_map<int, int> newIndices;
        indices.reserve(column_size * column_size);
        // diagonal state variables
        for (int j = 0; j < species_start; j++)
        {
            // non diagonal state variable derivative elements
            for (int k = 0; k < column_size; k++)
            {
                indices.push_back(k + j * column_size); // traverse row
                indices.push_back(column_size * (k) + j); // traverse column
            }
        }
        // get reaction derivative indices
        reaction_derv_mgr.getDerivativeIndices(&indices);
        // sort
        std::sort(indices.begin(), indices.end(), std::less<int>());
        // only keep unique elements
        vector_int::iterator unique_itr = std::unique(indices.begin(), indices.end());
        // resize based on unique elements
        indices.resize( std::distance(indices.begin(), unique_itr));
        // reduce indices capacity to it's size
        indices.shrink_to_fit();
        // convert indices to network based
        int currCol = currStart;
        int counter = 0;
        for (auto j = indices.begin(); j != indices.end(); j++)
        {
            // current index
            int cidx = *j;
            // associated preconditioner column
            int nextCol = cidx / column_size + currStart;
            int nextRow = cidx % column_size + currStart;
            int flatIdx = nextRow + nextCol * totalColLen;
            // check that index is not already in the state map
            if (m_state_indices.find(flatIdx) == m_state_indices.end())
            {
                // associated preconditioner row
                m_inner.push_back(nextRow);
                // global state index
                m_state_indices.insert(std::pair<int, int>(flatIdx, m_nnz));
                // remapped index
                newIndices.insert(std::pair<int, int> (cidx, counter));
                // find outer index locations
                if (nextCol != currCol)
                {
                    for (int k = 0; k < std::abs(nextCol-currCol); k++)
                    {
                        m_outer.push_back(m_nnz);
                    }
                    currCol = nextCol;
                }
                m_nnz++;
                counter++;
            }
        }
        // add final size to m_outer for current reactor
        m_outer.push_back(m_nnz);
        // set next size of reactor indices array
        m_sizes.push_back(m_nnz);
        // remap derivative indices to sparse structure
        reaction_derv_mgr.remapDerivativeIndices(&newIndices);
    }
    // shrink appropriate vectors
    m_inner.shrink_to_fit();
    m_outer.shrink_to_fit();
    // create space for working vector
    m_values.resize(m_inner.size());
    // reserve space for preconditioner
    m_precon_matrix.resize(m_dimensions[0], m_dimensions[1]);
    m_precon_matrix.reserve(m_nnz);
    // creating sparse identity matrix
    m_identity.resize(m_dimensions[0], m_dimensions[1]);
    m_identity.setIdentity();
    m_identity.makeCompressed();
    // update initialized status
    m_init = true;
}

void AdaptivePreconditioner::acceptReactor(MoleReactor& reactor, double t, double* N, double* Ndot, double* params)
{
    reactor.reactorPreconditionerSetup(*this, t, N, Ndot, params);
}

void AdaptivePreconditioner::acceptReactor(IdealGasMoleReactor& reactor, double t, double* N, double* Ndot, double* params)
{
     reactor.reactorPreconditionerSetup(*this, t, N, Ndot, params);
}

void AdaptivePreconditioner::acceptReactor(IdealGasConstPressureMoleReactor& reactor, double t, double* N, double* Ndot, double* params)
{
     reactor.reactorPreconditionerSetup(*this, t, N, Ndot, params);
}

void AdaptivePreconditioner::setup()
{
    // make into preconditioner as P = (I - gamma * J_bar)
    transformJacobianToPreconditioner();
    // compressing sparse matrix structure
    m_precon_matrix.makeCompressed();
    // analyze and factorize
    m_solver.analyzePattern(m_precon_matrix);
    m_solver.factorize(m_precon_matrix);
    // check for errors
    preconditionerErrorCheck();
}

void AdaptivePreconditioner::transformJacobianToPreconditioner()
{
    Eigen::Map<Eigen::SparseMatrix<double>> jacobian(m_dimensions[0], m_dimensions[1], m_nnz, m_outer.data(), m_inner.data(), m_values.data());
    m_precon_matrix = m_identity - m_gamma * jacobian;
    for (int k=0; k<m_precon_matrix.outerSize(); ++k)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(m_precon_matrix, k); it; ++it)
        {
            if (std::abs(it.value()) < m_threshold && it.row() != it.col())
            {
                m_precon_matrix.coeffRef(it.row(), it.col()) = 0;
            }
        }
    }
}

void AdaptivePreconditioner::preconditionerErrorCheck()
{
    if (m_solver.info() != Eigen::Success)
    {
        throw CanteraError("AdaptivePreconditioner::solve", m_solver.lastErrorMessage());
    }
}

void AdaptivePreconditioner::solve(const size_t state_len, double *rhs_vector, double* output)
{
    // creating vectors in the form of Ax=b
    Eigen::Map<Eigen::VectorXd> bVector(rhs_vector, state_len);
    Eigen::Map<Eigen::VectorXd> xVector(output, state_len);
    // solve for xVector
    xVector = m_solver.solve(bVector);
    preconditionerErrorCheck();
}

}
