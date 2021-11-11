//! @file AdaptivePreconditioner.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/numerics/AdaptivePreconditioner.h"
#include "cantera/base/utilities.h"
#include "cantera/base/global.h"
#include <iostream>

namespace Cantera
{

    AdaptivePreconditioner::AdaptivePreconditioner(const AdaptivePreconditioner &externalPrecon)
    {
        m_threshold = externalPrecon.m_threshold;
        m_values = externalPrecon.m_values;
        m_inner = externalPrecon.m_inner;
        m_outer = externalPrecon.m_outer;
        m_nnz = externalPrecon.m_nnz;
        m_precon_matrix = externalPrecon.m_precon_matrix;
        m_dimensions.clear();
        m_dimensions.push_back(externalPrecon.m_dimensions.at(0));
        m_dimensions.push_back(externalPrecon.m_dimensions.at(1));
    }

    bool AdaptivePreconditioner::operator== (const AdaptivePreconditioner &externalPrecon)
    {
        // only compare preconditioner matrix, will fail if it has not been formed
        return m_precon_matrix.isApprox(externalPrecon.m_precon_matrix);
    }

    void AdaptivePreconditioner::operator= (const AdaptivePreconditioner &externalPrecon)
    {
        m_threshold = externalPrecon.m_threshold;
        m_values = externalPrecon.m_values;
        m_inner = externalPrecon.m_inner;
        m_outer = externalPrecon.m_outer;
        m_nnz = externalPrecon.m_nnz;
        m_precon_matrix = externalPrecon.m_precon_matrix;
        m_atol = externalPrecon.m_atol;
        m_dimensions.clear();
        m_dimensions.push_back(externalPrecon.m_dimensions.at(0));
        m_dimensions.push_back(externalPrecon.m_dimensions.at(1));
        return;
    }

    double& AdaptivePreconditioner::operator() (size_t row, size_t col)
    {
        size_t index = __gidx(row, col);
        if (m_state_indices.find(index) != m_state_indices.end())
        {
            size_t idx = m_state_indices.at(index);
            return m_values[idx];
        }
        else
        {
            m_zero[0] = 0; // reset to zero just incase it was assigned
            return m_zero[0];
        }
    }

    void AdaptivePreconditioner::initialize(ReactorNet* network)
    {
        //! Don't use legacy rate constants
        use_legacy_rate_constants(false);
        // Reset arrays in case of re-initialization
        m_dimensions.clear();
        m_outer.clear();
        m_inner.clear();
        m_thermo_indices.clear();
        m_values.clear();
        m_sizes.clear();
        m_reaction_derv_mgrs.clear();
        m_state_indices.clear();
        m_nnz = 0;
        // Set dimensions of preconditioner from network
        size_t totalColLen = network->m_nv;
        m_dimensions.push_back(totalColLen);
        m_dimensions.push_back(totalColLen);
        if (m_dimensions[0] != m_dimensions[1])
        {
            throw CanteraError("initialize", "specified matrix dimensions are not square");
        }
        // Set zero parameter
        m_zero[0] = 0;
        // Set absolute tolerance from network
        m_atol = network->m_atols;
        // Reserve maximum space for vectors making up SparseMatrix
        m_inner.reserve(totalColLen * totalColLen);
        m_outer.reserve(totalColLen * totalColLen);
        m_outer.push_back(0); // first outer starts at zero
        m_thermo_indices.reserve(3 * network->m_reactors.size() * totalColLen);
        // Allocate spaces for sizes and set first start to zero
        m_sizes.reserve(network->m_reactors.size() + 1);
        m_sizes.push_back(0);
        // Set up sparse patterns
        m_reaction_derv_mgrs.resize(network->m_reactors.size());
        for (size_t i = 0; i < m_reaction_derv_mgrs.size(); i++)
        {
            auto currReactor = network->m_reactors[i];
            size_t currStart = network->m_start[i];
            m_reaction_derv_mgrs[i].initialize(currReactor);
            // Checking manager
            int column_size = currReactor->neq();
            int species_start = column_size-currReactor->getThermoMgr()->nSpecies();
            // Create temporary indices vector and reserve space for it
            std::vector<int> indices;
            std::unordered_map<int, int> newIndices;
            indices.reserve(column_size * column_size);
            // Diagonal state variables
            for (int j = 0; j < species_start; j++)
            {
                // Non diagonal state variable derivative elements
                for (int k = 0; k < column_size; k++)
                {
                    indices.push_back(k + j * column_size); // traverse row
                    indices.push_back(column_size * (k) + j); // traverse column
                }
            }
            // Get reaction derivative indices
            m_reaction_derv_mgrs[i].getDerivativeIndices(&indices);
            // Sort
            std::sort(indices.begin(), indices.end(), std::less<int>());
            // Only keep unique elements
            std::unique(indices.begin(), indices.end());
            // Reduce indices capacity to it's size
            indices.shrink_to_fit();
            // Convert indices to network based
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
                    //remapped index
                    newIndices.insert(std::pair<int, int> (cidx, counter));
                    // Find outer index locations
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
            // Add final size to m_outer for current reactor
            m_outer.push_back(m_nnz);
            // Set next size of reactor indices array
            m_sizes.push_back(m_nnz);
            // remap derivative indices to sparse structure
            m_reaction_derv_mgrs[i].remapDerivativeIndices(&newIndices);
        }
        // Shrink appropriate vectors
        m_inner.shrink_to_fit();
        m_outer.shrink_to_fit();
        // Create space for working vector
        m_values.resize(m_inner.size());
        // Reserve space for preconditioner
        m_precon_matrix.resize(m_dimensions[0], m_dimensions[1]);
        m_precon_matrix.reserve(m_nnz);
        // Creating sparse identity matrix
        m_identity.resize(m_dimensions[0], m_dimensions[1]);
        m_identity.setIdentity();
        m_identity.makeCompressed();
    }

    void AdaptivePreconditioner::setup(ReactorNet* network, double t, double* N, double* Ndot, double gamma)
    {
        // Set gamma value for M =I - gamma*J
        m_gamma = gamma;
        // Setting to zero to refill
        reset();
        // Calling
        for (size_t n = 0; n < network->m_reactors.size(); n++)
        {
            m_ri = n; // set reactor counter
            m_current_start = network->m_start[n]; // set reactor start for indexing
            // Reactor start is not added to N and Ndot because the
            // preconditioner is not broken into units like reactors are
            (network->m_reactors[n])->acceptPreconditioner(this, t, N + m_current_start, Ndot + m_current_start, (network->m_sens_params).data());
        }
        // Make into preconditioner as P = (I - gamma * J_bar)
        transformJacobianToPreconditioner();
        // Compressing sparse matrix structure
        m_precon_matrix.makeCompressed();
        // Analyze and factorize
        m_solver.analyzePattern(m_precon_matrix);
        m_solver.factorize(m_precon_matrix);
        // Check for errors
        preconditionerErrorCheck();
    }

    void AdaptivePreconditioner::reactorLevelSetup(IdealGasConstPressureMoleReactor* reactor, double t, double* N, double* Ndot, double* params)
    {
        // strictly positive composition
        std::vector<double> NCopy(reactor->neq());
        getStrictlyPositiveComposition(reactor->neq(), N, NCopy.data());
        reactor->updateState(NCopy.data());
        // Species Derivatives
        SpeciesSpeciesDerivatives(reactor, NCopy.data());
        // Temperature Derivatives
        if (reactor->energyEnabled())
        {
            TemperatureDerivatives(reactor, t, NCopy.data(), Ndot, params);
        }
    }

    void AdaptivePreconditioner::reactorLevelSetup(IdealGasMoleReactor* reactor, double t, double* N, double* Ndot, double* params)
    {
        // strictly positive composition
        std::vector<double> NCopy(reactor->neq());
        getStrictlyPositiveComposition(reactor->neq(), N, NCopy.data());
        reactor->updateState(NCopy.data());
        // Species Derivatives
        SpeciesSpeciesDerivatives(reactor, NCopy.data());
        // Temperature Derivatives
        if (reactor->energyEnabled())
        {
            TemperatureDerivatives(reactor, t, NCopy.data(), Ndot, params);
        }
    }

    void AdaptivePreconditioner::SpeciesSpeciesDerivatives(MoleReactor* reactor, double* N)
    {
        // Getting rate constant data
        auto kinetics = reactor->getKineticsMgr();
        auto thermo = reactor->getThermoMgr();
        size_t numberOfReactions = kinetics->nReactions();
        std::vector<double> kForward (numberOfReactions, 0.0);
        std::vector<double> kBackward (numberOfReactions, 0.0);
        kinetics->getFwdRateConstants(kForward.data());
        kinetics->getRevRateConstants(kBackward.data());
        kinetics->thirdbodyConcMultiply(kForward.data());
        kinetics->thirdbodyConcMultiply(kBackward.data());
        std::vector<double> concs(reactor->neq(), 0.0);
        thermo->getConcentrations(concs.data() + reactor->m_sidx);
        scale(kForward.begin(), kForward.end(), kForward.begin(), 1/reactor->volume());
        scale(kBackward.begin(), kBackward.end(), kBackward.begin(), 1/reactor->volume());
        // Calculating derivatives with reaction manager
        m_reaction_derv_mgrs[m_ri].getDerivatives(concs.data(), m_values.data() + m_sizes[m_ri], kForward.data(), kBackward.data());
    }

    void AdaptivePreconditioner::TemperatureDerivatives(MoleReactor* reactor, double t, double* N, double* Ndot, double* params)
    {
        auto kinetics = reactor->getKineticsMgr();
        auto thermo = reactor->getThermoMgr();
        // Important sizes to the determination of values
        size_t numberOfSpecies = kinetics->nTotalSpecies();
        // Size of state for current reactor
        size_t stateSize = reactor->neq();
        // Starting idx for species in reactor
        size_t speciesStart  = stateSize - numberOfSpecies;
        // Temperature Index in reactor
        size_t tempIndex = reactor->componentIndex("temperature");
        // Getting perturbed state for finite difference
        double deltaTemp = N[tempIndex] * m_perturb;
        /**
         *
         * Temperature Finite Difference
         *
         **/
        // net production rates and enthalpies for each state
        std::vector<double> NNext (stateSize);
        std::vector<double> NdotNext (stateSize);
        std::vector<double> NCurrent (stateSize);
        std::vector<double> NdotCurrent (stateSize);
        // Copy N to current and next
        copy(N, N + stateSize, NCurrent.begin());
        copy(N, N + stateSize, NNext.begin());
        // perturb temperature
        NNext[tempIndex] += deltaTemp;
        // Getting perturbed state
        reactor->updateState(NNext.data());
        reactor->evalEqs(t, NNext.data(), NdotNext.data(), params);
        // Reset and get original state
        reactor->updateState(NCurrent.data());
        reactor->evalEqs(t, NCurrent.data(), NdotCurrent.data(), params);
        // d T_dot/dT
        (*this)[__gidx(tempIndex, tempIndex)] = (NdotNext[tempIndex] - NdotCurrent[tempIndex]) / deltaTemp;
        /**
         *
         * Species Finite Difference
         *
         **/
        // d omega_dot_j/dT
        for (size_t j = speciesStart; j < stateSize; j++)
        {
            (*this)[__gidx(j, tempIndex)] = (NdotNext[j] - NdotCurrent[j]) / deltaTemp;
        }
        // d T_dot/dcj
        std::vector<double> specificHeat (numberOfSpecies);
        std::vector<double> netProductionRates (numberOfSpecies);
        std::vector<double> enthalpy (numberOfSpecies);
        std::vector<double> concentrations (numberOfSpecies);
        // Getting species concentrations
        thermo->getConcentrations(concentrations.data());
        thermo->getPartialMolarCp(specificHeat.data());
        thermo->getPartialMolarEnthalpies(enthalpy.data());
        kinetics->getNetProductionRates(netProductionRates.data());
        // Getting perturbed changes w.r.t temperature
        double CkCpkSum = 0;
        double hkwkSum = 0;
        double inverseVolume = 1/reactor->volume();
        for (size_t i = 0; i < numberOfSpecies; i++)
        {
            hkwkSum += enthalpy[i] * netProductionRates[i];
            CkCpkSum += concentrations[i] * specificHeat[i];
        }
        for (size_t j = 0; j < numberOfSpecies; j++) // Spans columns
        {
            double hkdwkdnjSum = 0;
            for (size_t k = 0; k < numberOfSpecies; k++) // Spans rows
            {
                hkdwkdnjSum += enthalpy[k] * (*this)[__gidx(k+speciesStart, j + speciesStart)];
            }
            // Set appropriate column of preconditioner
            (*this)[__gidx(tempIndex, j + speciesStart)] = (hkdwkdnjSum * CkCpkSum - specificHeat[j] * inverseVolume * hkwkSum) / (CkCpkSum * CkCpkSum);
        }
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

    void AdaptivePreconditioner::solve(ReactorNet* network, double *rhs_vector, double* output)
    {
        // Creating vectors in the form of Ax=b
        Eigen::Map<Eigen::VectorXd> bVector(rhs_vector, network->m_nv);
        Eigen::Map<Eigen::VectorXd> xVector(output, network->m_nv);
        // Solve for xVector
        xVector = m_solver.solve(bVector);
        preconditionerErrorCheck();
    }

}
