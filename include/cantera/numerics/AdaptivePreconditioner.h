/**
 *  @file AdaptivePreconditioner.h
 *   Declarations for the class AdaptivePreconditioner
 *   which is a child class of PreconditionerBase for preconditioners used by sundials
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef ADAPTIVEPRECONDITIONER_H
#define ADAPTIVEPRECONDITIONER_H

#if CT_USE_SYSTEM_EIGEN
#include <Eigen/Sparse>
#else
#include "cantera/ext/Eigen/Sparse"
#endif

#include "cantera/numerics/PreconditionerBase.h"
#include "cantera/kinetics/ReactionDerivativeManager.h"
#include "cantera/zerodim.h"
#include "cantera/thermo.h"
#include "cantera/kinetics.h"
#include "cantera/kinetics/Reaction.h"

#include "float.h"
#include <unordered_map>

namespace Cantera
{

//! Flag to indicate adaptive preconditioner is set
const int ADAPTIVE_MECHANISM_PRECON_MATRIX = 1;

class AdaptivePreconditioner : public PreconditionerBase
{
    protected:
        //! @param m_values is the container for the values that are mapped to m_jacobian
        std::vector<double> m_values;

        //! @param m_inner is the container for inner (row) indices in CSC format
        std::vector<int> m_inner;

        //! @param m_outer is the container for outer starts (nnz elements per row)
        std::vector<int> m_outer;

        //! @param m_sizes is the container for starts of each reactor in double
        std::vector<size_t> m_sizes;

        //! @param m_reaction_derv_mgrs vector of reaction derivative managers for each kinetic object
        std::vector<ReactionDerivativeManager> m_reaction_derv_mgrs;

        //! @param m_thermo_indices indices of thermo variables temp, mass, vol, etc.
        std::vector<size_t> m_thermo_indices;

        //! @param m_state_indices total index map of system, the key is a flattened index and the output is the corresponding derivative index
        std::unordered_map<int, int> m_state_indices;

        //! @param m_nnz is the number of non zeros
        int m_nnz = 0;

        //! @param m_identity is the container that is the sparse preconditioner
        Eigen::SparseMatrix<double> m_identity;

        //! @param m_precon_matrix is the container that is the sparse preconditioner
        Eigen::SparseMatrix<double> m_precon_matrix;

        //! @param m_solver is the solver used in solving the linear
        //! system
        Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> m_solver;

        //! @param m_ri index of the current reactor
        size_t m_ri = 0;

        //! @param m_threshold a double value to selectively fill the matrix structure based on this threshold
        double m_threshold = DBL_EPSILON; // default

        //! @param m_current_start an index value for the starting index
        //! of the current reactor in the state
        size_t m_current_start;

        //! @param m_atol absolute tolerance of the ODE solver
        double m_atol = 0;

        //! @param m_zero a value of zero;
        double m_zero[1];

        //! @param m_perturb perturbation constant that is multiplied by temperature for perturbation
        double m_perturb = std::sqrt(DBL_EPSILON);

    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW // Required for mis-alignment of EIGEN matrix
        AdaptivePreconditioner(/* args */){};
        ~AdaptivePreconditioner(){};
        AdaptivePreconditioner(const AdaptivePreconditioner &externalPrecon);

        //! This function is called during setup for any processes that need to be completed prior to setup functions used in sundials.
        //! @param network A pointer to the reactor net object associated with the integration
        void initialize(ReactorNet* network);

        //! Use this function to reset arrays within preconditioner object
        inline void reset()
        {
            m_precon_matrix.reserve(m_nnz);
            m_precon_matrix.setZero();
            std::fill(m_values.begin(), m_values.end(), 0);
        };

        //! This function performs the setup of the preconditioner for the specified reactor type and should be overloaded for each different reactor time
        //! @param ReactorNet object for integrating
        //! @param[in] t time.
        //! @param[in] N solution vector, length neq()
        //! @param[out] Ndot rate of change of solution vector, length
        //! neq()
        //! @param gamma in M = I - gamma*J
        void setup(ReactorNet* network, double t, double* N, double* Ndot, double gamma);

        //! Function used to complete individual reactor setups
        //! @param reactor A IdealGasConstPressureMoleReactor pointer
        void reactorLevelSetup(IdealGasConstPressureMoleReactor* reactor, double t, double* N, double* Ndot, double* params);

        //! Function used to complete individual reactor setups
        //! @param reactor A IdealGasMoleReactor pointer
        void reactorLevelSetup(IdealGasMoleReactor* reactor, double t, double* N, double* Ndot, double* params);

        //! This function determines rate law derivatives of species
        //! with respect to other species specifically it determines the
        //! derivatives of the rate laws of all species with respect to
        //! other species in terms of moles.
        //! @param *reactor A pointer to the current reactor being used
        //! for preconditioning
        void SpeciesSpeciesDerivatives(MoleReactor* reactor, double* N);

        //! This function determines derivatives of Species and Temperature with respect to Temperature for jacobian preconditioning with a finite difference.
        //! @param reactor A pointer to the current reactor being precondition
        //! @parama t A double value of the current time
        //! @param N A pointer to the current state passed from CVODES
        //! @param Ndot A pointer to the current state derivatives
        //! passed from CVODES
        //! @param params A double pointer to sensitivty parameters.
        void TemperatureDerivatives(MoleReactor* reactor, double t, double* N, double* Ndot, double* params);

        //! This function checks if there was an error with eigen and
        //! throws it if so.
        void preconditionerErrorCheck();

        //!Use this function to transform Jacobian into preconditioner
        void transformJacobianToPreconditioner();

        //! Function to solve a linear system Ax=b where A is the preconditioner contained in this matrix
        //! @param ReactorNet object for integrating
        //! @param[in] t time.
        //! @param[in] N solution vector, length neq()
        //! @param[out] Ndot rate of change of solution vector, length neq()
        //! @param[in] rhs right hand side vector used in linear system
        //! @param[out] output guess vector used by GMRES
        void solve(ReactorNet* network, double *rhs_vector, double* output);

        //! This function returns the current type of preconditioner as an integer
        size_t getPreconditionerType(){return ADAPTIVE_MECHANISM_PRECON_MATRIX;};

        //! global flat index from local
        //! @param row local row index
        //! @param col local column index
        inline size_t __gidx(size_t row, size_t col)
        {
            return (row + m_current_start) + (col + m_current_start) * m_dimensions[1];
        }

        // Use this function to get the current Jacobian
        Eigen::SparseMatrix<double> getJacobian()
        {
            Eigen::Map<Eigen::SparseMatrix<double>> jacobian(m_dimensions[0], m_dimensions[1], m_nnz, m_outer.data(), m_inner.data(), m_values.data());
            return jacobian;
        };

        //! Function used to get index start of the current reactor variable
        double getReactorStart(){return m_current_start;};

        //! Function used to get a specific element of the preconditioner matrix
        //! @param row size_t specifying the row location
        //! @param col size_t specifying the column location
        double getElement(size_t row, size_t col){return m_precon_matrix.coeffRef(row, col);}

        //! Function used to return compressed version of the matrix structure
        Eigen::SparseMatrix<double>* getMatrix(){return &(m_precon_matrix);};

        //! Use this function to get the threshold value for setting elements
        double getThreshold(){return m_threshold;};

        //! Use this function to return the used absolute tolerance
        double getAbsoluteTolerance(){return m_atol;};

        //! Use this function to get the pertubation constant
        double getPerturbationConst(){return m_perturb;};

        //! Use this function to get the ratio of nonzero preconditioner elements to the maximum number of elements
        double getSparsityPercentage()
        {
            size_t totalElements = m_dimensions[0] * m_dimensions[0];
            size_t p_nnz = (m_precon_matrix.nonZeros() > 0) ? m_precon_matrix.nonZeros() : m_nnz;
            return 1.0 - p_nnz/totalElements;
        };

        //! Use this function to get a strictly positive composition
        void getStrictlyPositiveComposition(size_t vlen, double* in, double* out)
        {
            for (size_t i = 0; i < vlen; i++)
            {
                out[i] = std::max(in[i], m_atol);
            }
        }
        //! Function used to set index start of the current reactor variable
        void setReactorStart(size_t reactorStart){m_current_start = reactorStart;};

        //! Function used to set a specific element of the preconditioner matrix
        //! @param row size_t specifying the row location
        //! @param col size_t specifying the column location
        //! @param element double value to be inserted into matrix structure
        void setElement(size_t row, size_t col, double element)
        {
            if (std::abs(element) >= m_threshold || row == col)
            {
                m_precon_matrix.coeffRef(row,col) = element;
            }
        }

        //! Use this function to set the threshold value to compare elements against
        //! @param threshold double value used in setting by threshold
        void setThreshold(double threshold){m_threshold = threshold;};

        //! Use this function to set the absolute tolerance in the
        //! solver outside of the network initialization
        //! @param atol the specified tolerance
        void setAbsoluteTolerance(double atol){m_atol = atol;};

        //! Use this function to set the perturbation constant
        //! @param perturb the new pertubation constant
        void setPerturbationConst(double perturb){m_perturb = perturb;};

        // Below this point are overloaded operators for certain operations regarding the preconditioner.

        //! Overloading of the == operator to compare values inside preconditioners
        //! @param externalPrecon - == comparison with this object
        bool operator== (const AdaptivePreconditioner &externalPrecon);

        //! Overloading of the = operator to copy one preconditioner to another
        //! @param externalPrecon. Preconditioner becoming this object
        void operator= (const AdaptivePreconditioner &externalPrecon);

        //! Overloading of the [] operator to assign values to the jacobian
        //! this function assumes that the index is in the index map
        //! @param index the flattened index of the point to be accessed
        inline double& operator[] (int index){return m_values[m_state_indices[index]];};

        //! Overloading of the () operator to assign values to the jacobian
        //! this function does not assume that index is index map
        //! @param row row index of jacobian
        //! @param col column index of jacobian
        double& operator() (size_t row, size_t col);

        //! Below this point are functions related to debugging the preconditioner and subclasses.

        //! @param reactor - the contents of this reactor will be printed
        inline void printReactorComponents(Reactor* reactor)
        {
            for (size_t i = 0; i < reactor->neq(); i++)
            {
                std::cout << reactor->componentName(i) << std::endl;
            }
        };

        //! Print preconditioner contents
        inline void printPreconditioner()
        {
            Eigen::IOFormat HeavyFmt(Eigen::FullPrecision, 0, ", ", ";\n", "[", "]", "[", "]");
            std::cout<<Eigen::MatrixXd(m_precon_matrix).format(HeavyFmt)<<std::endl;
        };

        //! Print jacobian contents
        inline void printJacobian()
        {
            Eigen::Map<Eigen::SparseMatrix<double>> jacobian(m_dimensions[0], m_dimensions[1], m_nnz, m_outer.data(), m_inner.data(), m_values.data());
            std::cout<<Eigen::MatrixXd(jacobian)<<std::endl;
        };

        //! Use this function to fill with ones to check the pattern
        inline void ones()
        {
            for (size_t i = 0; i < m_reaction_derv_mgrs.size(); i++)
            {
                m_reaction_derv_mgrs[i].setOnes(m_values.data() + m_sizes[i]);
            }
        }

};

}

#endif
