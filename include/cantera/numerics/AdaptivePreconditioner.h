/**
 *  @file AdaptivePreconditioner.h Declarations for the class
 *   AdaptivePreconditioner which is a child class of PreconditionerBase
 *   for preconditioners used by sundials
 */

// This file is part of Cantera. See License.txt in the top-level
// directory or at https://cantera.org/license.txt for license and
// copyright information.

#ifndef ADAPTIVEPRECONDITIONER_H
#define ADAPTIVEPRECONDITIONER_H

#include "cantera/numerics/eigen_sparse.h"
#include "cantera/numerics/PreconditionerBase.h"
#include "cantera/zeroD/ReactorNet.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/kinetics/ReactionDerivativeManager.h"
#include "cantera/kinetics/Kinetics.h"
#include "float.h"
#include <unordered_map>

namespace Cantera
{

//! Flag to indicate adaptive preconditioner is set
const int ADAPTIVE_MECHANISM_PRECON_MATRIX = 1;

//! AdaptivePreconditioner a preconditioner designed for use with large
//! mechanisms that leverages sparse solvers. It does this by pruning
//! the preconditioner by a threshold value. It also neglects pressure
//! dependence and thirdbody contributions in its formation and has a
//! finite difference approximation for temperature.
class AdaptivePreconditioner : public PreconditionerBase
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW // Required for mis-alignment of EIGEN matrix
    AdaptivePreconditioner(){};
    ~AdaptivePreconditioner(){};
    AdaptivePreconditioner(const AdaptivePreconditioner &externalPrecon){};

    //! Friend classes
    friend class MoleReactor;
    friend class IdealGasConstPressureMoleReactor;
    friend class IdealGasMoleReactor;

    //! This function is called during setup for any processes that need
    //! to be completed prior to setup functions used in sundials.
    //! @param network A pointer to the reactor net object associated
    //! with the integration
    void initialize(ReactorNet& network);

    //! Use this function to reset arrays within preconditioner object
    void reset(){
        m_precon_matrix.reserve(m_nnz);
        m_precon_matrix.setZero();
        std::fill(m_values.begin(), m_values.end(), 0);
    };

    //! This function performs preconditioner specific post-reactor
    //! setup operations such as factorize.
    void setup();

    //! This function is for a visitor design pattern to determine
    //! preconditioner type with MoleReactor
    void acceptReactor(MoleReactor& reactor, double t, double* N, double* Ndot, double* params);

    //! This function is for a visitor design pattern to determine
    //! preconditioner type with IdealGasMoleReactor
    void acceptReactor(IdealGasMoleReactor& reactor, double t, double* N, double* Ndot, double* params);

    //! This function is for a visitor design pattern to determine
    //! preconditioner type with IdealGasConstPressureMoleReactor
    void acceptReactor(IdealGasConstPressureMoleReactor& reactor, double t, double* N, double* Ndot, double* params);

    //! This function checks if there was an error with eigen and throws
    //! it if so.
    void preconditionerErrorCheck();

    //!Use this function to transform Jacobian vector and write into
    //!preconditioner
    void transformJacobianToPreconditioner();

    //! Function to solve a linear system Ax=b where A is the
    //! preconditioner contained in this matrix
    //! @param[in] state_len length of vectors supplied
    //! @param[in] rhs_vector right hand side vector supplied by cvodes
    //! @param[out] output output vector "z" sent back to cvodes
    void solve(const size_t state_len, double *rhs_vector, double* output);

    //! This function returns the preconditioning method as an integer
    size_t getPreconditionerMethod(){return ADAPTIVE_MECHANISM_PRECON_MATRIX;};

    //! This function returns preconditioning type as an integer
    PreconditionerType getPreconditionerType(){return LEFT_PRECONDITION;};

    //! Use this function to get a global flat index from local
    //! @param row local row index
    //! @param col local column index
    size_t globalIndex(size_t row, size_t col){
        return (row + m_starts[m_ctr]) + (col + m_starts[m_ctr]) * m_dimensions[1];
    };

    // Use this function to get the current Jacobian
    Eigen::SparseMatrix<double> getJacobian(){
        Eigen::Map<Eigen::SparseMatrix<double>> jacobian(m_dimensions[0], m_dimensions[1], m_nnz, m_outer.data(), m_inner.data(), m_values.data());
        return jacobian;
    };

    //! Function used to return pointer to preconditioner matrix
    Eigen::SparseMatrix<double>* getMatrix(){return &(m_precon_matrix);};

    //! Use this function to get the threshold value for setting
    //! elements
    double getThreshold(){return m_threshold;};

    //! Use this function to get the pertubation constant
    double getPerturbationConst(){return m_perturb;};

    //! Use this function to get the ratio of nonzero preconditioner
    //! elements to the maximum number of elements
    double getSparsityPercentage(){
        size_t totalElements = m_dimensions[0] * m_dimensions[0];
        size_t p_nnz = (m_precon_matrix.nonZeros() > 0) ? m_precon_matrix.nonZeros() : m_nnz;
        return 1.0 - ((double) p_nnz)/((double) totalElements);
    };

    //! Use this function to get a strictly positive composition
    void getStrictlyPositiveComposition(size_t vlen, double* in, double* out){
        for (size_t i = 0; i < vlen; i++)
        {
            out[i] = std::max(in[i], m_atol);
        }
    };

    //! Use this function to set the threshold value to compare elements
    //! against
    //! @param threshold double value used in setting by threshold
    void setThreshold(double threshold){m_threshold = threshold;};

    //! Use this function to set the perturbation constant used in
    //! finite difference calculations.
    //! @param perturb the new pertubation constant
    void setPerturbationConst(double perturb){m_perturb = perturb;};

    //! Overloading of the == operator to compare values strictly inside
    //! preconditioner matrix
    //! @param externalPrecon == comparison with this object
    bool operator== (const AdaptivePreconditioner &externalPrecon);
    //! Overloading of the = operator to copy one preconditioner to
    //! another
    //! @param externalPrecon the preconditioner becoming this object
    void operator= (const AdaptivePreconditioner &externalPrecon);

    //! Overloading of the [] operator to assign values to the jacobian
    //! this function assumes that the index is in the index map
    //! @param index the flattened index of the point to be accessed
    double& operator[] (int index){return m_values[m_state_indices[index]];};

    //! Overloading of the () operator to assign values to the jacobian
    //! this function does not assume that index is index map
    //! @param row row index of jacobian
    //! @param col column index of jacobian
    double& operator() (size_t row, size_t col);

    //! Print preconditioner contents
    void printPreconditioner(){
        Eigen::IOFormat HeavyFmt(Eigen::FullPrecision, 0, ", ", ";\n", "[", "]", "[", "]");
        std::cout<<Eigen::MatrixXd(m_precon_matrix).format(HeavyFmt)<<std::endl;
    };

    //! Print jacobian contents
    void printJacobian(){
        Eigen::Map<Eigen::SparseMatrix<double>> jacobian(m_dimensions[0], m_dimensions[1], m_nnz, m_outer.data(), m_inner.data(), m_values.data());
        std::cout<<Eigen::MatrixXd(jacobian)<<std::endl;
    };
protected:
    //! Container for the values that are mapped to m_jacobian
    vector_fp m_values;

    //! Container for inner (row) indices in CSC format
    vector_int m_inner;

    //! Container for outer starts (nnz elements per row)
    vector_int m_outer;

    //! Container for starts of each reactor in double m_values
    std::vector<size_t> m_sizes;

    //! Container for starts of each reactor in normal state
    std::vector<size_t> m_starts;

    //! Indices of thermo variables temp, mass, vol, etc.
    std::vector<size_t> m_thermo_indices;

    //! Total index map of system, the key is a flattened index and the
    //! output is the corresponding derivative index
    std::unordered_map<int, int> m_state_indices;

    //! The number of non zeros
    int m_nnz = 0;

    //! Container that is the sparse preconditioner
    Eigen::SparseMatrix<double> m_identity;

    //! Container that is the sparse preconditioner
    Eigen::SparseMatrix<double> m_precon_matrix;

    //! Solver used in solving the linear system
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> m_solver;

    //! Minimum value a non-diagonal element must be to be included in
    //! the preconditioner
    double m_threshold = DBL_EPSILON; // default

    //! Index value for the starting index A value of zero for returns
    //! and assignments to non-existent values
    double m_zero{0};

    //! Perturbation constant that is multiplied by temperature for
    //! perturbation
    double m_perturb = std::sqrt(DBL_EPSILON);
};

}

#endif
