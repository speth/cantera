/**
 *  @file PreconditionerBase.h
 *   Declarations for the class PreconditionerBase which is a virtual base class for preconditioning systems.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef PRECONDITIONERBASE_H
#define PRECONDITIONERBASE_H

#include "cantera/base/ctexceptions.h"


namespace Cantera
{

//! Forward Declarations
class ReactorNet;
class Reactor;
class IdealGasConstPressureMoleReactor;
class IdealGasMoleReactor;

//! Flag to indicate preconditioner is not set
const int PRECONDITIONER_NOT_SET = 0;

class PreconditionerBase
{
    protected:
        //! @param dimensions an size_t pointer to store dimensions
        std::vector<size_t> m_dimensions;

        //! @param m_gamma gamma value used in M = I - gamma*J
        double m_gamma = 1.0;

    public:
        PreconditionerBase(/* args */){}
        ~PreconditionerBase(){}
        //! This function returns zero for preconditioner not set
        virtual size_t getPreconditionerType(){return PRECONDITIONER_NOT_SET;};

        //! Function to solve a linear system Ax=b where A is the
        //! preconditioner contained in this matrix
        //! @param ReactorNet object for integrating
        //! @param[in] t time.
        //! @param[in] N solution vector, length neq()
        //! @param[out] Ndot rate of change of solution vector, length neq()
        //! @param[in] rhs right hand side vector used in linear system
        //! @param[out] output guess vector used by GMRES
        virtual void solve(ReactorNet* network, double *rhs_vector, double* output){
            throw NotImplementedError("PreconditionerBase::solve");
        };

        //! This function performs the setup of the preconditioner for
        //! the specified reactor type and should be overloaded for each
        //! different reactor time
        //! @param ReactorNet object for integrating
        //! @param[in] t time.
        //! @param[in] N solution vector, length neq()
        //! @param[out] Ndot rate of change of solution vector, length neq()
        virtual void setup(ReactorNet* network, double t, double* N, double* Ndot, double gamma){
            throw NotImplementedError("PreconditionerBase::setup");
        };

        //! This function is called during setup for any processes that need to be completed prior to setup functions
        //! @param network A pointer to the reactor net object associated with the integration
        virtual void initialize(ReactorNet* network){
            throw NotImplementedError("PreconditionerBase::initialize");
        };

        //! Function used to set a specific element of the matrix structure
        //! @param row size_t specifying the row location
        //! @param col size_t specifying the column location
        //! @param element double value to be inserted into matrix structure
        virtual void setElement(size_t row, size_t col, double element){
            throw NotImplementedError("PreconditionerBase::setElement");
        };

        //! Function used to get a specific element of the matrix structure
        //! @param row size_t specifying the row location
        //! @param col size_t specifying the column location
        virtual double getElement(size_t row, size_t col){
            throw NotImplementedError("PreconditionerBase::getElement");
        };

        //! Function used to set gamma
        //! @param gamma used in M = I - gamma*J
        virtual void setGamma(double gamma){
            m_gamma = gamma;
        };

        //! Function used to get gamma
        virtual double getGamma(){
            return m_gamma;
        };

        //! Function used to complete individual reactor setups
        //! @param reactor A Reactor pointer
        virtual void reactorLevelSetup(Reactor* reactor, double t, double* N, double* Ndot, double* params){
        throw NotImplementedError("PreconditionerBase-Reactor::reactorLevelSetup");
        };

        //! Function used to complete individual reactor setups
        //! @param reactor A IdealGasConstPressureMoleReactor pointer
        //! @param reactorStart an size_t providing the index location in which the state of the given reactor starts
        virtual void reactorLevelSetup(IdealGasConstPressureMoleReactor* reactor, double t, double* N, double* Ndot, double* params){
        throw NotImplementedError("PreconditionerBase-IdealGasConstPressureMoleReactor::reactorLevelSetup");
        };

        //! Function used to complete individual reactor setups
        //! @param reactor A IdealGasConstPressureMoleReactor pointer
        //! @param reactorStart an size_t providing the index location in which the state of the given reactor starts
        virtual void reactorLevelSetup(IdealGasMoleReactor* reactor, double t, double* N, double* Ndot, double* params){
        throw NotImplementedError("PreconditionerBase-IdealGasMoleReactor::reactorLevelSetup");
        };

        //! Function used to set dimensions of the preconditioner
        //! @param dims A pointer to an array of the dimensions
        void setDimensions(std::vector<size_t> *dims)
        {
            this->m_dimensions.clear();
            for (auto it = dims->begin(); it != dims->end(); ++it)
            {
                this->m_dimensions.push_back(*it);
            }
        }

        //! Function to return pointer to dimensions
        std::vector<size_t>* getDimensions()
        {
            return &(this->m_dimensions);
        }
};

}
#endif
