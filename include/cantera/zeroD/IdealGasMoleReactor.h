//! @file IdealGasMoleReactor.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_IDEALGASMOLE_REACTOR_H
#define CT_IDEALGASMOLE_REACTOR_H

#include "cantera/zeroD/MoleReactor.h"

namespace Cantera
{

/*!
 * IdealGasMoleReactor is a class for ideal gas
 * constant-volume reactors which use a state of moles. The reactor
 * may have an arbitrary number of inlets and outlets, each of which may
 * be connected to a "flow device" such as a mass flow controller, a
 * pressure regulator, etc. Additional reactors may be connected to the
 * other end of the flow device, allowing construction of arbitrary
 * reactor networks.
 */
class IdealGasMoleReactor : public MoleReactor
{
public:
    IdealGasMoleReactor() {}

    //! Deprecated function for returning type as a string
    virtual std::string typeStr() const {
        warn_deprecated("IdealGasMoleReactor::typeStr",
                        "To be removed after Cantera 2.6. Use type() instead.");
        return "IdealGasMoleReactor";
    }

    //! Use this function to return a string of reactor type
    virtual std::string type() const {
        return "IdealGasMoleReactor";
    }

    //! Use this function to set the thermo manager of this reactor
    virtual void setThermoMgr(ThermoPhase& thermo);

    //! Use this function to get the state in moles
    virtual void getState(double* N);

    //! Use this function to initialize the reactor
    virtual void initialize(double t0 = 0.0);

    //! Right hand side function used to integrate by CVODES
    //! @param t current time of the simulation
    //! @param N state vector in moles
    //! @param Ndot derivative vector in moles per second
    //! @param params sensitivity parameters
    virtual void evalEqs(double t, double* N, double* Ndot, double* params);

    //! Use to update state vector N
    virtual void updateState(double* N);

    //! Use this function to precondition the supplied preconditioner
    //! with state variable related derivatives. It can be overloaded
    //! for multiple derivative types.
    //! @param preconditioner the preconditioner being used by cvodes
    //! @param t current time of the simulation
    //! @param N state vector in moles
    //! @param Ndot derivative vector in moles per second
    //! @param params sensitivity parameters
    virtual void StateDerivatives(AdaptivePreconditioner& preconditioner, double t, double* N, double* Ndot, double* params);

protected:
    vector_fp m_uk; //!< Species molar enthalpies
};

}

#endif
