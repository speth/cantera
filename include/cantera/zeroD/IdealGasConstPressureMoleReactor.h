//! @file ConstPressureReactor.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_IDEALGASCONSTPRESSMOLE_REACTOR_H
#define CT_IDEALGASCONSTPRESSMOLE_REACTOR_H

#include "ConstPressureReactor.h"

namespace Cantera
{

/**
 * Class ConstPressureReactor is a class for constant-pressure reactors. The
 * reactor may have an arbitrary number of inlets and outlets, each of which may
 * be connected to a "flow device" such as a mass flow controller, a pressure
 * regulator, etc. Additional reactors may be connected to the other end of the
 * flow device, allowing construction of arbitrary reactor networks.
 */
class IdealGasConstPressureMoleReactor : public ConstPressureReactor
{
public:
    IdealGasConstPressureMoleReactor() {}

    virtual std::string typeStr() const {
        warn_deprecated("IdealGasConstPressureMoleReactor::typeStr",
                        "To be removed after Cantera 2.6. Use type() instead.");
        return "IdealGasConstPressureMoleReactor";
    }
    virtual std::string type() const {
        return "IdealGasConstPressureMoleReactor";
    }
    virtual void setThermoMgr(ThermoPhase& thermo);
    virtual void getState(double* N);
    virtual void initialize(double t0 = 0.0);
    virtual void evalEqs(double t, double* N,
                         double* Ndot, double* params);
    //! Use to update state vector N
    virtual void updateState(double* N);

    //! Evaluate terms related to surface reactions. Calculates #m_sdot and rate
    //! of change in surface species coverages.
    //! @param t          the current time
    //! @param[out] Ndot  array of d(coverage)/dt for surface species
    //! @returns          Net mass flux from surfaces
    virtual double evalSurfaces(double t, double* Ndot);

    //! Update the state of SurfPhase objects via moles
    virtual void updateSurfaceState(double* N);

    //! Get initial conditions for SurfPhase objects attached to this reactor
    virtual void getSurfaceInitialConditions(double* N);

    //! Return the index in the solution vector for this reactor of the
    //! component named *nm*. Possible values for *nm* are "mass",
    //! "temperature", the name of a homogeneous phase species, or the name of a
    //! surface species.
    virtual size_t componentIndex(const std::string& nm) const;
    std::string componentName(size_t k);

    //!This is a function to accept a preconditioner and perform an action based on reactor type.
    //!@param preconditioner a preconditioner base subclass for preconditioning the system
    //!@param t, @param c, @param cdot, @param params double pointers used in integration
    virtual void acceptPreconditioner(PreconditionerBase *preconditioner, double t, double* N, double* Ndot, double* params);
    //! const value for the species start index
    const int m_sidx = 1;
protected:
    std::vector<double> m_hk; //!< Species molar enthalpies
};
}

#endif
