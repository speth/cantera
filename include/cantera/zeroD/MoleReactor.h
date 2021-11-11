//! @file MoleReactor.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_MOLEREACTOR_H
#define CT_MOLEREACTOR_H

#include "Reactor.h"

namespace Cantera
{

class MoleReactor : public Reactor
{
public:
    MoleReactor(){};

    virtual std::string typeStr() const {
        warn_deprecated("MoleReactor::typeStr",
                        "To be removed after Cantera 2.6. Use type() instead.");
        return "MoleReactor";
    }

    virtual std::string type() const {
        return "MoleReactor";
    }

    virtual void initialize(doublereal t0 = 0.0);

    //! Return the index in the solution vector for this MoleReactor of the
    //! component named *nm*. Possible values for *nm* are "mass", "volume",
    //! "int_energy", the name of a homogeneous phase species, or the name of a
    //! surface species.
    virtual size_t componentIndex(const std::string& nm) const;

    //! Return the name of the solution component with index *i*.
    //! @see componentIndex()
    virtual std::string componentName(size_t k);

    //!This is a function to accept a preconditioner and perform an action based on MoleReactor type.
    //!@param preconditioner a preconditioner base subclass for preconditioning the system
    //!@param reactorStart start of the MoleReactor within the network
    //!@param t, @param c, @param cdot, @param params double pointers used in integration
    virtual void acceptPreconditioner(PreconditionerBase *preconditioner, double t, double* c, double* cdot, double* params);

    //! const value for the species start index
    const int m_sidx = 1;

protected:
    //! Evaluate terms related to surface reactions. Calculates #m_sdot and rate
    //! of change in surface species coverages.
    //! @param t          the current time
    //! @param[out] ydot  array of d(coverage)/dt for surface species
    //! @returns          Net mass flux from surfaces
    virtual double evalSurfaces(double t, double* ydot);

    //! Update the state of SurfPhase objects attached to this MoleReactor
    virtual void updateSurfaceState(double* y);

    //! Get initial conditions for SurfPhase objects attached to this MoleReactor
    virtual void getSurfaceInitialConditions(double* y);
};
}

#endif
