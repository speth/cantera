//! @file MoleReactor.h

// This file is part of Cantera. See License.txt in the top-level
// directory or at https://cantera.org/license.txt for license and
// copyright information.

#ifndef CT_MOLEREACTOR_H
#define CT_MOLEREACTOR_H

#include "Reactor.h"
#include "cantera/numerics/AdaptivePreconditioner.h"
#include "cantera/kinetics/ReactionDerivativeManager.h"

namespace Cantera
{

/*!
 * MoleReactor is a class base class similar which use a state of moles.
 * The reactor may have an arbitrary number of inlets and outlets, each
 * of which may be connected to a "flow device" such as a mass flow
 * controller, a pressure regulator, etc. Additional reactors may be
 * connected to the other end of the flow device, allowing construction
 * of arbitrary reactor networks.
 */
class MoleReactor : public Reactor
{
public:
    MoleReactor(){};

    //! Deprecated function for returning type as a string
    virtual std::string typeStr() const {
        warn_deprecated("MoleReactor::typeStr",
                        "To be removed after Cantera 2.6. Use type() instead.");
        return "MoleReactor";
    }

    //! Use this function to return a string of reactor type
    virtual std::string type() const {
        return "MoleReactor";
    }

    //! Use this function to run reactor specific initialization
    virtual void initialize(double t0 = 0.0);

    //! Return the index in the solution vector for this MoleReactor of
    //! the component named *nm*. Possible values for *nm* are "mass",
    //! "volume", "int_energy", the name of a homogeneous phase species,
    //! or the name of a surface species.
    virtual size_t componentIndex(const std::string& nm) const;

    //! Return the name of the solution component with index *i*.
    //! @see componentIndex()
    virtual std::string componentName(size_t k);

    //! This is a function to setup an associated portion of the
    //! preconditioner for this reactor and is part of the visitor
    //! design pattern.
    //! @param preconditioner the preconditioner being used by cvodes
    //! @param t current time of the simulation
    //! @param N state vector in moles
    //! @param Ndot derivative vector in moles per second
    //! @param params sensitivity parameters
    virtual void preconditionerSetup(PreconditionerBase& preconditioner, double t, double* N, double* Ndot, double* params){
        preconditioner.acceptReactor(*this, t, N, Ndot, params);
    }

    //! This function is the next level of preconditioner setup used in
    //! the visitor design pattern. This is necessary for determining
    //! specific types of both the reactor and preconditioner object
    //! @param preconditioner the preconditioner being used by cvodes
    //! @param t current time of the simulation
    //! @param N state vector in moles
    //! @param Ndot derivative vector in moles per second
    //! @param params sensitivity parameters
    virtual void reactorPreconditionerSetup(AdaptivePreconditioner& preconditioner, double t, double* N, double* Ndot, double* params);

    //! Use this function to precondition the supplied preconditioner
    //! with species variable related derivatives. It can be overloaded
    //! for multiple derivative types.
    //! @param preconditioner the preconditioner being used by cvodes
    //! @param t current time of the simulation
    //! @param N state vector in moles
    //! @param Ndot derivative vector in moles per second
    //! @param params sensitivity parameters
    virtual void SpeciesSpeciesDerivatives(AdaptivePreconditioner& preconditioner, double* N);

    //! Use this function to precondition the supplied preconditioner
    //! with state variable related derivatives. It can be overloaded
    //! for multiple derivative types.
    //! @param preconditioner the preconditioner being used by cvodes
    //! @param t current time of the simulation
    //! @param N state vector in moles
    //! @param Ndot derivative vector in moles per second
    //! @param params sensitivity parameters
    virtual void StateDerivatives(AdaptivePreconditioner& preconditioner, double t, double* N, double* Ndot, double* params)
    {
        throw NotImplementedError("MoleReactor::StateDerivatives");
    }

    //! Return the associated ReactionDerivativeManager object
    virtual ReactionDerivativeManager& reaction_derivative_manager() {return m_reaction_derivative_mgr;};

    //! Return species start in state
    virtual size_t species_start() {
        return m_sidx;};

    //! Return number of nonzero jacobian elements
    virtual size_t nonzero_jacobian_elements();

protected:
    //! TODO: Fix for mole based reactors
    //! Evaluate terms related to surface reactions. Calculates #m_sdot
    //! and rate of change in surface species coverages.
    //! @param t          the current time
    //! @param[out] ydot  array of d(coverage)/dt for surface species
    //! @returns          Net mass flux from surfaces
    virtual double evalSurfaces(double t, double* ydot);

    //! Update the state of SurfPhase objects attached to this
    //! MoleReactor
    virtual void updateSurfaceState(double* y);

    //! Get initial conditions for SurfPhase objects attached to this
    //! MoleReactor
    virtual void getSurfaceInitialConditions(double* y);

    //! const value for the species start index
    const int m_sidx = 1;

    //! Reaction derivative manager
    ReactionDerivativeManager m_reaction_derivative_mgr;
};

}

#endif
