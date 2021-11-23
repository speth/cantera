//! @file MoleReactor.cpp A zero-dimensional reactor with a moles as the
//! state

// This file is part of Cantera. See License.txt in the top-level
// directory or at https://cantera.org/license.txt for license and
// copyright information.

#include "cantera/zeroD/MoleReactor.h"
#include "cantera/zeroD/FlowDevice.h"
#include "cantera/zeroD/ReactorNet.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/zeroD/ReactorSurface.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/base/utilities.h"

using namespace std;

namespace Cantera
{

void MoleReactor::getSurfaceInitialConditions(double* N)
{
    size_t loc = 0;
    for (auto& S : m_surfaces) {
        auto currPhase = S->thermo();
        currPhase->getConcentrations(N + loc);
        double area = S->area();
        for (size_t i = loc; i < loc + currPhase->nSpecies(); i++)
        {
            N[i] *= area;
        }
        loc += currPhase->nSpecies();
    }
}

void MoleReactor::initialize(double t0)
{
    Reactor::initialize(t0);
    m_nv -= 2; // moles gives the state one fewer variables
    m_reaction_derivative_mgr.initialize(*this);
}

void MoleReactor::updateSurfaceState(double* N)
{
    size_t loc = 0;
    for (auto& S : m_surfaces) {
        auto currPhase = S->thermo();
        currPhase->setMolesNoNorm(N+loc);
        loc += currPhase->nSpecies();
    }
}

double MoleReactor::evalSurfaces(double t, double* Ndot)
{
    const vector_fp& mw = m_thermo->molecularWeights();
    fill(m_sdot.begin(), m_sdot.end(), 0.0);
    size_t loc = 0; // offset into ydot
    double mdot_surf = 0.0; // net mass flux from surface

    for (auto S : m_surfaces) {
        Kinetics* kin = S->kinetics();
        SurfPhase* surf = S->thermo();
        double wallarea = S->area();
        size_t nk = surf->nSpecies();
        surf->setTemperature(m_state[0]);
        S->syncCoverages();
        kin->getNetProductionRates(&m_work[0]);
        size_t ns = kin->surfacePhaseIndex();
        size_t surfloc = kin->kineticsSpeciesIndex(0,ns);
        for (size_t k = 0; k < nk; k++) {
            Ndot[loc + k] = m_work[surfloc+k] * wallarea;
        }
        loc += nk;
        size_t bulkloc = kin->kineticsSpeciesIndex(m_thermo->speciesName(0));
        for (size_t k = 0; k < m_nsp; k++) {
            m_sdot[k] += m_work[bulkloc + k] * wallarea;
            mdot_surf += m_sdot[k] * mw[k];
        }
    }
    return mdot_surf;
}

size_t MoleReactor::componentIndex(const string& nm) const
{
    size_t k = speciesIndex(nm);
    if (k != npos) {
        return k + m_sidx;
    } else if (nm == "temperature") {
        return 0;
    } else {
        return npos;
    }
}

std::string MoleReactor::componentName(size_t k) {
    if (k == 0) {
        return "temperature";
    } else if (k >= 1 && k < neq()) {
        k -= 1;
        if (k < m_thermo->nSpecies()) {
            return m_thermo->speciesName(k);
        } else {
            k -= m_thermo->nSpecies();
        }
        for (auto& S : m_surfaces) {
            ThermoPhase* th = S->thermo();
            if (k < th->nSpecies()) {
                return th->speciesName(k);
            } else {
                k -= th->nSpecies();
            }
        }
    }
    throw CanteraError("MoleReactor::componentName",
                       "Index is out of bounds.");
}

void MoleReactor::reactorPreconditionerSetup(AdaptivePreconditioner& preconditioner, double t, double* N, double* Ndot, double* params)
{
    // strictly positive composition
    vector_fp NCopy(m_nv);
    preconditioner.getStrictlyPositiveComposition(m_nv, N, NCopy.data());
    updateState(NCopy.data());
    // species derivatives
    SpeciesSpeciesDerivatives(preconditioner, NCopy.data());
    // state derivatives
    if (m_energy)
    {
        StateDerivatives(preconditioner, t, NCopy.data(), Ndot, params);
    }
}

void MoleReactor::SpeciesSpeciesDerivatives(AdaptivePreconditioner& preconditioner, double* N)
{
    // getting rate constant data
    size_t numberOfReactions = m_kin->nReactions();
    std::vector<double> kForward (numberOfReactions, 0.0);
    std::vector<double> kBackward (numberOfReactions, 0.0);
    m_kin->getFwdRateConstants(kForward.data());
    m_kin->getRevRateConstants(kBackward.data());
    m_kin->thirdbodyConcMultiply(kForward.data());
    m_kin->thirdbodyConcMultiply(kBackward.data());
    // getting concentrations for derivatives
    std::vector<double> concs(m_nv, 0.0);
    m_thermo->getConcentrations(concs.data() + m_sidx);
    // scale by volume for molar derivative
    scale(kForward.begin(), kForward.end(), kForward.begin(), 1/m_vol);
    scale(kBackward.begin(), kBackward.end(), kBackward.begin(), 1/m_vol);
    // calculating derivatives with reaction manager
    double* derivatives = preconditioner.m_values.data() +  preconditioner.m_sizes[preconditioner.m_ctr];
    m_reaction_derivative_mgr.getDerivatives(concs.data(), derivatives, kForward.data(), kBackward.data());
}

size_t MoleReactor::nonzero_jacobian_elements()
{
    size_t nonzeros = 2 * m_sidx * m_nv - m_sidx;
    return nonzeros + m_reaction_derivative_mgr.getNumNonzeros();
}

}
