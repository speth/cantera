//! @file IdealGasConstPressureMoleReactor.cpp A constant pressure zero-dimensional reactor

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

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

void MoleReactor::acceptPreconditioner(PreconditionerBase *preconditioner, double t, double* N, double* Ndot, double* params)
{
    preconditioner->reactorLevelSetup(this, t, N, Ndot, params);
}

}
