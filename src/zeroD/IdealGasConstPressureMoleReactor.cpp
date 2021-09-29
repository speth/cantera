//! @file IdealGasConstPressureMoleReactor.cpp A constant pressure zero-dimensional reactor

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/IdealGasConstPressureMoleReactor.h"
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

void IdealGasConstPressureMoleReactor::setThermoMgr(ThermoPhase& thermo)
{
    //! @TODO: Add a method to ThermoPhase that indicates whether a given
    //! subclass is compatible with this reactor model
    if (thermo.type() != "IdealGas") {
        throw CanteraError("IdealGasConstPressureMoleReactor::setThermoMgr",
                           "Incompatible phase type provided");
    }
    Reactor::setThermoMgr(thermo);
}

void IdealGasConstPressureMoleReactor::getState(double* N)
{
    if (m_thermo == 0) {
        throw CanteraError("IdealGasConstPressureMoleReactor::getState",
                           "Error: reactor is empty.");
    }
    m_thermo->restoreState(m_state);

    // get mass for calculations
    m_mass = m_thermo->density() * m_vol;
    // set the second component to the temperature
    N[0] = m_thermo->temperature();
    // Use inv. molecular weights and mass to get moles from mass fractionsconst std::vector<double>& imw = m_thermo->inverseMolecularWeights();
    const double* Y = m_thermo->massFractions();
    const std::vector<double>& imw = m_thermo->inverseMolecularWeights();
    double *Ns = N + m_sidx;
    for (size_t i = 0; i < m_nsp; i++)
    {
        Ns[i] = m_mass * imw[i] * Y[i];
    }
    // set the remaining components to the surface species
    // coverages on the walls
    getSurfaceInitialConditions(N + m_nsp + m_sidx);
}

void IdealGasConstPressureMoleReactor::getSurfaceInitialConditions(double* N)
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

void IdealGasConstPressureMoleReactor::initialize(double t0)
{
    ConstPressureReactor::initialize(t0);
    m_nv -= 1; // moles gives the state one fewer variables
    m_hk.resize(m_nsp, 0.0);
}

void IdealGasConstPressureMoleReactor::updateState(double* N)
{
    // The components of N are: [0] the temperature,
    // [1...K+1) are the moles of each species, and [K+1...] are the
    // coverages of surface species on each wall.
    m_thermo->setMolesNoNorm(N + m_sidx);
    m_thermo->setState_TP(N[0], m_pressure);
    // get mass
    std::vector<double> mass(m_nv-m_sidx);
    const std::vector<double>& mw = m_thermo->molecularWeights();
    copy(N + m_sidx, N + m_nv, mass.begin());
    transform(mass.begin(), mass.end(), mw.begin(),
              mass.begin(), multiplies<double>());
    m_mass = accumulate(mass.begin(), mass.end(), 0.0);
    m_vol = m_mass / m_thermo->density();
    updateSurfaceState(N + m_nsp + m_sidx);
    updateConnected(false);
}

void IdealGasConstPressureMoleReactor::updateSurfaceState(double* N)
{
    size_t loc = 0;
    for (auto& S : m_surfaces) {
        auto currPhase = S->thermo();
        currPhase->setMolesNoNorm(N+loc);
        loc += currPhase->nSpecies();
    }
}

void IdealGasConstPressureMoleReactor::evalEqs(double time, double* N,
                                   double* Ndot, double* params)
{
    double mcpdTdt = 0.0; // m * c_p * dT/dt
    double* dNdt = Ndot + m_sidx; //kmol per s

    applySensitivity(params);
    evalWalls(time);

    m_thermo->restoreState(m_state);
    evalSurfaces(time, dNdt + m_nsp);

    m_thermo->getPartialMolarEnthalpies(&m_hk[0]);
    const std::vector<double>& imw = m_thermo->inverseMolecularWeights();

    if (m_chem) {
        m_kin->getNetProductionRates(&m_wdot[0]); // "omega dot"
    }

    // external heat transfer
    mcpdTdt -= m_Q;

    for (size_t n = 0; n < m_nsp; n++) {
        // heat release from gas phase and surface reactions
        mcpdTdt -= m_wdot[n] * m_hk[n] * m_vol;
        mcpdTdt -= m_sdot[n] * m_hk[n];
        // production in gas phase and from surfaces
        dNdt[n] = (m_wdot[n] * m_vol + m_sdot[n]);
    }

    // add terms for outlets
    for (auto outlet : m_outlet) {
        for (size_t n = 0; n < m_nsp; n++) {
            // flow of species into system and dilution by other species
            dNdt[n] -= outlet->outletSpeciesMolarFlowRate(n);
        }
    }

    // add terms for inlets
    for (auto inlet : m_inlet) {
        double mdot = inlet->massFlowRate();
        mcpdTdt += inlet->enthalpy_mass() * mdot;
        for (size_t n = 0; n < m_nsp; n++) {
            double mdot_spec = inlet->outletSpeciesMassFlowRate(n);
            // flow of species into system and dilution by other species
            dNdt[n] += inlet->outletSpeciesMolarFlowRate(n);
            mcpdTdt -= m_hk[n] * imw[n] * mdot_spec;
        }
    }

    if (m_energy) {
        Ndot[0] = mcpdTdt / (m_mass * m_thermo->cp_mass());
    } else {
        Ndot[0] = 0.0;
    }

    resetSensitivity(params);
}

double IdealGasConstPressureMoleReactor::evalSurfaces(double t, double* Ndot)
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

size_t IdealGasConstPressureMoleReactor::componentIndex(const string& nm) const
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

std::string IdealGasConstPressureMoleReactor::componentName(size_t k) {
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
    throw CanteraError("IdealGasConstPressureMoleReactor::componentName",
                       "Index is out of bounds.");
}

void IdealGasConstPressureMoleReactor::acceptPreconditioner(PreconditionerBase *preconditioner, double t, double* N, double* Ndot, double* params)
{
    preconditioner->reactorLevelSetup(this, t, N, Ndot, params);
}

}
