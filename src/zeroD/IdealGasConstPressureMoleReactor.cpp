//! @file IdealGasConstPressureMoleReactor.cpp A constant pressure
//! zero-dimensional reactor with moles as the state

// This file is part of Cantera. See License.txt in the top-level
// directory or at https://cantera.org/license.txt for license and
// copyright information.

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
    // Use inverse molecular weights
    const double* Y = m_thermo->massFractions();
    const vector_fp& imw = m_thermo->inverseMolecularWeights();
    double *Ns = N + m_sidx;
    for (size_t i = 0; i < m_nsp; i++)
    {
        Ns[i] = m_mass * imw[i] * Y[i];
    }
    // set the remaining components to the surface species coverages on
    // the walls
    getSurfaceInitialConditions(N + m_nsp + m_sidx);
}

void IdealGasConstPressureMoleReactor::initialize(double t0)
{
    MoleReactor::initialize(t0);
    m_hk.resize(m_nsp, 0.0);
}

void IdealGasConstPressureMoleReactor::updateState(double* N)
{
    // the components of N are: [0] the temperature, [1...K+1) are the
    // moles of each species, and [K+1...] are the coverages of surface
    // species on each wall.
    m_thermo->setMolesNoNorm(N + m_sidx);
    m_thermo->setState_TP(N[0], m_pressure);
    // get mass
    vector_fp mass(m_nv-m_sidx);
    const vector_fp& mw = m_thermo->molecularWeights();
    copy(N + m_sidx, N + m_nv, mass.begin());
    transform(mass.begin(), mass.end(), mw.begin(),
              mass.begin(), multiplies<double>());
    m_mass = accumulate(mass.begin(), mass.end(), 0.0);
    m_vol = m_mass / m_thermo->density();
    updateSurfaceState(N + m_nsp + m_sidx);
    updateConnected(false);
}

void IdealGasConstPressureMoleReactor::evalEqs(double time, double* N,
                                   double* Ndot, double* params)
{
    double mcpdTdt = 0.0; // m * c_p * dT/dt
    double* dNdt = Ndot + m_sidx; // kmol per s

    applySensitivity(params);
    evalWalls(time);

    m_thermo->restoreState(m_state);
    evalSurfaces(time, dNdt + m_nsp);

    m_thermo->getPartialMolarEnthalpies(&m_hk[0]);
    const vector_fp& imw = m_thermo->inverseMolecularWeights();

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

void IdealGasConstPressureMoleReactor::StateDerivatives(AdaptivePreconditioner& preconditioner, double t, double* N, double* Ndot, double* params)
{
    // getting perturbed state for finite difference
    double deltaTemp = N[0] * preconditioner.m_perturb;
    // finite difference temperature derivatives
    vector_fp NNext (m_nv);
    vector_fp NdotNext (m_nv);
    vector_fp NCurrent (m_nv);
    vector_fp NdotCurrent (m_nv);
    // copy N to current and next
    copy(N, N + m_nv, NCurrent.begin());
    copy(N, N + m_nv, NNext.begin());
    // perturb temperature
    NNext[0] += deltaTemp;
    // getting perturbed state
    updateState(NNext.data());
    evalEqs(t, NNext.data(), NdotNext.data(), params);
    // reset and get original state
    updateState(NCurrent.data());
    evalEqs(t, NCurrent.data(), NdotCurrent.data(), params);
    // d T_dot/dT
    preconditioner(0, 0) = (NdotNext[0] - NdotCurrent[0]) / deltaTemp;
    // d omega_dot_j/dT
    for (size_t j = m_sidx; j < m_nv; j++)
    {
        preconditioner(j, 0) = (NdotNext[j] - NdotCurrent[j]) / deltaTemp;
    }
    // d T_dot/dcj
    vector_fp specificHeat (m_nsp);
    vector_fp netProductionRates (m_nsp);
    vector_fp enthalpy (m_nsp);
    vector_fp concentrations (m_nsp);
    // getting species concentrations
    m_thermo->getConcentrations(concentrations.data());
    m_thermo->getPartialMolarCp(specificHeat.data());
    m_thermo->getPartialMolarEnthalpies(enthalpy.data());
    m_kin->getNetProductionRates(netProductionRates.data());
    // getting perturbed changes w.r.t temperature
    double CkCpkSum = 0;
    double hkwkSum = 0;
    double inverseVolume = 1/volume();
    double cp_mole = m_thermo->cp_mole();
    for (size_t i = 0; i < m_nsp; i++)
    {
        hkwkSum += enthalpy[i] * netProductionRates[i];
        CkCpkSum += concentrations[i] * specificHeat[i];
    }
    for (size_t j = 0; j < m_nsp; j++) // spans columns
    {
        double hkdwkdnjSum = 0;
        for (size_t k = 0; k < m_nsp; k++) // spans rows
        {
            hkdwkdnjSum += enthalpy[k] * preconditioner(k+m_sidx, j + m_sidx);
        }
        // set appropriate column of preconditioner
        preconditioner(0, j + m_sidx) = -(hkdwkdnjSum * CkCpkSum - (specificHeat[j] - cp_mole) * inverseVolume * hkwkSum) / (CkCpkSum * CkCpkSum);
    }
}

}
