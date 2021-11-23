//! @file IdealGasMoleReactor.cpp A constant volume zero-dimensional
//! reactor with moles as the state

// This file is part of Cantera. See License.txt in the top-level
// directory or at https://cantera.org/license.txt for license and
// copyright information.

#include "cantera/zeroD/IdealGasMoleReactor.h"
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

void IdealGasMoleReactor::setThermoMgr(ThermoPhase& thermo)
{
    if (thermo.type() != "IdealGas") {
        throw CanteraError("IdealGasMoleReactor::setThermoMgr",
                           "Incompatible phase type provided");
    }
    Reactor::setThermoMgr(thermo);
}

void IdealGasMoleReactor::getState(double* N)
{
    if (m_thermo == 0) {
        throw CanteraError("IdealGasMoleReactor::getState",
                           "Error: reactor is empty.");
    }
    m_thermo->restoreState(m_state);

    // get mass for calculations
    m_mass = m_thermo->density() * m_vol;
    // set the second component to the temperature
    N[0] = m_thermo->temperature();
    // use inverse molecular weights
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

void IdealGasMoleReactor::initialize(double t0)
{
   MoleReactor::initialize(t0);
    m_uk.resize(m_nsp, 0.0);
}

void IdealGasMoleReactor::updateState(double* N)
{
    // the components of N are: [0] the temperature, [1...K+1) are the
    // moles of each species, and [K+1...] are the coverages of surface
    // species on each wall. get mass
    vector_fp mass(m_nv-m_sidx);
    const vector_fp& mw = m_thermo->molecularWeights();
    copy(N + m_sidx, N + m_nv, mass.begin());
    transform(mass.begin(), mass.end(), mw.begin(),
              mass.begin(), multiplies<double>());
    m_mass = accumulate(mass.begin(), mass.end(), 0.0);
    // set state
    m_thermo->setMolesNoNorm(N + m_sidx);
    m_thermo->setState_TR(N[0], m_mass / m_vol);
    updateSurfaceState(N + m_nsp + m_sidx);
    updateConnected(true);
}

void IdealGasMoleReactor::evalEqs(double time, double* N,
                                   double* Ndot, double* params)
{
    double mcvdTdt = 0.0; // m * c_v * dT/dt
    double* dNdt = Ndot + m_sidx; // kmol per s

    applySensitivity(params);
    evalWalls(time);

    m_thermo->restoreState(m_state);
    evalSurfaces(time, dNdt + m_nsp);

    m_thermo->getPartialMolarIntEnergies(&m_uk[0]);
    const vector_fp& imw = m_thermo->inverseMolecularWeights();

    if (m_chem) {
        m_kin->getNetProductionRates(&m_wdot[0]); // "omega dot"
    }

    // external heat transfer
    mcvdTdt += - m_pressure * m_vdot - m_Q;

    for (size_t n = 0; n < m_nsp; n++) {
        // heat release from gas phase and surface reactions
        mcvdTdt -= m_wdot[n] * m_uk[n] * m_vol;
        mcvdTdt -= m_sdot[n] * m_uk[n];
        // production in gas phase and from surfaces
        dNdt[n] = (m_wdot[n] * m_vol + m_sdot[n]);
    }

    // add terms for outlets
    for (auto outlet : m_outlet) {
        for (size_t n = 0; n < m_nsp; n++) {
            // flow of species into system and dilution by other species
            dNdt[n] -= outlet->outletSpeciesMolarFlowRate(n);
        }
        double mdot = outlet->massFlowRate();
        mcvdTdt -= mdot * m_pressure * m_vol / m_mass; // flow work
    }

    // add terms for inlets
    for (auto inlet : m_inlet) {
        double mdot = inlet->massFlowRate();
        mcvdTdt += inlet->enthalpy_mass() * mdot;
        for (size_t n = 0; n < m_nsp; n++) {
            double mdot_spec = inlet->outletSpeciesMassFlowRate(n);
            // flow of species into system and dilution by other species
            dNdt[n] += inlet->outletSpeciesMolarFlowRate(n);
            mcvdTdt -= m_uk[n] * imw[n] * mdot_spec;
        }
    }

    if (m_energy) {
        Ndot[0] = mcvdTdt / (m_mass * m_thermo->cv_mass());
    } else {
        Ndot[0] = 0.0;
    }

    resetSensitivity(params);
}

void IdealGasMoleReactor::StateDerivatives(AdaptivePreconditioner& preconditioner, double t, double* N, double* Ndot, double* params)
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
    // find derivatives d T_dot/dNj
    vector_fp specificHeat (m_nsp);
    vector_fp netProductionRates (m_nsp);
    vector_fp internal_energy (m_nsp);
    vector_fp concentrations (m_nsp);
    // getting species concentrations
    m_thermo->getConcentrations(concentrations.data());
    m_thermo->getPartialMolarCp(specificHeat.data());
    m_thermo->getPartialMolarIntEnergies(internal_energy.data());
    m_kin->getNetProductionRates(netProductionRates.data());
    // convert Cp to Cv for ideal gas as Cp - Cv = R
    for (size_t i = 0; i < specificHeat.size(); i++)
    {
        specificHeat[i] -= GasConstant;
    }
    // getting perturbed changes w.r.t temperature
    double CkCvkSum = 0;
    double ukwkSum = 0;
    double inverseVolume = 1/volume();
    for (size_t i = 0; i < m_nsp; i++)
    {
        ukwkSum += internal_energy[i] * netProductionRates[i];
        CkCvkSum += concentrations[i] * specificHeat[i];
    }
    for (size_t j = 0; j < m_nsp; j++) // spans columns
    {
        double ukdwkdnjSum = 0;
        for (size_t k = 0; k < m_nsp; k++) // spans rows
        {
            ukdwkdnjSum += internal_energy[k] * preconditioner(k+m_sidx, j + m_sidx);
        }
        // set appropriate column of preconditioner
        preconditioner(0, j + m_sidx) = (ukdwkdnjSum * CkCvkSum - specificHeat[j] * inverseVolume * ukwkSum) / (CkCvkSum * CkCvkSum);
    }
}

}
