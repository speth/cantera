//! @file ReactionData.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/ReactionData.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/base/ctexceptions.h"

namespace Cantera
{

void ReactionData::update(double T, double extra)
{
    throw NotImplementedError("ReactionData::update",
        "ReactionData type does not use extra argument.");
}

bool ReactionData::perturbT(double deltaT)
{
    // only perturb if there is no buffered value
    if (m_temperature_buf > 0.) {
        return false;
    }
    m_temperature_buf = temperature;
    ReactionData::update(temperature * (1. + deltaT));
    return true;
}

bool ReactionData::restore()
{
    // only restore if there is a valid buffered value
    if (m_temperature_buf < 0.) {
        return false;
    }
    ReactionData::update(m_temperature_buf);
    m_temperature_buf = -1.;
    return true;
}

bool ArrheniusData::update(const ThermoPhase& bulk, const Kinetics& kin)
{
    double T = bulk.temperature();
    if (T == temperature) {
        return false;
    }
    update(T);
    return true;
}

bool TwoTempPlasmaData::update(const ThermoPhase& bulk, const Kinetics& kin)
{
    double T = bulk.temperature();
    double Te = bulk.electronTemperature();
    bool changed = false;
    if (T != temperature) {
        ReactionData::update(T);
        changed = true;
    }
    if (Te != electronTemp) {
        updateTe(Te);
        changed = true;
    }
    return changed;
}

void TwoTempPlasmaData::update(double T)
{
    throw CanteraError("TwoTempPlasmaData::update",
        "Missing state information: 'TwoTempPlasmaData' requires electron temperature.");
}

void TwoTempPlasmaData::update(double T, double Te)
{
    ReactionData::update(T);
    updateTe(Te);
}

void TwoTempPlasmaData::updateTe(double Te)
{
    electronTemp = Te;
    logTe = std::log(Te);
    recipTe = 1./Te;
}

BlowersMaselData::BlowersMaselData()
    : ready(false)
    , density(NAN)
    , m_state_mf_number(-1)
{
    dH.resize(1, NAN);
}

void BlowersMaselData::update(double T)
{
    throw CanteraError("BlowersMaselData::update",
        "Missing state information: 'BlowersMaselData' requires enthalpy change.");
}

void BlowersMaselData::update(double T, double deltaH)
{
    ReactionData::update(T);
    dH[0] = deltaH;
}

bool BlowersMaselData::update(const ThermoPhase& bulk, const Kinetics& kin)
{
    double rho = bulk.density();
    int mf = bulk.stateMFNumber();
    double T = bulk.temperature();
    bool changed = false;
    if (T != temperature) {
        ReactionData::update(T);
        changed = true;
    }
    if (changed || rho != density || mf != m_state_mf_number) {
        density = rho;
        m_state_mf_number = mf;
        bulk.getPartialMolarEnthalpies(m_grt.data());
        kin.getReactionDelta(m_grt.data(), dH.data());
        changed = true;
    }
    return changed;
}

FalloffData::FalloffData()
    : ready(false)
    , molar_density(NAN)
    , m_state_mf_number(-1)
    , m_perturbed(false)
{
    conc_3b.resize(1, NAN);
    m_conc_3b_buf.resize(1, NAN);
}

void FalloffData::update(double T)
{
    throw CanteraError("FalloffData::update",
        "Missing state information: 'FalloffData' requires third-body concentration.");
}

void FalloffData::update(double T, double M)
{
    ReactionData::update(T);
    conc_3b[0] = M;
}

bool FalloffData::update(const ThermoPhase& bulk, const Kinetics& kin)
{
    double rho_m = bulk.molarDensity();
    int mf = bulk.stateMFNumber();
    double T = bulk.temperature();
    bool changed = false;
    if (T != temperature) {
        ReactionData::update(T);
        changed = true;
    }
    if (rho_m != molar_density || mf != m_state_mf_number) {
        molar_density = rho_m;
        m_state_mf_number = mf;
        conc_3b = kin.thirdBodyConcentrations();
        changed = true;
    }
    return changed;
}

bool FalloffData::perturbM(double deltaM)
{
    // only perturb if there is no buffered value
    if (m_perturbed) {
        return false;
    }
    m_conc_3b_buf = conc_3b;
    for (auto& c3b : conc_3b) {
        c3b *= 1. + deltaM;
    }
    m_perturbed = true;
    return true;
}

bool FalloffData::restore()
{
    bool ret = ReactionData::restore();
    // only restore if there is a valid buffered value
    if (!m_perturbed) {
        return ret;
    }
    conc_3b = m_conc_3b_buf;
    m_perturbed = false;
    return true;
}

void PlogData::update(double T)
{
    throw CanteraError("PlogData::update",
        "Missing state information: 'PlogData' requires pressure.");
}

bool PlogData::update(const ThermoPhase& bulk, const Kinetics& kin)
{
    double T = bulk.temperature();
    double P = bulk.pressure();
    if (P != pressure || T != temperature) {
        update(T, P);
        return true;
    }
    return false;
}

bool PlogData::perturbP(double deltaP)
{
    // only perturb if there is no buffered value
    if (m_pressure_buf > 0.) {
        return false;
    }
    m_pressure_buf = pressure;
    update(temperature, pressure * (1. + deltaP));
    return true;
}

bool PlogData::restore()
{
    bool ret = ReactionData::restore();
    // only restore if there is a valid buffered value
    if (m_pressure_buf < 0.) {
        return ret;
    }
    update(temperature, m_pressure_buf);
    m_pressure_buf = -1.;
    return true;
}

void ChebyshevData::update(double T)
{
    throw CanteraError("ChebyshevData::update",
        "Missing state information: 'ChebyshevData' requires pressure.");
}

bool ChebyshevData::update(const ThermoPhase& bulk, const Kinetics& kin)
{
    double T = bulk.temperature();
    double P = bulk.pressure();
    if (P != pressure || T != temperature) {
        update(T, P);
        return true;
    }
    return false;
}

bool ChebyshevData::perturbP(double deltaP)
{
    // only perturb if there is no buffered value
    if (m_pressure_buf > 0.) {
        return false;
    }
    m_pressure_buf = pressure;
    update(temperature, pressure * (1. + deltaP));
    return true;
}

bool ChebyshevData::restore()
{
    bool ret = ReactionData::restore();
    // only restore if there is a valid buffered value
    if (m_pressure_buf < 0.) {
        return ret;
    }
    update(temperature, m_pressure_buf);
    m_pressure_buf = -1.;
    return true;
}

}
