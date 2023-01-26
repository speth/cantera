//! @file ReactionData.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/ReactionData.h"

namespace Cantera
{

void ReactionData::update(double T)
{
    temperature = T;
    logT = std::log(T);
    recipT = 1./T;
}

void ReactionData::update(double T, double extra)
{
    throw NotImplementedError("ReactionData::update",
        "ReactionData type does not use extra scalar argument.");
}

void ReactionData::update(double T, const vector_fp& extra)
{
    throw NotImplementedError("ReactionData::update",
        "ReactionData type does not use extra vector argument.");
}

void ReactionData::perturbTemperature(double deltaT)
{
    if (m_temperature_buf > 0.) {
        throw CanteraError("ReactionData::perturbTemperature",
            "Cannot apply another perturbation as state is already perturbed.");
    }
    m_temperature_buf = temperature;
    ReactionData::update(temperature * (1. + deltaT));
}

void ReactionData::restore()
{
    // only restore if there is a valid buffered value
    if (m_temperature_buf < 0.) {
        return;
    }
    ReactionData::update(m_temperature_buf);
    m_temperature_buf = -1.;
}

void ReactionData::resize(size_t nSpecies, size_t nReactions, size_t nPhases)
{
}

void ReactionData::invalidateCache()
{
    temperature = NAN;
}

}
