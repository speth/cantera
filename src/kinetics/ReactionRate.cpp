//! @file ReactionRate.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/ReactionRate.h"

namespace Cantera
{

void ReactionRate::setParameters(const AnyMap& node, const UnitStack& units)
{
    m_input = node;
}

void ReactionRate::check(const string& equation)
{
}

void ReactionRate::validate(const string& equation, const Kinetics& kin)
{
}

void ReactionRate::getParameters(AnyMap& node) const {
    throw NotImplementedError("ReactionRate::getParameters",
                              "Not implemented by '{}' object.", type());
}

}
