//! @file LmrRate.h
// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_LMRRATE_H
#define CT_LMRRATE_H
#include "cantera/kinetics/Arrhenius.h"

namespace Cantera{
struct LmrData : public ReactionData{
    // LmrData() = default;
    LmrData();
    virtual bool update(const ThermoPhase& phase, const Kinetics& kin) override;
    using ReactionData::update;
    void perturbPressure(double deltaP);
    virtual void restore() override;

    virtual void resize(size_t nSpecies, size_t nReactions, size_t nPhases) override {
        moleFractions.resize(nReactions, NAN);
        ready = true;
    }

    virtual void invalidateCache() override {
        ReactionData::invalidateCache();
        pressure = NAN;
    }
    double pressure = NAN; //!< pressure
    double logP = 0.0; //!< logarithm of pressure
    bool ready = false; //!< boolean indicating whether vectors are accessible
    vector<double> moleFractions;
    int mfNumber; 
protected:
    double m_pressure_buf = -1.0; //!< buffered pressure
};

class LmrRate final : public ReactionRate
{
public:
    LmrRate() = default;//! Default constructor.
    LmrRate(const AnyMap& node, const UnitStack& rate_units={});
    unique_ptr<MultiRateBase> newMultiRate() const {
        return make_unique<MultiRate<LmrRate, LmrData>>();
    } 
    const string type() const { return "LMR_R"; } //! Identifier of reaction rate type
    
    void setParameters(const AnyMap& node, const UnitStack& rate_units);
    void validate(const string& equation, const Kinetics& kin);
    double speciesPlogRate(const LmrData& shared_data);
    double evalFromStruct(const LmrData& shared_data);
    void getParameters(AnyMap& rateNode, const Units& rate_units) const;
    

    vector<string> allSpecies_; //list of all yaml species (not just those for which LMRR data exists)
    map<string, map<double, pair<size_t, size_t>>> pressures_;
    map<string, vector<ArrheniusRate>> rates_;
    map<string,ArrheniusRate> eig0_;
    map<double, pair<size_t, size_t>> pressures_s_;
    vector<ArrheniusRate> rates_s_;
    double logPeff_;
    double eig0_mix_ = 0.0;
    double log_eig0_mix_ = 0.0;
    double k_LMR_ = 0.0;

protected:
    double logP_ = -1000;
    double logP1_ = 1000;
    double logP2_ = -1000;
    size_t ilow1_, ilow2_, ihigh1_, ihigh2_;
    double rDeltaP_ = -1.0; //!< reciprocal of (logP2 - logP1)
};
}
#endif