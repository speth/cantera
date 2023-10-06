//! @file LmrRate.h
// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_LMRRATE_H
#define CT_LMRRATE_H
#include "cantera/kinetics/Arrhenius.h"

namespace Cantera{
struct LmrData : public ReactionData{
    LmrData() = default;
    virtual void update(double T) override;
    virtual void update(double T, double P) override {
        ReactionData::update(T);
        pressure = P;
    }
    virtual bool update(const ThermoPhase& phase);
    using ReactionData::update;
    void perturbPressure(double deltaP);
    virtual void restore() override;
    virtual void invalidateCache() override {
        ReactionData::invalidateCache();
        pressure = NAN;
    }
    // double moleFraction = NAN;
    double pressure = NAN; //!< pressure
    // double logP = 0.0; //!< logarithm of pressure
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
    // void getParameters(AnyMap& rateNode, const Units& rate_units) const;
    // void getParameters(AnyMap& rateNode) const {
    //     return getParameters(rateNode, Units(0));
    // }
    double computeSpeciesRate(const LmrData& shared_data);
    double evalFromStruct(const LmrData& shared_data);
    void validate(const string& equation);

    map<string, map<double, pair<size_t, size_t>>> pressures_;
    map<string, vector<ArrheniusRate>> rates_;
    map<string,ArrheniusRate> eig0_;
    double Xs_;
    string s_;
    double eig0_mix_ = 0.0;
    double k_LMR_ = 0.0;

protected:
    //species_ = FXNTOEXTRACTSPECIESLISTFROMSOLUTIONOBJECT
    double P_ = 1e-300;
    double P1_ = 1e20;
    double P2_ = 1e-300;
    size_t ilow1_, ilow2_, ihigh1_, ihigh2_;
};
}
#endif