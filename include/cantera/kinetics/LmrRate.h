//! @file LmrRate.h
// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_LMRRATE_H
#define CT_LMRRATE_H
#include "cantera/kinetics/Arrhenius.h"

namespace Cantera{
struct LmrData : public ReactionData{
    LmrData() = default;
    ReactionData::ThermoPhase::Phase::speciesNames();
    virtual void update(double T) override;
    virtual void update(double T, double P) override {
        ReactionData::update(T);
        pressure = P;
        // moleFraction = X;
        logP = std::log(P);
    }
    virtual bool update(const ThermoPhase& phase, const Kinetics& kin) override; //This is the important one
    using ReactionData::update; //this is where we get logT
    void perturbPressure(double deltaP);
    virtual void restore() override;
    virtual void invalidateCache() override {
        ReactionData::invalidateCache();
        pressure = NAN;
    }
    // double moleFraction = NAN;
    double pressure = NAN; //!< pressure
    double logP = 0.0; //!< logarithm of pressure
protected:
    double m_pressure_buf = -1.0; //!< buffered pressure
    //somedatatype m_molefraction = someinitializedvalue
    // vector<string> m_speciesNames;
};

class LmrRate final : public ReactionRate
{
public:
    LmrRate() = default;//! Default constructor.
    explicit LmrRate(const std::multimap<double, ArrheniusRate>& rates);//! Constructor from Arrhenius rate expressions at a set of pressures
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
    double computeSpeciesRate(const LmrData& shared_data, string s, double eig0_mix);
    double evalFromStruct(const LmrData& shared_data);
    void validate(const string& equation, const Kinetics& kin);

protected:
    map<string, map<double, pair<size_t, size_t>>> pressures_;
    map<string, vector<ArrheniusRate>> rates_;
    map<string,ArrheniusRate> eig0_;

    //species_ = FXNTOEXTRACTSPECIESLISTFROMSOLUTIONOBJECT

    double logP_ = -1000; //!< log(p) at the current state
    double logP1_ = 1000; //!< log(p) at the lower pressure reference
    double logP2_ = -1000; //!< log(p) at the upper pressure reference
    //! Indices to the ranges within rates_ for the lower / upper pressure, such
    //! that rates_[ilow1_] through rates_[ilow2_] (inclusive) are the rates
    //! expressions which are combined to form the rate at the lower reference
    //! pressure.
    size_t ilow1_, ilow2_, ihigh1_, ihigh2_;
    double rDeltaP_ = -1.0; //!< reciprocal of (logP2 - logP1)
};
}
#endif