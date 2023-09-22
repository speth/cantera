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
        logP = std::log(P);
    }
    virtual bool update(const ThermoPhase& phase, const Kinetics& kin) override; //This is the important one
    using ReactionData::update;
    void perturbPressure(double deltaP);
    virtual void restore() override;
    virtual void invalidateCache() override {
        ReactionData::invalidateCache();
        pressure = NAN;
    }
    double pressure = NAN; //!< pressure
    double logP = 0.0; //!< logarithm of pressure
protected:
    double m_pressure_buf = -1.0; //!< buffered pressure
    //somedatatype m_molefraction = someinitializedvalue
};

class LmrRate{
public:
    //! Default constructor.
    LmrRate() = default;
    //! Constructor from Arrhenius rate expressions at a set of pressures
    explicit LmrRate(const std::multimap<double, ArrheniusRate>& rates);
    LmrRate(const AnyMap& node, const UnitStack& rate_units={});
    unique_ptr<MultiRateBase> newMultiRate() const {
        return make_unique<MultiRate<LmrRate, LmrData>>();
    } 
    //! Identifier of reaction rate type
    const string type() const { return "LMR_R"; }
    //! Perform object setup based on AnyMap node information
    /*!
     *  @param node  AnyMap containing rate information
     *  @param rate_units  Unit definitions specific to rate information
     */
    void setParameters(const AnyMap& node, const UnitStack& rate_units);
    void getParameters(AnyMap& rateNode, const Units& rate_units) const;
    void getParameters(AnyMap& rateNode) const {
        return getParameters(rateNode, Units(0));
    }
    //! Update information specific to reaction
    /*!
     *  @param shared_data  data shared by all reactions of a given type
     */
    void updateFromStruct(const LmrData& shared_data) {
        if (shared_data.logP != logP_) {
            logP_ = shared_data.logP;
            if (logP_ > logP1_ && logP_ < logP2_) {
                return;
            }
            auto iter = pressures_.upper_bound(logP_);
            AssertThrowMsg(iter != pressures_.end(), "LmrRate::updateFromStruct",
                           "Pressure out of range: {}", logP_);
            AssertThrowMsg(iter != pressures_.begin(), "LmrRate::updateFromStruct",
                           "Pressure out of range: {}", logP_);
            // upper interpolation pressure
            logP2_ = iter->first;
            ihigh1_ = iter->second.first;
            ihigh2_ = iter->second.second;
            // lower interpolation pressure
            logP1_ = (--iter)->first;
            ilow1_ = iter->second.first;
            ilow2_ = iter->second.second;
            rDeltaP_ = 1.0 / (logP2_ - logP1_);
        }
    }
    //! Evaluate reaction rate
    /*!
     *  @param shared_data  data shared by all reactions of a given type
     */
    double evalFromStruct(const LmrData& shared_data);

    //! Set up LMR object
    void setRates(const std::multimap<string, std::multimap<double,ArrheniusRate>>& all_rates)
    // void setRates(const std::multimap<double, ArrheniusRate>& rates);
    //! Check to make sure that the rate expression is finite over a range of
    //! temperatures at each interpolation pressure. This is potentially an
    //! issue when one of the Arrhenius expressions at a particular pressure
    //! has a negative pre-exponential factor.
    void validate(const string& equation, const Kinetics& kin);
    //! Return the pressures and Arrhenius expressions which comprise this
    //! reaction.
    std::multimap<double, ArrheniusRate> getRates() const;

protected:
    //! log(p) to (index range) in the rates_ vector
    map<double, pair<size_t, size_t>> pressures_;
    // Rate expressions which are referenced by the indices stored in pressures_
    vector<ArrheniusRate> rates_;
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