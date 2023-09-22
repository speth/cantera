//! @file LmrRate.cpp
// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/LmrRate.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/kinetics/Kinetics.h" //! @todo remove after Cantera 3.0 (only used for deprecation)

namespace Cantera{
void LmrData::update(double T){
    throw CanteraError("LmrData::update",
        "Missing state information: 'LmrData' requires pressure.");
}
bool LmrData::update(const ThermoPhase& phase, const Kinetics& kin){
    double T = phase.temperature();
    double P = phase.pressure();
    double X = phase.stateMFNumber(); //Need to figure out the correct way to write this
    if (P != pressure || T != temperature) {
        update(T, P);
        return true;
    }
    return false;
}
void LmrData::perturbPressure(double deltaP){
    if (m_pressure_buf > 0.) {
        throw CanteraError("LmrData::perturbPressure",
            "Cannot apply another perturbation as state is already perturbed.");
    }
    m_pressure_buf = pressure;
    update(temperature, pressure * (1. + deltaP));
}
void LmrData::restore(){
    ReactionData::restore();
    // only restore if there is a valid buffered value
    if (m_pressure_buf < 0.) {
        return;
    }
    update(temperature, m_pressure_buf);
    m_pressure_buf = -1.;
}

// Methods of class LmrRate
LmrRate::LmrRate(const AnyMap& node, const UnitStack& rate_units){
    setParameters(node, rate_units);
}

void LmrRate::setParameters(const AnyMap& node, const UnitStack& rate_units){
    ReactionRate::setParameters(node, rate_units);
    std::map<string,ArrheniusRate> all_eig0; // key = colliderName, value = eig0
    std::map<string, std::multimap<double,ArrheniusRate>> all_multiRates; //key=colliderName, value=(1xN list of rates)
    if (node.hasKey("collider-list")) {
        auto& colliders = node["collider-list"].asVector<AnyMap>();
        for (const auto& collider : colliders) { //iterate through the list (vector) of collider species
            if (collider.hasKey("name") && collider.hasKey("low-P-rate-constant") && collider.hasKey("rate-constants")) {
                std::string key = collider["name"].as<std::string>();
                all_eig0.insert({key, ArrheniusRate(AnyValue(collider["low-P-rate-constant"]), node.units(), rate_units)});

                std::multimap<double, ArrheniusRate> multi_rates;
                auto& rates = collider["rate-constants"].asVector<AnyMap>();
                for (const auto& rate : rates){
                    multi_rates.insert({rate.convert("P","Pa"),ArrheniusRate(AnyValue(rate), node.units(), rate_units)});
                }
                all_multiRates.insert({key,multi_rates});
            }
        }
    }
    // setRates(all_multiRates);
    ////CODE BELOW IS PASTED FROM SETRATES, BUT HAS NOT YET BEEN MODIFIED
    // size_t j = 0;
    // rates_.clear();
    // pressures_.clear();
    // m_valid = !rates.empty();
    // rates_.reserve(rates.size());
    // // Insert intermediate pressures
    // for (const auto& [pressure, rate] : rates) {
    //     double logp = std::log(pressure);
    //     if (pressures_.empty() || pressures_.rbegin()->first != logp) {
    //         // starting a new group
    //         pressures_[logp] = {j, j+1};
    //     } else {
    //         // another rate expression at the same pressure
    //         pressures_[logp].second = j+1;
    //     }
    //     j++;
    //     rates_.push_back(rate);
    // }
    // if (!m_valid) {
    //     // ensure that reaction rate can be evaluated (but returns NaN)
    //     rates_.reserve(1);
    //     pressures_[std::log(OneBar)] = {0, 0};
    //     rates_.push_back(ArrheniusRate());
    // }
    // // Duplicate the first and last groups to handle P < P_0 and P > P_N
    // pressures_.insert({-1000.0, pressures_.begin()->second});
    // pressures_.insert({1000.0, pressures_.rbegin()->second});
}

void LmrRate::getParameters(AnyMap& rateNode, const Units& rate_units) const{ //COMBINES GETPARAMETERS AND GETRATES - CAN MAKE SIMPLIFICATIONS TO THIS CODE
    vector<AnyMap> rateList;
    if (!valid()) {
        // object not fully set up
        return;
    }
    std::multimap<double, ArrheniusRate> rateMap;
    // initial preincrement to skip rate for P --> 0
    for (auto iter = ++pressures_.begin();
            iter->first < 1000; // skip rates for (P --> infinity)
            ++iter) {
        for (size_t i = iter->second.first; i < iter->second.second; i++) {
            rateMap.insert({std::exp(iter->first), rates_[i]});
        }
    }
    for (const auto& [pressure, rate] : rateMap) {
        AnyMap rateNode_;
        rateNode_["P"].setQuantity(pressure, "Pa");
        rate.getRateParameters(rateNode_);
        rateList.push_back(std::move(rateNode_));
    }
    rateNode["rate-constants"] = std::move(rateList);
}

void LmrRate::validate(const string& equation, const Kinetics& kin){ //[HAVEN'T YET MODIFIED THIS]
//     if (!valid()) {
//         throw InputFileError("LmrRate::validate", m_input,
//             "Rate object for reaction '{}' is not configured.", equation);
//     }
//     fmt::memory_buffer err_reactions;
//     double T[] = {300.0, 500.0, 1000.0, 2000.0, 5000.0, 10000.0};
//     LmrData data;
//     for (auto iter = ++pressures_.begin(); iter->first < 1000; iter++) {
//         data.update(T[0], exp(iter->first)); // iter->first contains log(p)
//         updateFromStruct(data);
//         for (size_t i=0; i < 6; i++) {
//             double k = 0;
//             for (size_t p = ilow1_; p < ilow2_; p++) {
//                 k += rates_.at(p).evalRate(log(T[i]), 1.0 / T[i]);
//             }
//             if (!(k > 0)) {
//                 fmt_append(err_reactions,
//                     "at P = {:.5g}, T = {:.1f}\n", std::exp(iter->first), T[i]);
//             }
//         }
// v
//     }
//     if (err_reactions.size()) {
//         throw InputFileError("LmrRate::validate", m_input,
//             "\nInvalid rate coefficient for reaction '{}'\n{}",
//             equation, to_string(err_reactions));
//     }
}

double LmrRate::evalFromStruct(const LmrData& shared_data) {
    X=shared_data.X;
    eig0=shared_data.eig0; 
    rates=shared_data.rates;
    double eig_0_mix = 0.0; //eig_0_mix = np.sum(np.array(X) * np.array(eig_0))
    for (int i = 0; i < eig0.size(); i++) {
        eig_0_mix += X[i]*eig0[i];
    }   
    double k_LMR = 0.0;
    for (int i=0; i<eig0.size(); i++){ //each iteration performs calculations for a different species
        Xtilde=eig0[i]*X[i]/eig0_mix;
        // Check to see if pressure is within range
        if (shared_data.logP != logP_) {
            logP_=log10(shared_data.P*eig0_mix/eig0[i]); //need to use natl log instead to align with code later on
            if (logP_ > logP1_ && logP_ < logP2_) {
                return;
            }
            auto iter = pressures_.upper_bound(logP_); //locate the pressure of interest
            AssertThrowMsg(iter != pressures_.end(), "LmrRate::evalFromStruct",
                            "Pressure out of range: {}", logP_);
            AssertThrowMsg(iter != pressures_.begin(), "LmrRate::evalFromStruct",
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
            double log_k1, log_k2;
            if (ilow1_ == ilow2_) {
                log_k1 = rates_[ilow1_].evalLog(shared_data.logT, shared_data.recipT);
            } else {
                double k = 1e-300; // non-zero to make log(k) finite
                for (size_t i = ilow1_; i < ilow2_; i++) {
                    k += rates_[i].evalRate(shared_data.logT, shared_data.recipT);
                }
                log_k1 = std::log(k);
            }
            if (ihigh1_ == ihigh2_) {
                log_k2 = rates_[ihigh1_].evalLog(shared_data.logT, shared_data.recipT);
            } else {
                double k = 1e-300; // non-zero to make log(k) finite
                for (size_t i = ihigh1_; i < ihigh2_; i++) {
                    k += rates_[i].evalRate(shared_data.logT, shared_data.recipT);
                }
                log_k2 = std::log(k);
            }
            k_LMR += std::exp(log_k1 + (log_k2-log_k1)/(logP2_-logP1_) * (logP_ - logP1_ + std::log(eig0_mix) - std::log(eig0[i])));
        }
    }
    return k_LMR;
}
}

