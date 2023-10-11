//! @file LmrRate.cpp
// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/LmrRate.h"
#include "cantera/thermo/ThermoPhase.h"
// @todo remove after Cantera 3.0 (only used for deprecation)
#include "cantera/kinetics/Kinetics.h"

namespace Cantera{
void LmrData::update(double T){
    throw CanteraError("LmrData::update",
        "Missing state information: 'LmrData' requires pressure.");
}
bool LmrData::update(const ThermoPhase& phase, const Kinetics& kin){
    double T = phase.temperature();
    double P = phase.pressure(); //find out what units this is in
    int X = phase.stateMFNumber();
    if (P != pressure || T != temperature || X != mfNumber) {
        update(T,P);
        phase.getMoleFractions(moleFractions.data());
        mfNumber=X;
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
    if (node.hasKey("collider-list")) {
        auto& colliders = node["collider-list"].asVector<AnyMap>();
        for (const auto& collider : colliders) { //iterate through the list (vector) of collider species
            if (collider.hasKey("name") && collider.hasKey("low-P-rate-constant") && collider.hasKey("rate-constants")) {
                string species_i_ = collider["name"].as<std::string>();
                ArrheniusRate eig0_i_ = ArrheniusRate(AnyValue(collider["low-P-rate-constant"]), node.units(), rate_units);       
                map<double, pair<size_t, size_t>> pressures_i_;
                vector<ArrheniusRate> rates_i_; 
                std::multimap<double, ArrheniusRate> multi_rates;
                auto& rates = collider["rate-constants"].asVector<AnyMap>();
                for (const auto& rate : rates){
                    multi_rates.insert({rate.convert("P","Pa"),ArrheniusRate(AnyValue(rate), node.units(), rate_units)});
                }
                rates_i_.reserve(multi_rates.size());
                m_valid = !multi_rates.empty(); //if rates object empty, m_valid==FALSE. if rates is not empty, m_valid==TRUE
                size_t j = 0;
                for (const auto& [pressure, rate] : multi_rates) { 
                    double logp = std::log(pressure);
                    if (pressures_i_.empty() || pressures_i_.rbegin()->first != logp) {
                        pressures_i_[logp] = {j, j+1};
                    } else {
                        pressures_i_[logp].second = j+1;
                    }
                    j++;
                    rates_i_.push_back(rate); 
                }
                if (!m_valid) { //runs if multi_rates is empty
                    // ensure that reaction rate can be evaluated (but returns NaN)
                    rates_i_.reserve(1);
                    pressures_i_[std::log(OneBar)] = {0, 0};
                    rates_i_.push_back(ArrheniusRate());
                }
                // Duplicate the first and last groups to handle P < P_0 and P > P_N
                pressures_i_.insert({-1000.0, pressures_i_.begin()->second});
                pressures_i_.insert({1000.0, pressures_i_.rbegin()->second});

                rates_.insert({species_i_, rates_i_});
                pressures_.insert({species_i_, pressures_i_});
                eig0_.insert({species_i_, eig0_i_});
            }
        }
    }
}

// //adapted from chebyshevrate
// void LmrRate::getParameters(AnyMap& rateNode) const{
//     if (!valid()) { //valid==FALSE makes if statement == TRUE
//         // object not fully set up
//         return;
//     }
//     auto& colliders = rateNode["collider-list"].asVector<AnyMap>();

//     rateNode["collider-list"]["name"].setQuantity("Species_Name"); //incorrect syntax
//     rateNode["temperature-range"].setQuantity({Tmin(), Tmax()}, "K");
//     rateNode["pressure-range"].setQuantity({Pmin(), Pmax()}, "Pa");
//     size_t nT = m_coeffs.nRows();
//     size_t nP = m_coeffs.nColumns();
//     vector<vector<double>> coeffs2d(nT, vector<double>(nP));
//     for (size_t i = 0; i < nT; i++) {
//         for (size_t j = 0; j < nP; j++) {
//             coeffs2d[i][j] = m_coeffs(i, j);
//         }
//     }
//     // Unit conversions must take place later, after the destination unit system
//     // is known. A lambda function is used here to override the default behavior
//     Units rate_units2 = conversionUnits();
//     auto converter = [rate_units2](AnyValue& coeffs, const UnitSystem& units) {
//         if (rate_units2.factor() != 0.0) {
//             coeffs.asVector<vector<double>>()[0][0] += \
//                 std::log10(units.convertFrom(1.0, rate_units2));
//         } else if (units.getDelta(UnitSystem()).size()) {
//             throw CanteraError("ChebyshevRate::getParameters lambda",
//                 "Cannot convert rate constant with unknown dimensions to a "
//                 "non-default unit system");
//         }
//     };
//     AnyValue coeffs;
//     coeffs = std::move(coeffs2d);
//     rateNode["data"].setQuantity(coeffs, converter);
// }

void LmrRate::validate(const string& equation, const Kinetics& kin){
    // STILL NEED TO FIGURE OUT HOW TO GET SPECIES NAME LIST FROM THERMOPHASE OBJECT
    if (!valid()) {
        throw InputFileError("LmrRate::validate", m_input,
            "Rate object for reaction '{}' is not configured.", equation);
    }
    fmt::memory_buffer err_reactions1; //for k-related errors
    fmt::memory_buffer err_reactions2; //for eig0-related errors
    double T[] = {300.0, 500.0, 1000.0, 2000.0, 5000.0, 10000.0};
    LmrData data;
    // Iterate through the outer map (string to inner map)
    for (const auto& outer_pair : pressures_) {
        const std::string& s = outer_pair.first; //s refers only to the species for which LMR data is provided in yaml (e.g. 'H2O', 'M')
        const std::map<double, std::pair<size_t, size_t>>& inner_map = outer_pair.second;
        for (auto iter = ++inner_map.begin(); iter->first < 1000; iter++) { 
            data.update(T[0], exp(iter->first));
            ilow1_ = iter->second.first;
            ilow2_ = iter->second.second;       
            for (size_t i=0; i < 6; i++) {
                double k = 0;
                for (size_t p = ilow1_; p < ilow2_; p++) {
                    k += rates_[s].at(p).evalRate(log(T[i]), 1.0 / T[i]);
                }
                double eig0 = eig0_[s].evalRate(log(T[i]), 1.0/T[i]);

                if (!(k > 0)){ //flags error if k at a given T, P is not > 0
                    fmt_append(err_reactions1,"at P = {:.5g}, T = {:.1f}\n", std::exp(iter->first), T[i]);
                }
                else if (!(eig0>0)){ //flags error if eig0 at a given T is not > 0
                    fmt_append(err_reactions2,"at T = {:.1f}\n", T[i]);
                }
            }
        }
        if (err_reactions1.size()) {
            throw InputFileError("LmrRate::validate", m_input,
                    "\nInvalid rate coefficient, k, for reaction '{}'\n{}",equation, to_string(err_reactions1));
        }
        else if (err_reactions2.size()) {
            throw InputFileError("LmrRate::validate", m_input,
                    "\nInvalid rate coefficient, eig0 (k at low-pressure limit), for reaction '{}'\n{}",equation, to_string(err_reactions2));
        }
    }

}

//Similar to updateFromStruct, evalFromStruct in PlogRate.h but with minor modifications
double LmrRate::speciesPlogRate(const LmrData& shared_data){ 
    if (logPeff_ > logP1_ && logPeff_ < logP2_) {
        return;
    }
    auto iter = pressures_s_.upper_bound(logPeff_);
    AssertThrowMsg(iter != pressures_s_.end(), "LmrRate::speciesPlogRate","Reduced-pressure out of range: {}", logPeff_);
    AssertThrowMsg(iter != pressures_s_.begin(), "LmrRate::speciesPlogRate","Reduced-pressure out of range: {}", logPeff_); 
    logP2_ = iter->first;
    ihigh1_ = iter->second.first;
    ihigh2_ = iter->second.second;
    logP1_ = (--iter)->first;
    ilow1_ = iter->second.first;
    ilow2_ = iter->second.second;
    rDeltaP_ = 1.0 / (logP2_ - logP1_);
    double log_k1, log_k2;
    if (ilow1_ == ilow2_) {
        log_k1 = rates_s_[ilow1_].evalLog(shared_data.logT, shared_data.recipT);
    } else {
        double k = 1e-300;
        for (size_t i = ilow1_; i < ilow2_; i++) {
            k += rates_s_[i].evalRate(shared_data.logT, shared_data.recipT);
        }
        log_k1 = std::log(k);
    }
    if (ihigh1_ == ihigh2_) {
        log_k2 = rates_s_[ihigh1_].evalLog(shared_data.logT, shared_data.recipT);
    } else {
        double k = 1e-300;
        for (size_t i = ihigh1_; i < ihigh2_; i++) {
            k += rates_s_[i].evalRate(shared_data.logT, shared_data.recipT);
        }
        log_k2 = std::log(k);
    }
    return exp(log_k1 + (log_k2-log_k1)*rDeltaP_*(logPeff_-logP1_));
}

double LmrRate::evalFromStruct(const LmrData& shared_data){
    double eig0_mix;
    double eig0_M=eig0_["M"].evalRate(shared_data.logT, shared_data.recipT);
    for (size_t i=0; i<shared_data.allSpecies.size(); i++){ //testing each species listed at the top of yaml file
        double Xi = shared_data.moleFractions[i];
        std::map<string, ArrheniusRate>::iterator it = eig0_.find(allSpecies[i]);
        if (it != eig0_.end()) {//key found, i.e. species has corresponding LMR data  
            eig0_mix += Xi*eig0_[allSpecies[i]].evalRate(shared_data.logT, shared_data.recipT);
        } else {//no LMR data for this species, so use M data as default
            eig0_mix += Xi*eig0_M;
        }
    }
    double log_eig0_mix = std::log(eig0_mix);
    for (size_t i=0; i<shared_data.allSpecies.size(); i++){ //testing each species listed at the top of yaml file
        double Xi = shared_data.moleFractions[i];
        double eig0; //eig0 val of a single species
        std::map<string, ArrheniusRate>::iterator it = eig0_.find(allSpecies[i]);
        if (it != eig0_.end()) {
            eig0 = eig0_[allSpecies[i]].evalRate(shared_data.logT, shared_data.recipT);
            pressures_s_=pressures_[allSpecies[i]];
            rates_s_=rates_[allSpecies[i]];
        } else {
            eig0 = eig0_M;
            pressures_s_=pressures_["M"];
            rates_s_=rates_["M"];
        }
        if (shared_data.logP != logP_) { //WHAT IS THE PURPOSE OF THIS STEP?
            logP_=shared_data.logP; 
            logPeff_=logP_+log_eig0_mix-log(eig0); //Peff is the effective pressure, formerly called "Ptilde"
            k_LMR_ += LmrRate::speciesPlogRate(shared_data)*eig0*Xi/eig0_mix; //Xtilde=eig0_s*Xi/eig0_mix
        }
    }
    return k_LMR_;
}
}

