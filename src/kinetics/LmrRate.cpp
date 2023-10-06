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

// void LmrRate::getParameters(AnyMap& rateNode, const Units& rate_units) const{
//     vector<AnyMap> rateList;
//     if (!valid()) {
//         // object not fully set up
//         return;
//     }
//     std::multimap<double, ArrheniusRate> rateMap;
//     // initial preincrement to skip rate for P --> 0
//     for (auto iter = ++pressures_.begin();
//             iter->first < 1000; // skip rates for (P --> infinity)
//             ++iter) {
//         for (size_t i = iter->second.first; i < iter->second.second; i++) {
//             rateMap.insert({std::exp(iter->first), rates_[i]});
//         }
//     }
//     for (const auto& [pressure, rate] : rateMap) {
//         AnyMap rateNode_;
//         rateNode_["P"].setQuantity(pressure, "Pa");
//         rate.getRateParameters(rateNode_);
//         rateList.push_back(std::move(rateNode_));
//     }
//     rateNode["rate-constants"] = std::move(rateList);
// }

void LmrRate::validate(const string& equation, const Kinetics& kin){
    if (!valid()) {
        throw InputFileError("LmrRate::validate", m_input,
            "Rate object for reaction '{}' is not configured.", equation);
    }
    fmt::memory_buffer err_reactions;
    double T[] = {300.0, 500.0, 1000.0, 2000.0, 5000.0, 10000.0};
    LmrData data;
    // Iterate through the outer map (string to inner map)
    for (const auto& outer_pair : pressures_) { // recall that pressures_ is a map within a map
        const std::string& specieskey = outer_pair.first; //specieskey refers only to the species for which LMR data is provided in yaml (e.g. 'H2O', 'M')
        const std::map<double, std::pair<size_t, size_t>>& inner_map = outer_pair.second;
        vector<ArrheniusRate> rates_i_ = rates_[specieskey];
        for (auto iter = ++inner_map.begin(); iter->first < 1000; iter++) {
            data.update(T[0], exp(iter->first));
            ilow1_ = iter->second.first;
            ilow2_ = iter->second.second;       
            for (size_t i=0; i < 6; i++) {
                double k = 0;
                for (size_t p = ilow1_; p < ilow2_; p++) {
                    k += rates_i_.at(p).evalRate(log(T[i]), 1.0 / T[i]);
                }
                if (!(k > 0)) {
                    fmt_append(err_reactions,
                        "at P = {:.5g}, T = {:.1f}\n", std::exp(iter->first), T[i]);
                }
            }
        }
        if (err_reactions.size()) {
            throw InputFileError("LmrRate::validate", m_input,
                "\nInvalid rate coefficient for reaction '{}'\n{}",
                equation, to_string(err_reactions));
        }
    }
}

double LmrRate::computeSpeciesRate(const LmrData& shared_data){
    double eig0 = eig0_[s_].evalRate(shared_data.logT, shared_data.recipT);
    double Xtilde=eig0*Xs_/eig0_mix_; //DOESN'T WORK YET BC WE DON'T HAVE MF DATA FOR EACH SPECIES
    //STILL NEED TO ACCOUND FOR THE REDUCED PRESSURE, P_i *************************************************************************
    if (shared_data.logP != logP_) { 
        logP_=log(shared_data.P*eig0_mix_/eig0); //need to use natl log instead to align with code later on
        if (logP_ > logP1_ && logP_ < logP2_) {
            return;
        }
        auto iter = pressures_[s_].upper_bound(logP_); //locate the pressure of interest
        AssertThrowMsg(iter != pressures_[s_].end(), "LmrRate::evalFromStruct",
                        "Pressure out of range: {}", logP_);
        AssertThrowMsg(iter != pressures_[s_].begin(), "LmrRate::evalFromStruct",
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
            log_k1 = rates_[s_][ilow1_].evalLog(shared_data.logT, shared_data.recipT);
        } else {
            double k = 1e-300; // non-zero to make log(k) finite
            for (size_t i = ilow1_; i < ilow2_; i++) {
                k += rates_[s_][i].evalRate(shared_data.logT, shared_data.recipT);
            }
            log_k1 = std::log(k);
        }
        if (ihigh1_ == ihigh2_) {
            log_k2 = rates_[s_][ihigh1_].evalLog(shared_data.logT, shared_data.recipT);
        } else {
            double k = 1e-300; // non-zero to make log(k) finite
            for (size_t i = ihigh1_; i < ihigh2_; i++) {
                k += rates_[s_][i].evalRate(shared_data.logT, shared_data.recipT);
            }
            log_k2 = std::log(k);
        }
        return Xtilde*exp(log_k1 + (log_k2-log_k1)/(logP2_-logP1_) * (logP_-logP1_+log(eig0_mix_)-log(eig0)));
    }
}

double LmrRate::evalFromStruct(const LmrData& shared_data){
    for (size_t i=0; i<shared_data.allSpecies.size(); i++){ //testing each species listed at the top of yaml file
        Xs_ = shared_data.moleFractions[i];
        s_ = allSpecies[i];
        std::map<string, ArrheniusRate>::iterator it = eig0_.find(s_);
        if (it != eig0_.end()) {//key found, i.e. s has corresponding LMR data  
            eig0_mix_ += Xs_*eig0_[s_].evalRate(shared_data.logT, shared_data.recipT); 
        } else {//no LMR data for this species s, so use M data as default
            s_ = "M";
            eig0_mix_ += Xs_*eig0_[s_].evalRate(shared_data.logT, shared_data.recipT);
        }
    }
    for (size_t i=0; i<shared_data.allSpecies.size(); i++){ //testing each species listed at the top of yaml file
        Xs_ = shared_data.moleFractions[i];
        s_ = allSpecies[i];
        std::map<string, ArrheniusRate>::iterator it = eig0_.find(s_);
        if (it != eig0_.end()) {//key found, i.e. s has corresponding LMR data  
            k_LMR_ += LmrRate::computeSpeciesRate(shared_data);
        } else {//no LMR data for this species s, so use M data as default
            s_ = "M";
            k_LMR_ += LmrRate::computeSpeciesRate(shared_data);
        }
    }
    return k_LMR_;
}
}

