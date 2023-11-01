//! @file LmrRate.cpp
// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/LmrRate.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/kinetics/Kinetics.h"

namespace Cantera{

LmrData::LmrData(){ //THIS METHOD WAS ADAPTED SOMEWHAT BLINDLY FROM FALLOFF.CPP, PLEASE VERIFY IF CORRECT
    moleFractions.resize(1, NAN);
}

bool LmrData::update(const ThermoPhase& phase, const Kinetics& kin){
    double T = phase.temperature();
    double P = phase.pressure(); //find out what units this is in
    int X = phase.stateMFNumber();
    if (P != pressure || T != temperature || X != mfNumber) {
        ReactionData::update(T);
        pressure = P;
        logP = std::log(P);
        mfNumber=X;
        phase.getMoleFractions(moleFractions.data()); //IS THIS THE CORRECT WAY TO UPDATE THE MOLE FRACTION VECTOR
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

void LmrRate::validate(const string& equation, const Kinetics& kin){
    //Get the list of all species in yaml (not just the ones for which LMRR data exists)
    ThermoPhase::Phase phase;
    allSpecies_ = phase.speciesNames();
    //Validate the LMRR input data for each species
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

double LmrRate::speciesPlogRate(const LmrData& shared_data){ 
    if (logPeff_ > logP1_ && logPeff_ < logP2_) {
        // return 10;
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
    // double val = exp(log_k1 + (log_k2-log_k1)*rDeltaP_*(logPeff_-logP1_));
    return exp(log_k1 + (log_k2-log_k1)*rDeltaP_*(logPeff_-logP1_));
}

double LmrRate::evalFromStruct(const LmrData& shared_data){
    double eig0_mix;
    double eig0_M=eig0_["M"].evalRate(shared_data.logT, shared_data.recipT);
    for (size_t i=0; i<allSpecies_.size(); i++){ //testing each species listed at the top of yaml file
        double Xi = shared_data.moleFractions[i];
        std::map<string, ArrheniusRate>::iterator it = eig0_.find(allSpecies_[i]);
        if (it != eig0_.end()) {//key found, i.e. species has corresponding LMR data  
            eig0_mix += Xi*eig0_[allSpecies_[i]].evalRate(shared_data.logT, shared_data.recipT);
        } else {//no LMR data for this species, so use M data as default
            eig0_mix += Xi*eig0_M;
        }
    }
    double log_eig0_mix = std::log(eig0_mix);
    for (size_t i=0; i<allSpecies_.size(); i++){ //testing each species listed at the top of yaml file
        double Xi = shared_data.moleFractions[i];
        double eig0; //eig0 val of a single species
        std::map<string, ArrheniusRate>::iterator it = eig0_.find(allSpecies_[i]);
        if (it != eig0_.end()) {
            eig0 = eig0_[allSpecies_[i]].evalRate(shared_data.logT, shared_data.recipT);
            pressures_s_=pressures_[allSpecies_[i]];
            rates_s_=rates_[allSpecies_[i]];
        } else {
            eig0 = eig0_M;
            pressures_s_=pressures_["M"];
            rates_s_=rates_["M"];
        }
        if (shared_data.logP != logP_) { //WHAT IS THE PURPOSE OF THIS STEP?
            logP_=shared_data.logP; 
            logPeff_=logP_+log_eig0_mix-log(eig0); //Peff is the effective pressure, formerly called "Ptilde"
            k_LMR_ += LmrRate::speciesPlogRate(shared_data)*eig0*Xi/eig0_mix; //Xtilde=eig0_s*Xi/eig0_mix
            //k_LMR_ += 1; //placeholder, just used as a test
        }
    }
    return k_LMR_;
}

void LmrRate::getParameters(AnyMap& rateNode, const Units& rate_units) const{
    //Create copies of existing variables
    map<string, map<double, pair<size_t, size_t>>> pressures__ = pressures_;
    map<string, vector<ArrheniusRate>> rates__ = rates_;
    map<string,ArrheniusRate> eig0__ = eig0_;
    vector<string> colliderList;
    for (const auto& entry : eig0__) {
        colliderList.push_back(entry.first);
    }
    vector<AnyMap> topLevelList;
    for (size_t i=0; i<colliderList.size(); i++) {
        AnyMap speciesNode; //will be filled with all LMR data for a single collider (species)
        if (!valid()) { //WHERE TO PUT THIS CONDITIONAL STATEMENT?
            return;
        }
        string& s = colliderList[i];
        //1) Save name of species to "name"
        speciesNode["name"]=s;
        //2) Save single set of arrhenius params for eig0 to "low-P-rate-constant"
        AnyMap tempNode;
        eig0__[s].getRateParameters(tempNode);
        if (!tempNode.empty()){
            speciesNode["low-P-rate-constant"]=std::move(tempNode);
        }
        tempNode.clear();
        //3) Save list of rate constant params to "rate-constants"
        std::multimap<double, ArrheniusRate> rateMap;
        for (auto iter = ++pressures__[s].begin(); iter->first < 1000; ++iter) {
            for (size_t i = iter->second.first; i < iter->second.second; i++) {
                rateMap.insert({std::exp(iter->first), rates__[s][i]});
            }
        }
        vector<AnyMap> rateList;
        for (const auto& [pressure, rate] : rateMap) {
            tempNode["P"].setQuantity(pressure, "atm");
            rate.getRateParameters(tempNode);
            rateList.push_back(std::move(tempNode));
            tempNode.clear();
        }    
        speciesNode["rate-constants"] = std::move(rateList);
        topLevelList.push_back(std::move(speciesNode));
    }
    rateNode["collider-list"]=std::move(topLevelList);
}
}

