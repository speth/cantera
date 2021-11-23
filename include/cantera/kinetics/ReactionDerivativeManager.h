/**
 * @file ReactionDerivativeManager.h this file contains classes used in
 *  determining rate of progress derivatives and is currently
 *  implemented through the reactor system.
 */

// This file is part of Cantera. See License.txt in the top-level
// directory or at https://cantera.org/license.txt for license and
// copyright information.

#ifndef REACTION_DERIVATIVE_MANAGER_H
#define REACTION_DERIVATIVE_MANAGER_H


#include "cantera/kinetics/Kinetics.h"
#include <vector>
#include <unordered_map>
#include <utility>
#include "cantera/base/ct_defs.h"

namespace Cantera
{
//! Forward Declarations
class Reactor;
//! ReactionDerivative represents a single derivative associated with
//! each chemical reaction. These derivatives are managed by the
//! ReactionDerivativeManager class
class ReactionDerivative
{
public:
    ReactionDerivative(size_t compLen, size_t dervIdx, size_t reactIdx, size_t* indices, double* coeffs, bool rev, double multiplier){
        m_rev = rev;
        m_rxn = reactIdx;
        m_didx = dervIdx;
        m_multiplier = multiplier;
        m_indices.resize(compLen);
        m_coeffs.resize(compLen);
        for (size_t i = 0; i < compLen; i++)
        {
            m_indices[i] = indices[i];
            m_coeffs[i] = coeffs[i];
        }
    };
    ~ReactionDerivative(){};
    //! Use this function to get the value of the derivative written to
    //! a buffer
    virtual void getDerivative(double rateConst, double * state, double* derivatives){
        double derv = m_multiplier * rateConst * std::pow(state[m_indices[0]], m_coeffs[0]-1);
        for (size_t i = 1; i < m_indices.size(); i++)
        {
            derv *= std::pow(state[m_indices[i]], m_coeffs[i]);
        }
        derivatives[m_didx] += derv;
    };
    //! Use this function to get the associated reaction number
    size_t reactionNumber(){return m_rxn;};
    //! Use this function to determine if the associated reaction is
    //! reversible
    bool isRev(){return m_rev;};
    //! Use this function to get the index of the derivative
    size_t getDerivativeIndex(){return m_didx;};
    //! Use this function to remap the index of the derivative
    void setDerivativeIndex(size_t nidx){m_didx = nidx;};
protected:
    // first element is always the independent variable
    std::vector<size_t> m_indices;
    // coefficients used in derivatives
    std::vector<double> m_coeffs;
    size_t m_didx; // derivative index
    bool m_rev; // uses reverse reaction rate or not
    size_t m_rxn; // reaction idx
    double m_multiplier;
};
//! ReactionOneDerivative is a ReactionDerivative with purely ones for
//! coefficients. This makes the calculation of the derivatives more
//! efficient.
class ReactionOneDerivative : public ReactionDerivative
{
public:
    using ReactionDerivative::ReactionDerivative;
    ~ReactionOneDerivative(){};
    //! Use this function to get the value of the derivative written to
    //! a buffer
    void getDerivative(double rateConst, double * state, double* derivatives){
        double derv = m_multiplier * rateConst;
        for (size_t i = 1; i < m_indices.size(); i++)
        {
            derv *= state[m_indices[i]];
        }
        derivatives[m_didx] += derv;
    };
};
//! Use this class to manager and interface with all reaction
//! derivatives associated with a mechanism or reactor.
class ReactionDerivativeManager
{
public:
    ReactionDerivativeManager(){};
    ~ReactionDerivativeManager(){};
    //! Use this to generate all reaction derivatives
    void initialize(Reactor& reactor);
    //! Use this to return the number of nonzero derivatives
    size_t getNumNonzeros();
    //! Use this function to obtain all derivatives writting to the
    //! derivatives buffer
    void getDerivatives(double* state, double* derivatives, double* kFwd, double* kRev);
    //! Use this function to obtain all derivative indices
    void getDerivativeIndices(vector_int* indexVector);
    //! Use this function to remap derivative indices based on the given
    //! indexMap
    void remapDerivativeIndices(std::unordered_map<int, int>* indexMap, bool warnRemapped=true);
    //! Use this function to reset this object which requires
    //! reinitialization
    void reset();
    //! Check if this object has been initialized
    bool isInit(){return !m_init;};
protected:
    //! Adds a ReactionDerivative element to vector
    void addDerivative(bool isOne, size_t compLen, size_t dervIdx, size_t reactIdx, size_t* indices, double* coeffs, bool rev, double multiplier);
    //! Checks reaction for all ones as coefficients
    bool checkOneCoeffs(Composition currComp);
    //! Loop to go through all reactions
    void derivativeLoop(size_t nReaction, Reactor& reactor, Composition* c1, Composition* c2, Composition* c3, bool isOne, bool rev);
    //! Vector storing all reaction derivatives
    std::vector<ReactionDerivative> m_reaction_derivatives;
    //! Jacobian species dimensions
    size_t m_column_size;
    //! Says whether object needs initialized or not
    bool m_init = true;
    //! Says whether object has been remapped or not
    bool m_remapped = true;
    //! Number of nonzero elements
    size_t m_nnz = 0;
};

}
#endif
