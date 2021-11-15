#ifndef REACTION_DERIVATIVE_MANAGER_H
#define REACTION_DERIVATIVE_MANAGER_H

#include "cantera/zerodim.h"
#include "cantera/kinetics.h"
#include <vector>
#include <unordered_map>
#include <utility>

namespace Cantera
{

class ReactionDerivative
{
protected:
    // first element is always the independent variable
    std::vector<size_t> m_indices;
    std::vector<double> m_coeffs;
    size_t m_didx; // derivative index
    bool m_rev; // uses reverse reaction rate or not
    size_t m_rxn; // reaction idx
    double m_multiplier;

public:
    ReactionDerivative(size_t compLen, size_t dervIdx, size_t reactIdx, size_t* indices, double* coeffs, bool rev, double multiplier)
    {
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

    // Get derivatives
    inline virtual void __getDerivative(double rateConst, double * state, double* derivatives)
    {
        double derv = m_multiplier * rateConst * std::pow(state[m_indices[0]], m_coeffs[0]-1);
        for (size_t i = 1; i < m_indices.size(); i++)
        {
            derv *= std::pow(state[m_indices[i]], m_coeffs[i]);
        }
        derivatives[m_didx] += derv;
    };

    // Get derivatives
    inline virtual void __printDerivative(Reactor* reactor, double rateConst, double* state, double* derivatives)
    {
        std::cout<<"-------ReactionDerivative-------"<<std::endl;
        std::cout<<"m_rxn: "<<m_rxn<<std::endl;
        std::cout<<"m_multiplier: "<<m_multiplier<<std::endl;
        std::cout<<"rate constant:"<<rateConst<<std::endl;
        std::string kstr = m_rev ? "kr" : "kf";
        std::cout << m_multiplier * m_coeffs[0] << kstr << "_" << m_rxn;

        if (m_coeffs[0]-1 != 0)
        {
            std::cout<<"[" << reactor->componentName(m_indices[0]) << "]^" << "(" << m_coeffs[0]-1 << ")";
        }

        for (size_t i = 1; i < m_indices.size(); i++)
        {
            std::cout << "[" << reactor->componentName(m_indices[i]) << "]^" << "(" << m_coeffs[i]<<")";
        }
        std::cout<<std::endl;
    };

    inline size_t __reactionNumber(){return m_rxn;};

    inline bool __isRev(){return m_rev;};

    inline size_t __getDerivativeIndex(){return m_didx;};

    inline void __setDerivativeIndex(size_t nidx){m_didx = nidx;};

    inline void __ones(double* derivatives){derivatives[m_didx] = 1;};
};

class ReactionOneDerivative : public ReactionDerivative
{
public:
    ReactionOneDerivative(size_t compLen, size_t dervIdx, size_t reactIdx, size_t* indices, double* coeffs, bool rev, double multiplier) : ReactionDerivative(compLen, dervIdx, reactIdx, indices, coeffs, rev, multiplier){};
    ~ReactionOneDerivative(){};
    // Get derivatives
    inline void __getDerivative(double rateConst, double * state, double* derivatives)
    {
        double derv = m_multiplier * rateConst;
        for (size_t i = 1; i < m_indices.size(); i++)
        {
            derv *= state[m_indices[i]];
        }
        derivatives[m_didx] += derv;
    }
};

class ReactionDerivativeManager
{
private:
    std::vector<ReactionDerivative> m_reaction_derivatives;
    size_t m_column_size;
public:
    ReactionDerivativeManager(){};
    ~ReactionDerivativeManager(){};
    void initialize(Reactor* reactor);
    size_t getNumNonzeros();
    void getDerivatives(double* state, double* derivatives, double* kFwd, double* kRev);
    void getDerivativeIndices(std::vector<int>* indexVector);
    void remapDerivativeIndices(std::unordered_map<int, int>* indexMap);
    void addDerivative(bool isOne, size_t compLen, size_t dervIdx, size_t reactIdx, size_t* indices, double* coeffs, bool rev, double multiplier);
    bool checkOneCoeffs(Composition currComp);
    void derivativeLoop(size_t nReaction, Reactor* reactor, Composition* c1, Composition* c2, Composition* c3, bool isOne, bool rev);
    void printDerivativeContributions(Reactor* reactor, double* state, double* derivatives, double* kFwd, double* kRev);

    // Debugging function to print indices
    inline void printIndices()
    {
        for (size_t i = 0; i < m_reaction_derivatives.size(); i++)
        {
            std::cout<<m_reaction_derivatives[i].__getDerivativeIndex()<<std::endl;
        }
    }

    // Debugging function to set ones
    inline void setOnes(double* derivatives)
    {
        for (size_t i = 0; i < m_reaction_derivatives.size(); i++)
        {
            m_reaction_derivatives[i].__ones(derivatives);
        }
    }
};

}
#endif
