#include "cantera/kinetics/ReactionDerivativeManager.h"
#include "cantera/zeroD/Reactor.h"

using namespace Cantera;

void ReactionDerivativeManager::getDerivativeIndices(vector_int* indexVector)
{
    for (size_t i = 0; i < m_reaction_derivatives.size(); i++)
    {
        indexVector->push_back(m_reaction_derivatives[i].getDerivativeIndex());
    }
}

void ReactionDerivativeManager::getDerivatives(double* state, double* derivatives, double* kFwd, double* kRev)
{
    for (size_t i = 0; i < m_reaction_derivatives.size(); i++)
    {
        size_t rateIdx = m_reaction_derivatives[i].reactionNumber();
        double rateConst = m_reaction_derivatives[i].isRev() ? kRev[rateIdx] : kFwd[rateIdx];
        m_reaction_derivatives[i].getDerivative(rateConst, state, derivatives);
    }
}

void ReactionDerivativeManager::addToDerivatives(double* derivatives, double* additions)
{
    for (size_t i = 0; i < m_reaction_derivatives.size(); i++)
    {
        m_reaction_derivatives[i].addToDerivative(derivatives, additions);
    }
}

void ReactionDerivativeManager::remapDerivativeIndices(std::unordered_map<int, int>* indexMap, bool warnRemapped)
{
    if (m_remapped)
    {
        for (size_t i = 0; i < m_reaction_derivatives.size(); i++)
        {
            int newIndex = m_reaction_derivatives[i].getDerivativeIndex();
            newIndex = indexMap->at(newIndex);
            m_reaction_derivatives[i].setDerivativeIndex(newIndex);
        }
        m_remapped = false;
    }
    else if (warnRemapped)
    {
        warn_user("ReactionDerivativeManager::remapDerivativeIndices","ReactionDerivativeManager has already been remapped, use ReactionDerivativeManager::reset() and reinitialize to remap again.");
    }

}

void ReactionDerivativeManager::addDerivative(bool isOne, size_t compLen, size_t dervIdx, size_t reactIdx, size_t* indices, double* coeffs, bool rev, double multiplier)
{
    if (isOne)
    {
        m_reaction_derivatives.push_back(ReactionDerivative(compLen, dervIdx, reactIdx, indices, coeffs, rev, multiplier));
    }
    else
    {
        m_reaction_derivatives.push_back(ReactionOneDerivative(compLen, dervIdx, reactIdx, indices, coeffs, rev, multiplier));
    }
}

bool ReactionDerivativeManager::checkOneCoeffs(Composition currComp)
{
    bool isOne = true;
    for (auto i = currComp.begin(); i != currComp.end(); i++)
    {
        if (i->second != 1)
        {
            isOne = false;
            break;
        }
    }
    return isOne;
}

void ReactionDerivativeManager::derivativeLoop(size_t nReaction, Reactor& reactor, Composition* c1, Composition* c2, Composition* c3, bool isOne, bool rev)
{
    for (auto j = c1->begin(); j != c1->end(); j++) // column loop (independent)
    {
        bool contained = c2->end() != (c2->find(j->first));
        if (contained)
        {
            size_t col = reactor.componentIndex(j->first);
            std::vector<size_t> indices(1, col);
            std::vector<double> coeffs(1, j->second);
            // Find derivative indices
            for (auto k = c2->begin(); k != c2->end(); k++) // row loop (dependent)
            {
                if (j->first != k->first)
                {
                    indices.push_back(reactor.componentIndex(k->first));
                    coeffs.push_back(k->second);
                }
            }
            // Add derivative contributions of c2
            for (auto k = c2->begin(); k != c2->end(); k++) // row loop (dependent)
            {
                size_t didx = reactor.componentIndex(k->first) + col * m_column_size;
                double multiplier = -j->second * k->second;
                addDerivative(isOne, indices.size(), didx, nReaction, indices.data(), coeffs.data(), rev, multiplier);
            }
            // Add derivative contributions of c3
            for (auto k = c3->begin(); k != c3->end(); k++) // row loop (dependent)
            {
                size_t didx = reactor.componentIndex(k->first) + col * m_column_size;
                double multiplier = j->second * k->second;
                addDerivative(isOne, indices.size(), didx, nReaction, indices.data(), coeffs.data(), rev, multiplier);
            }
        }
    }
}

void ReactionDerivativeManager::initialize(Reactor& reactor)
{
    if (m_init)
    {
        auto kin = reactor.getKineticsMgr();
        size_t nreactions = kin->nReactions();
        size_t nspecies = kin->nTotalSpecies();
        m_reaction_derivatives.reserve(nreactions * nspecies);
        m_column_size = reactor.neq();
        for (size_t i = 0; i < nreactions; i++)
        {
            //current reaction
            auto currReaction = kin->reaction(i);
            // get reaction components
            Composition reactants = currReaction->reactants;
            Composition products = currReaction->products;
            Composition allSpecies;
            allSpecies.insert(reactants.begin(), reactants.end());
            allSpecies.insert(products.begin(), products.end());
            bool isOne = checkOneCoeffs(reactants) && checkOneCoeffs(products);

            derivativeLoop(i, reactor, &allSpecies, &reactants, &products, isOne, false);
            // Reversible components
            if (currReaction->reversible)
            {
                derivativeLoop(i, reactor, &allSpecies, &products, &reactants, isOne, true);
            }
        }
        // get number of unique nonzero elements
        m_reaction_derivatives.shrink_to_fit();
        vector_int indexVector;
        getDerivativeIndices(&indexVector);
        std::sort(indexVector.begin(), indexVector.end(), std::less<int>());
        vector_int::iterator unique_itr = std::unique(indexVector.begin(), indexVector.end());
        indexVector.resize( std::distance(indexVector.begin(), unique_itr));
        indexVector.shrink_to_fit();
        m_nnz = indexVector.size();
        m_init = false;
    }
}

size_t ReactionDerivativeManager::getNumNonzeros()
{
    return m_nnz;
}

void ReactionDerivativeManager::reset()
{
    m_reaction_derivatives.clear();
    m_column_size = 0;
    m_init = true;
    m_remapped = true;
}
