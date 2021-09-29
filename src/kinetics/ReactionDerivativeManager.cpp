#include "cantera/kinetics/ReactionDerivativeManager.h"

using namespace Cantera;

void ReactionDerivativeManager::getDerivativeIndices(std::vector<int>* indexVector)
{
    for (size_t i = 0; i < m_reaction_derivatives.size(); i++)
    {
        indexVector->push_back(m_reaction_derivatives[i].__getDerivativeIndex());
    }
}

void ReactionDerivativeManager::getDerivatives(double* state, double* derivatives, double* kFwd, double* kRev)
{
    for (size_t i = 0; i < m_reaction_derivatives.size(); i++)
    {
        size_t rateIdx = m_reaction_derivatives[i].__reactionNumber();
        double rateConst = m_reaction_derivatives[i].__isRev() ? kRev[rateIdx] : kFwd[rateIdx];
        m_reaction_derivatives[i].__getDerivative(rateConst, state, derivatives);
    }
}

void ReactionDerivativeManager::remapDerivativeIndices(std::unordered_map<int, int>* indexMap)
{
    for (size_t i = 0; i < m_reaction_derivatives.size(); i++)
    {
        int newIndex = m_reaction_derivatives[i].__getDerivativeIndex();
        newIndex = indexMap->at(newIndex);
        m_reaction_derivatives[i].__setDerivativeIndex(newIndex);
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

void ReactionDerivativeManager::derivativeLoop(size_t nReaction, Reactor* reactor, Composition* c1, Composition* c2, Composition* c3, bool isOne, bool rev)
{
    for (auto j = c1->begin(); j != c1->end(); j++) // column loop (independent)
    {
        bool contained = c2->end() != (c2->find(j->first));
        if (contained)
        {
            size_t col = reactor->componentIndex(j->first);
            std::vector<size_t> indices(1, col);
            std::vector<double> coeffs(1, j->second);
            // Find derivative indices
            for (auto k = c2->begin(); k != c2->end(); k++) // row loop (dependent)
            {
                if (j->first != k->first)
                {
                    indices.push_back(reactor->componentIndex(k->first));
                    coeffs.push_back(k->second);
                }
            }
            // Add derivative contributions of c2
            for (auto k = c2->begin(); k != c2->end(); k++) // row loop (dependent)
            {
                size_t didx = reactor->componentIndex(k->first) + col * m_column_size;
                double multiplier = -j->second * k->second;
                addDerivative(isOne, indices.size(), didx, nReaction, indices.data(), coeffs.data(), rev, multiplier);
            }
            // Add derivative contributions of c3
            for (auto k = c3->begin(); k != c3->end(); k++) // row loop (dependent)
            {
                size_t didx = reactor->componentIndex(k->first) + col * m_column_size;
                double multiplier = j->second * k->second;
                addDerivative(isOne, indices.size(), didx, nReaction, indices.data(), coeffs.data(), rev, multiplier);
            }
        }
    }
}

void ReactionDerivativeManager::initialize(Reactor* reactor)
{
    auto kin = reactor->getKineticsMgr();
    size_t nreactions = kin->nReactions();
    size_t nspecies = kin->nTotalSpecies();
    m_reaction_derivatives.reserve(nreactions * nspecies);
    m_column_size = reactor->neq();
    for (size_t i = 0; i < nreactions; i++)
    {
        //current reaction
        auto currReaction = kin->getReaction(i);
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
}

void ReactionDerivativeManager::printDerivativeContributions(Reactor* reactor, double* state, double* derivatives, double* kFwd, double* kRev)
{
    for (size_t i = 0; i < m_reaction_derivatives.size(); i++)
    {
        size_t rateIdx = m_reaction_derivatives[i].__reactionNumber();
        double rateConst = m_reaction_derivatives[i].__isRev() ? kRev[rateIdx] : kFwd[rateIdx];
        m_reaction_derivatives[i].__printDerivative(reactor, rateConst, state, derivatives);
    }
}
