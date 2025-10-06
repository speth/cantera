/*
 * Property calculation demo
 * =========================
 *
 * This demonstration program builds an object representing a reacting gas
 * mixture, and uses it to compute thermodynamic properties, chemical
 * equilibrium, and transport properties.
 *
 * .. tags:: C++, tutorial, thermodynamics, equilibrium, transport
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

// Include cantera header files. They should be included in the form
// "cantera/*.h". These headers are designed for use in C++ programs and
// provide a simplified interface to the Cantera header files. If you need
// to include core headers directly, use the format "cantera/module/*.h".

#include "cantera/core.h"
#include "cantera/kinetics/Reaction.h"
#include "cantera/base/global.h" // provides Cantera::writelog
#include <iostream>
#include <span>

// All Cantera kernel names are in namespace Cantera. You can either
// reference everything as Cantera::<name>, or include the following
// 'using namespace' line.
using namespace Cantera;
using std::span;

// The program is put into a function so that error handling code can
// be conveniently put around the whole thing. See main() below.

int f(span<double> x)
{
    int ret = 0;
    for (size_t i = 0; i < 13; i++) {
        writelog("    x[{}] = {}\n", i, x[i]);
        if (x[i] > 1.0) {
            ret = i;
        }
    }
    return ret;
}

void demoprog()
{
    writelog("\n**** C++ Test Program ****\n");



    auto sol = newSolution("h2o2.yaml", "ohmech");
    auto gas = sol->thermo();
    double temp = 1200.0;
    double pres = OneAtm;
    gas->setState_TPX(temp, pres, "H2:1, O2:1, AR:2");

    size_t K = gas->nSpecies();

    double* ss = new double[2*K];
    gas->getMoleFractions(ss);
    writelog("f(ss) = {}\n", f(span<double>(ss, K)));
    delete[] ss;

    vector<double> yy(K);
    gas->getPartialMolarCp(yy.data());
    writelog("f(yy) = {}\n", f(yy));

    Eigen::ArrayXd ee(K);
    gas->getPartialMolarEnthalpies(ee.data());
    writelog("f(ee) = {}\n", f(ee)); // Works starting in Eigen 5.0.0
    writelog("f(ee) = {}\n", f(span<double>(ee.data(), ee.data() + ee.size()))); // Works for all Eigen

    // Thermodynamic properties

    writelog("\n\nInitial state:\n\n");
    writelog(
        "Temperature:    {:14.5g} K\n"
        "Pressure:       {:14.5g} Pa\n"
        "Density:        {:14.5g} kg/m3\n"
        "Molar Enthalpy: {:14.5g} J/kmol\n"
        "Molar Entropy:  {:14.5g} J/kmol-K\n"
        "Molar cp:       {:14.5g} J/kmol-K\n",
        gas->temperature(), gas->pressure(), gas->density(),
        gas->enthalpy_mole(), gas->entropy_mole(), gas->cp_mole());

    // set the gas to the equilibrium state with the same specific
    // enthalpy and pressure
    gas->equilibrate("HP");

    writelog("\n\nEquilibrium state:\n\n");
    writelog(
        "Temperature:    {:14.5g} K\n"
        "Pressure:       {:14.5g} Pa\n"
        "Density:        {:14.5g} kg/m3\n"
        "Molar Enthalpy: {:14.5g} J/kmol\n"
        "Molar Entropy:  {:14.5g} J/kmol-K\n"
        "Molar cp:       {:14.5g} J/kmol-K\n",
        gas->temperature(), gas->pressure(), gas->density(),
        gas->enthalpy_mole(), gas->entropy_mole(), gas->cp_mole());


    // Reaction information

    auto kin = sol->kinetics();
    int irxns = kin->nReactions();
    vector<double> qf(irxns);
    vector<double> qr(irxns);
    vector<double> q(irxns);

    // since the gas has been set to an equilibrium state, the forward
    // and reverse rates of progress should be equal for all
    // reversible reactions, and the net rates should be zero.
    kin->getFwdRatesOfProgress(&qf[0]);
    kin->getRevRatesOfProgress(&qr[0]);
    kin->getNetRatesOfProgress(&q[0]);

    writelog("\n\n");
    for (int i = 0; i < irxns; i++) {
        const auto& rxn = kin->reaction(i);
        writelog("{:30s} {:14.5g} {:14.5g} {:14.5g}  kmol/m3/s\n",
               rxn->equation(), qf[i], qr[i], q[i]);
    }


    // transport properties

    // create a transport manager for the gas that computes
    // mixture-averaged properties
    sol->setTransportModel("mixture-averaged");
    auto tr = sol->transport();

    // print the viscosity, thermal conductivity, and diffusion
    // coefficients
    writelog("\n\nViscosity:            {:14.5g} Pa-s\n", tr->viscosity());
    writelog("Thermal conductivity: {:14.5g} W/m/K\n", tr->thermalConductivity());

    int nsp = gas->nSpecies();
    vector<double> diff(nsp);
    tr->getMixDiffCoeffs(&diff[0]);
    int k;
    writelog("\n\n{:20s}  {:21s}\n", "Species", "Diffusion Coefficient");
    for (k = 0; k < nsp; k++) {
        writelog("{:20s}  {:14.5g} m2/s\n", gas->speciesName(k), diff[k]);
    }
}



int main()
{
    try {
        demoprog();
        return 0;
    } catch (CanteraError& err) {
        std::cout << err.what() << std::endl;
        return -1;
    }
}
