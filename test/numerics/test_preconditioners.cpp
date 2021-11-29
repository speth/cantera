#include "gtest/gtest.h"
#include "cantera/base/Solution.h"
#include "cantera/thermo.h"
#include "cantera/kinetics.h"
#include "cantera/kinetics/GasKinetics.h"
#include "cantera/transport.h"
#include "cantera/numerics/AdaptivePreconditioner.h"
#include "cantera/zerodim.h"
#include <stdlib.h>
#include <time.h>
#include <queue>

using namespace Cantera;

void twoStepManualPrecondition(AdaptivePreconditioner& externalPrecon, MoleReactor& reactor, double* N, double* Ndot, double t=0.0, size_t loc=0)
    {
        auto thermo = reactor.getThermoMgr();
        auto kinetics = reactor.getKineticsMgr();
        double volInv = 1/reactor.volume();
        // Getting data for manual calculations
        size_t numberOfReactions = kinetics->nReactions();
        // vectors for manual calcs
        vector_fp kf(numberOfReactions, 0.0);
        vector_fp kr(numberOfReactions, 0.0);
        // getting actual data
        vector_fp C(thermo->nSpecies(), 0.0);
        thermo->getConcentrations(C.data());
        kinetics->getFwdRateConstants(kf.data());
        kinetics->getRevRateConstants(kr.data());
        // Assign concentrations to species
        double O2 = C[0];
        double CH4 = C[2];
        double CO = C[4];
        // Setting elements
        // O2
        externalPrecon(1 + loc, 1 + loc) = (- 2.25 * kf[0] * CH4 * std::pow(O2, 0.5) - 0.25 * kf[1] * CO * std::pow(O2, -0.5)) * volInv; // dO2/dO2
        externalPrecon(1 + loc, 3 + loc) = - 1.5 * kf[0] * std::pow(O2, 1.5) * volInv; // dO2/dCH4
        externalPrecon(1 + loc, 4 + loc) = 0.5 * kr[1] * volInv; // dO2/dCO2
        externalPrecon(1 + loc, 5 + loc) = - 0.5 * kf[1] * std::pow(O2, 0.5) * volInv; // dO2/dCO
        // H2O
        externalPrecon(2 + loc, 1 + loc) = 3 * kf[0] * CH4 * std::pow(O2, 0.5) * volInv; // dH2O/dO2
        externalPrecon(2 + loc, 3 + loc) = 2 * kf[0] * std::pow(O2, 1.5) * volInv; // dH2O/dCH4
        // CH4
        externalPrecon(3 + loc, 1 + loc) = - 1.5 * kf[0] * CH4 * std::pow(O2, 0.5) * volInv; // dCH4/dO2
        externalPrecon(3 + loc, 3 + loc) = - kf[0] * std::pow(O2, 1.5) * volInv; // dCH4/dCH4
        // CO2
        externalPrecon(4 + loc, 1 + loc) = (0.5 * kf[1] * CO * std::pow(O2, -0.5)) * volInv; // dCO2/dO2
        externalPrecon(4 + loc, 4 + loc) = -kr[1] * volInv; // dCO2/dCO2
        externalPrecon(4 + loc, 5 + loc) = kf[1] * std::pow(O2, 0.5) * volInv; // dCO2/CO
        //CO
        externalPrecon(5 + loc, 1 + loc) = (1.5 * kf[0] * CH4 * std::pow(O2, 0.5) - 0.5 * kf[1] * CO * std::pow(O2, -0.5)) * volInv; // dCO/dO2
        externalPrecon(5 + loc, 3 + loc) = kf[0] * std::pow(O2, 1.5) * volInv; // dCO/dCH4
        externalPrecon(5 + loc, 4 + loc) = kr[1] * volInv; // dCO/dCO2
        externalPrecon(5 + loc, 5 + loc) = -kf[1] * std::pow(O2,0.5) * volInv; // dCO/CO
        // temperature derivatives
        if (reactor.energyEnabled())
        {
            reactor.StateDerivatives(externalPrecon, t, N, Ndot, nullptr);
        }
    };

void preconditionWithClass(AdaptivePreconditioner& internalPrecon, IdealGasConstPressureMoleReactor& reactor, double* N, double* Ndot, double t=0.0)
{
    // class based
    internalPrecon.reset();
    reactor.preconditionerSetup(internalPrecon, t, N, Ndot, nullptr);
    internalPrecon.setup();
}

void preconditionerInitialize(ReactorNet& network, AdaptivePreconditioner& precon, size_t reactorStart=0, double gamma=1.0, double thresh=0)
{
    precon.initialize(network);
    precon.setGamma(gamma);
    precon.setThreshold(thresh);
}

void TwoStepMechanism(bool energy)
{
    // Constants
    double volume = 0.2;
    // Setting up solution object and thermo/kinetics pointers
    auto sol = newSolution("methane_twostep.yaml");
    auto thermo = sol->thermo();
    auto kinetics = sol->kinetics();
    thermo->setEquivalenceRatio(1, "CH4", "O2:1");
    thermo->setState_TP(1000, 101325);
    // Set up reactor object
    IdealGasConstPressureMoleReactor reactor;
    reactor.insert(sol);
    reactor.setInitialVolume(volume);
    reactor.initialize();
    reactor.setEnergy(energy);
    // Setup network
    ReactorNet network;
    network.addReactor(reactor);
    network.initialize();
    // Internal preconditioner
    AdaptivePreconditioner internalPrecon;
    preconditionerInitialize(network, internalPrecon);
    // State produced within CVODES for this example
    vector_fp N(reactor.neq(), 0.0);
    vector_fp Ndot(reactor.neq(), 0.0);
    vector_fp NCopy(reactor.neq(), 0.0);
    // Creating external preconditioner for comparison
    AdaptivePreconditioner externalPrecon;
    preconditionerInitialize(network, externalPrecon);
    // Compare the preconditioners while stepping
    double ct = 0.0;
    for (size_t i = 0; i < 1000; i++)
    {
        reactor.getState(N.data());
        // precondition from reactor object
        preconditionWithClass(internalPrecon, reactor, N.data(), Ndot.data(), ct);

        // manual precondition
        // set state to strictly positive composition for manual comparison
        externalPrecon.getStrictlyPositiveComposition(reactor.neq(), N.data(), NCopy.data());
        reactor.updateState(NCopy.data());
        // reset former values
        externalPrecon.reset();
        twoStepManualPrecondition(externalPrecon, reactor, NCopy.data(), Ndot.data(), ct);
        // post setup processes
        externalPrecon.setup();
        // check that the two are equal
        EXPECT_EQ(externalPrecon == internalPrecon, true);
        // reset state
        reactor.updateState(N.data());
        // step the network
        ct = network.step();
    }
}

TEST(AdaptivePreconditionerTestSet, test_two_step_mechanism)
{
    // testing with energy off
    TwoStepMechanism(false);
    // testing with energy on
    TwoStepMechanism(true);
}

TEST(AdaptivePreconditionerTestSet, test_two_step_multitype_network)
{
    // Constants
    double volume = 0.2;
    // Setting up solution object and thermo/kinetics pointers
    auto sol = newSolution("methane_twostep.yaml");
    auto thermo = sol->thermo();
    auto kinetics = sol->kinetics();
    thermo->setState_TPX(1000.0, 101325, "CH4:1.0, O2:1.0, H2O:0, CO2:0, CO:0");
    // Set up reactor object
    IdealGasConstPressureMoleReactor reactor1;
    IdealGasMoleReactor reactor2;
    reactor1.insert(sol);
    reactor1.setInitialVolume(volume);
    reactor2.insert(sol);
    reactor2.setInitialVolume(volume);
    // Setup network for initialization
    ReactorNet network;
    network.addReactor(reactor1);
    network.addReactor(reactor2);
    // Internal preconditioner
    AdaptivePreconditioner internalPrecon;
    network.setProblemType(GMRES);
    network.setPreconditioner(internalPrecon);
    network.initialize();
    // Creating external preconditioner for comparison
    AdaptivePreconditioner externalPrecon;
    preconditionerInitialize(network, externalPrecon);
    // State produced within CVODES for this example
    vector_fp N(network.neq(), 0.0);
    vector_fp Ndot(network.neq(), 0.0);
    vector_fp NCopy(network.neq(), 0.0);
    network.getState(N.data());
    // Get used gamma value
    double gamma = internalPrecon.getGamma();
    // Compare the preconditioners while stepping
    double ct = 0.0;
    for (size_t i = 0; i < 1000; i++)
    {
        // precondition from reactor object
        network.preconditionerSetup(ct, N.data(), Ndot.data(), nullptr, gamma);
        externalPrecon.getStrictlyPositiveComposition(network.neq(), N.data(), NCopy.data());
        network.updateState(NCopy.data());
        // manual precondition external
        externalPrecon.reset();
        externalPrecon.m_ctr = 0;
        twoStepManualPrecondition(externalPrecon, reactor1, NCopy.data(), Ndot.data(), ct);
        externalPrecon.m_ctr = 1;
        twoStepManualPrecondition(externalPrecon, reactor2, NCopy.data(), Ndot.data(), ct);
        externalPrecon.transformJacobianToPreconditioner();
        network.updateState(N.data());
        // Check that the two are equal
        EXPECT_EQ(externalPrecon == internalPrecon, true);
        // step the network
        ct = network.step();
    }
}

TEST(AdaptivePreconditionerTestSet, test_one_step_mechanism_network)
{
    // Constants
    double volume = 0.2;
    double volInv = 1/volume;
    double startTime = 0.0;
    double sharedThreshold = 1e-16;
    double gamma = 1.0;
    size_t numberOfReactors = 3;
    bool energy = true;
    // Network
    ReactorNet network;
    IdealGasConstPressureMoleReactor reactor;
    auto sol = newSolution("methane_onestep.yaml");
    auto thermo = sol->thermo();
    auto kinetics = sol->kinetics();
    thermo->setState_TPX(300.0, 101325, "CH4:1.0, O2:1.0, H2O:0, CO2:0");
    // Set up reactor object
    reactor.insert(sol);
    reactor.setInitialVolume(volume);
    // Setting up solution object and thermo/kinetics pointers one
    for (size_t i = 0; i < numberOfReactors; i++)
    {
        network.addReactor(reactor);
    }
    reactor.setEnergy(energy);
    // Create and add preconditioner
    AdaptivePreconditioner internalPrecon;
    internalPrecon.setThreshold(sharedThreshold);
    network.setProblemType(GMRES);
    network.setPreconditioner(internalPrecon);
    network.initialize();
    // State produced within CVODES for this example
    vector_fp Ndot(network.neq(), 0.0);
    vector_fp N(network.neq(), 0.0);
    vector_fp rhs(network.neq(), 0.0);
    vector_fp output(network.neq(), 0.0);
    network.getState(N.data());
    network.getState(rhs.data());
    // Creating external preconditioner for comparison
    AdaptivePreconditioner externalPrecon;
    externalPrecon = internalPrecon;
    // preconditioner internal preconditioner
    int flag = network.FuncEval::preconditioner_setup_nothrow(startTime, N.data(), Ndot.data(), gamma);
    // precondition external preconditioner
    externalPrecon.getStrictlyPositiveComposition(reactor.neq(), N.data(), N.data());
    reactor.updateState(N.data());
    // Getting data for manual calculations
    vector_fp kf(kinetics->nReactions(), 0.0);
    kinetics->getFwdRateConstants(kf.data());
    kinetics->thirdbodyConcMultiply(kf.data());
    vector_fp C(thermo->nSpecies(), 0.0);
    thermo->getConcentrations(C.data());
    // Assign concentrations to species
    double O2 = C[0];
    double CH4 = C[2];
    // Setting elements
    // O2
    externalPrecon(1, 1) = -4 * kf[0] * O2 * CH4 * volInv; // dO2/dO2
    externalPrecon(1, 3) = -2 * kf[0] * O2 * O2 * volInv; // dO2/dCH4
    // H2O
    externalPrecon(2, 1) = 4 * kf[0] * O2 * CH4 * volInv; // dH2O/dO2
    externalPrecon(2, 3) = 2 * kf[0] * O2 * O2 * volInv; // dH2O/dCH4
    // CH4
    externalPrecon(3, 1) = -2 * kf[0] * O2 * CH4 * volInv; // dCH4/dO2
    externalPrecon(3, 3) = -kf[0] * O2 * O2 * volInv; // dCH4/dCH4
    // CO2
    externalPrecon(4, 1) = 2 * kf[0] * O2 * CH4 * volInv; // dCO2/dO2
    externalPrecon(4, 3) = kf[0] * O2 * O2 * volInv; // dCO2/dCO2
    if (reactor.energyEnabled())
    {
        reactor.StateDerivatives(externalPrecon, startTime, N.data(), Ndot.data(), nullptr);
    }
    size_t neq = reactor.neq();
    for (size_t k = 1; k < numberOfReactors; k++)
    {
        for (size_t i = 0; i < neq; i++)
        {
            for (size_t j = 0; j < neq; j++)
            {
                externalPrecon(i + neq * k, j + neq * k) = externalPrecon(i, j);
            }
        }
    }
    // Make into preconditioner as P = (I - gamma * J_bar)
    externalPrecon.transformJacobianToPreconditioner();
    // Check that the two are equal
    EXPECT_EQ(externalPrecon==internalPrecon, true);
    internalPrecon.getMatrix()->setZero();
    network.preconditionerSetup(startTime, N.data(), Ndot.data(), nullptr, 1.0);
    EXPECT_EQ(externalPrecon==internalPrecon, true);
    EXPECT_EQ(flag, 0);
    EXPECT_EQ(network.FuncEval::preconditioner_solve_nothrow(startTime, N.data(), Ndot.data(), rhs.data(), output.data()), 0);
    EXPECT_NO_THROW(network.preconditionerSolve(startTime, N.data(), Ndot.data(), rhs.data(), output.data()));
    EXPECT_EQ(externalPrecon.getGamma(), gamma);
    // test copy
    externalPrecon.reset();
    externalPrecon = internalPrecon;
    EXPECT_EQ(externalPrecon==internalPrecon, true);
    // Reset preconditioner then compare again
    internalPrecon.reset();
    EXPECT_EQ(externalPrecon==internalPrecon, false);
    // Call assignment then compare again
    internalPrecon = externalPrecon;
    EXPECT_EQ(externalPrecon==internalPrecon, true);
}

TEST(AdaptivePreconditionerTestSet, test_run_sim)
{
    // Setting up solution object and thermo/kinetics pointers
    auto sol = newSolution("methane_twostep.yaml");
    sol->thermo()->setState_TPX(300.0, 101325, "CH4:1.0, O2:1.0, H2O:0, CO2:0");
    // Set up reactor object
    IdealGasConstPressureMoleReactor reactor;
    reactor.insert(sol);
    reactor.setInitialVolume(1.0);
    // Creating inlet reservoir object and adding gas
    Reservoir inlet;
    inlet.insert(sol);
    //Creating exhaust reservoir object and adding gas
    Reservoir exhaust;
    exhaust.insert(sol);
    // Creating mass flow controllers
    MassFlowController inletMassFlowController;
    PressureController outletMassFlowController;
    //Connecting reactors
    inletMassFlowController.install(inlet,reactor);
    outletMassFlowController.install(reactor,exhaust);
    outletMassFlowController.setMaster(&inletMassFlowController);
    outletMassFlowController.setPressureCoeff(0.01);
    // Set constant massflow rate
    inletMassFlowController.setMassFlowRate(1.0);
    // Creating reactor network
    ReactorNet network;
    // Create and add preconditioner
    AdaptivePreconditioner precon;
    network.addReactor(reactor); //Adding reactor to network
    network.setProblemType(GMRES);
    network.setPreconditioner(precon);
    // Setting up simulation
    network.setInitialTime(0.0);
    network.setMaxTimeStep(0.1);
    network.setMaxSteps(10000);
    network.setTolerances(1e-6,1e-6);
    network.setSensitivityTolerances(1e-6,1e-6);
    network.step();
}

TEST(AdaptivePreconditionerTestSet, test_preconditioned_hydrogen_auto_ignition)
{
    // create an ideal gas mixture that corresponds to GRI-Mech 3.0
    auto sol = newSolution("gri30.yaml", "gri30", "None");
    auto gas = sol->thermo();
    // set the state
    gas->setState_TPX(1001.0, OneAtm, "H2:2.0, O2:1.0, N2:4.0");
    // create a reactor
    IdealGasConstPressureMoleReactor r;
    // 'insert' the gas into the reactor and environment.
    r.insert(sol);
    // create preconditioner
    AdaptivePreconditioner precon;
    // create reactor network and set to use preconditioner
    ReactorNet sim;
    sim.addReactor(r);
    sim.setProblemType(GMRES);
    sim.setPreconditioner(precon);
    // main loop
    double dt = 1.e-5; // interval at which output is written
    int nsteps = 100; // number of intervals
    for (int i = 1; i <= nsteps; i++) {
        double tm = i*dt;
        sim.advance(tm);
    }
}

TEST(AdaptivePreconditionerTestSet, test_utilities_get_set)
{
    // Seed random number generator
    std::srand(std::time(NULL));
    // Generate random test dimensions between 10 and 50
    size_t base = 10;
    size_t limit = 50;
    size_t randomNum = std::rand() % (limit-base+1) + base;
    std::vector<size_t> dims{randomNum, randomNum};
    // Create preconditioner object
    AdaptivePreconditioner precon;
    auto precon_mat = precon.getMatrix();
    precon_mat->resize(dims[0], dims[1]);
    // Set threshold
    double thresh = (double) base+2;
    precon.setThreshold(thresh);
    // Set dimensions randomly
    precon.setDimensions(&dims);
    // Get dimensions newly
    std::vector<size_t> *preconDims = precon.getDimensions();
    // Test that the dimensions are set properly via setDimensions
    EXPECT_EQ (dims.at(0), preconDims->at(0));
    EXPECT_EQ (dims.at(1), preconDims->at(1));
    // Testing some get/set utilities
    double randSetGet = std::rand();
    precon.setThreshold(randSetGet);
    EXPECT_EQ(precon.getThreshold(), randSetGet);
    precon.setAbsoluteTolerance(randSetGet);
    EXPECT_EQ(precon.getAbsoluteTolerance(), randSetGet);
}

TEST(PreconditionerBaseTestSet, test_preconditioner_base)
{
    // Required variables
    PreconditionerBase precon;
    IdealGasConstPressureMoleReactor pressureReactor;
    Reactor generalReactor;
    ReactorNet network;
    // Test error throwing of preconditioner base
    EXPECT_THROW(precon.solve(0, nullptr, nullptr), NotImplementedError);
    EXPECT_THROW(precon.setup(), NotImplementedError);
    EXPECT_EQ(PRECONDITIONER_NOT_SET, precon.getPreconditionerMethod());
}

TEST(PreconditionerBaseTestSet, test_funceval_precon_funcs)
{
    // Setting up solution to insert in reactor
    auto sol = newSolution("h2o2.yaml");
    // Set up reactor object
    Reactor reactor;
    reactor.insert(sol);
    ReactorNet network;
    network.addReactor(reactor);
    network.initialize();
    EXPECT_THROW(network.FuncEval::preconditionerSetup(0, nullptr, nullptr, nullptr, 0.0), NotImplementedError);
    EXPECT_THROW(network.FuncEval::preconditionerSolve(0, nullptr, nullptr, nullptr, nullptr), NotImplementedError);
}

int main(int argc, char** argv)
{
    printf("Running main() from test_preconditioners.cpp\n");
    testing::InitGoogleTest(&argc, argv);
    Cantera::make_deprecation_warnings_fatal();
    int result = RUN_ALL_TESTS();
    Cantera::appdelete();
    return result;
}
