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

void twoStepManualPrecondition(AdaptivePreconditioner* externalPrecon, MoleReactor* reactor, double* N, double* Ndot, double t=0.0, size_t loc=0)
    {
        auto thermo = reactor->getThermoMgr();
        auto kinetics = reactor->getKineticsMgr();
        double volInv = 1/reactor->volume();
        // Getting data for manual calculations
        size_t numberOfReactions = kinetics->nReactions();
        // vectors for manual calcs
        std::vector<double> kf(numberOfReactions, 0.0);
        std::vector<double> kr(numberOfReactions, 0.0);
        // getting actual data
        std::vector<double> C(thermo->nSpecies(), 0.0);
        thermo->getConcentrations(C.data());
        kinetics->getFwdRateConstants(kf.data());
        kinetics->getRevRateConstants(kr.data());
        // Assign concentrations to species
        double O2 = C[0];
        double CH4 = C[2];
        double CO = C[4];
        // Setting elements
        // O2
        externalPrecon->operator()(1 + loc, 1 + loc) = (- 2.25 * kf[0] * CH4 * std::pow(O2, 0.5) - 0.25 * kf[1] * CO * std::pow(O2, -0.5)) * volInv; // dO2/dO2
        externalPrecon->operator()(1 + loc, 3 + loc) = - 1.5 * kf[0] * std::pow(O2, 1.5) * volInv; // dO2/dCH4
        externalPrecon->operator()(1 + loc, 4 + loc) = 0.5 * kr[1] * volInv; // dO2/dCO2
        externalPrecon->operator()(1 + loc, 5 + loc) = - 0.5 * kf[1] * std::pow(O2, 0.5) * volInv; // dO2/dCO
        // H2O
        externalPrecon->operator()(2 + loc, 1 + loc) = 3 * kf[0] * CH4 * std::pow(O2, 0.5) * volInv; // dH2O/dO2
        externalPrecon->operator()(2 + loc, 3 + loc) = 2 * kf[0] * std::pow(O2, 1.5) * volInv; // dH2O/dCH4
        // CH4
        externalPrecon->operator()(3 + loc, 1 + loc) = - 1.5 * kf[0] * CH4 * std::pow(O2, 0.5) * volInv; // dCH4/dO2
        externalPrecon->operator()(3 + loc, 3 + loc) = - kf[0] * std::pow(O2, 1.5) * volInv; // dCH4/dCH4
        // CO2
        externalPrecon->operator()(4 + loc, 1 + loc) = (0.5 * kf[1] * CO * std::pow(O2, -0.5)) * volInv; // dCO2/dO2
        externalPrecon->operator()(4 + loc, 4 + loc) = -kr[1] * volInv; // dCO2/dCO2
        externalPrecon->operator()(4 + loc, 5 + loc) = kf[1] * std::pow(O2, 0.5) * volInv; // dCO2/CO
        //CO
        externalPrecon->operator()(5 + loc, 1 + loc) = (1.5 * kf[0] * CH4 * std::pow(O2, 0.5) - 0.5 * kf[1] * CO * std::pow(O2, -0.5)) * volInv; // dCO/dO2
        externalPrecon->operator()(5 + loc, 3 + loc) = kf[0] * std::pow(O2, 1.5) * volInv; // dCO/dCH4
        externalPrecon->operator()(5 + loc, 4 + loc) = kr[1] * volInv; // dCO/dCO2
        externalPrecon->operator()(5 + loc, 5 + loc) = -kf[1] * std::pow(O2,0.5) * volInv; // dCO/CO
        // temperature derivatives
        if (reactor->energyEnabled())
        {
            externalPrecon->TemperatureDerivatives(reactor, t, N, Ndot, nullptr);
        }
    };

void preconditionWithClass(AdaptivePreconditioner* internalPrecon, IdealGasConstPressureMoleReactor* reactor, double* N, double* Ndot, double t=0.0)
{
    // class based
    reactor->getState(N);
    internalPrecon->reset();
    internalPrecon->reactorLevelSetup(reactor, t, N, Ndot, nullptr);
    internalPrecon->transformJacobianToPreconditioner();
}

void preconditionerInitialize(ReactorNet* network, AdaptivePreconditioner* precon, size_t reactorStart=0, double gamma=1.0, double thresh=0)
{
    precon->initialize(network);
    precon->setReactorStart(reactorStart);
    precon->setGamma(gamma);
    precon->setThreshold(thresh);
}

void TwoStepMechanism(bool energy)
{
    // Constants
    double volume = 0.2;
    // Setting up solution object and thermo/kinetics pointers
    auto sol = newSolution("methane_twostep.yaml");
    auto thermo = sol->thermo();
    auto kinetics = sol->kinetics();
    thermo->setState_TPX(1000.0, 101325, "CH4:1.0, O2:1.0, H2O:0, CO2:0, CO:0");
    // Set up reactor object
    IdealGasConstPressureMoleReactor reactor;
    reactor.insert(sol);
    reactor.setInitialVolume(volume);
    reactor.initialize();
    reactor.setEnergy(energy);
    // Setup network for initialization
    ReactorNet network;
    network.addReactor(reactor);
    network.initialize();
    // State produced within CVODES for this example
    std::vector<double> N(reactor.neq(), 0.0);
    std::vector<double> Ndot(reactor.neq(), 0.0);
    std::vector<double> NCopy(reactor.neq(), 0.0);
    reactor.getState(N.data());
    // Internal preconditioner
    AdaptivePreconditioner internalPrecon;
    preconditionerInitialize(&network, &internalPrecon);
    // Creating external preconditioner for comparison
    AdaptivePreconditioner externalPrecon;
    preconditionerInitialize(&network, &externalPrecon);
    // Compare the preconditioners while stepping
    double ct = 0.0;
    for (size_t i = 0; i < 1000; i++)
    {
        // precondition from reactor object
        preconditionWithClass(&internalPrecon, &reactor, N.data(), Ndot.data(), ct);
        externalPrecon.getStrictlyPositiveComposition(reactor.neq(), N.data(), NCopy.data());
        reactor.updateState(NCopy.data());
        // manual precondition
        externalPrecon.reset();
        twoStepManualPrecondition(&externalPrecon, &reactor, NCopy.data(), Ndot.data(), ct);
        externalPrecon.transformJacobianToPreconditioner();
        reactor.updateState(N.data());
        // Check that the two are equal
        EXPECT_EQ(externalPrecon == internalPrecon, true);
        // step the network
        ct = network.step();
    }
    // Reset Internal and test acceptPreconditioner
    internalPrecon.reset();
    reactor.acceptPreconditioner(&internalPrecon, ct, N.data(), Ndot.data(), nullptr);
    internalPrecon.transformJacobianToPreconditioner();
    EXPECT_EQ(externalPrecon == internalPrecon, true);
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
    network.setIntegratorType(&internalPrecon, GMRES);
    network.initialize();
    // Creating external preconditioner for comparison
    AdaptivePreconditioner externalPrecon;
    preconditionerInitialize(&network, &externalPrecon);
    // State produced within CVODES for this example
    std::vector<double> N(network.neq(), 0.0);
    std::vector<double> Ndot(network.neq(), 0.0);
    std::vector<double> NCopy(network.neq(), 0.0);
    network.getState(N.data());
    // Get used gamma value
    double gamma = internalPrecon.getGamma();
    // Compare the preconditioners while stepping
    double ct = 0.0;
    for (size_t i = 0; i < 1000; i++)
    {
        // precondition from reactor object
        network.preconditionerSetup(ct, N.data(), Ndot.data(), gamma);
        externalPrecon.getStrictlyPositiveComposition(network.neq(), N.data(), NCopy.data());
        network.updateState(NCopy.data());
        // manual precondition external
        externalPrecon.reset();
        externalPrecon.setReactorStart(0);
        twoStepManualPrecondition(&externalPrecon, &reactor1, NCopy.data(), Ndot.data(), ct);
        externalPrecon.setReactorStart(reactor1.neq());
        twoStepManualPrecondition(&externalPrecon, &reactor2, NCopy.data(), Ndot.data(), ct);
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
    network.setIntegratorType(&internalPrecon, GMRES);
    network.initialize();
    // State produced within CVODES for this example
    std::vector<double> Ndot(network.neq(), 0.0);
    std::vector<double> N(network.neq(), 0.0);
    std::vector<double> rhs(network.neq(), 0.0);
    std::vector<double> output(network.neq(), 0.0);
    network.getState(N.data());
    network.getState(rhs.data());
    int flag = network.FuncEval::preconditioner_setup_nothrow(startTime, N.data(), Ndot.data(), gamma);
    // Creating external preconditioner for comparison
    AdaptivePreconditioner externalPrecon;
    externalPrecon.initialize(&network);
    externalPrecon.setThreshold(sharedThreshold);
    externalPrecon.setReactorStart(0);
    externalPrecon.setGamma(gamma);
    externalPrecon.getStrictlyPositiveComposition(reactor.neq(), N.data(), N.data());
    reactor.updateState(N.data());
    // Getting data for manual calculations
    std::vector<double> kf(kinetics->nReactions(), 0.0);
    kinetics->getFwdRateConstants(kf.data());
    kinetics->thirdbodyConcMultiply(kf.data());
    std::vector<double> C(thermo->nSpecies(), 0.0);
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
        externalPrecon.TemperatureDerivatives(&reactor, startTime, N.data(), Ndot.data(), nullptr);
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
    network.preconditionerSetup(startTime, N.data(), Ndot.data(), 1.0);
    EXPECT_EQ(externalPrecon==internalPrecon, true);
    EXPECT_EQ(flag, 0);
    EXPECT_EQ(network.FuncEval::preconditioner_solve_nothrow(startTime, N.data(), Ndot.data(), rhs.data(), output.data()), 0);
    EXPECT_NO_THROW(network.preconditionerSolve(startTime, N.data(), Ndot.data(), rhs.data(), output.data()));
    EXPECT_EQ(externalPrecon.getGamma(), gamma);
    // test copy constructor
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
    network.setIntegratorType(&precon,GMRES);
    // network->setVerbose(); //Setting verbose to be true
    network.addReactor(reactor); //Adding reactor to network
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
    sim.setIntegratorType(&precon,GMRES);
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
    // Ensure all elements are zero
    precon_mat->setZero();
    // Set threshold
    double thresh = (double) base+2;
    precon.setThreshold(thresh);
    // Random values to put in matrix
    std::queue<size_t> values;
    for(size_t i=0; i<base; i++)
    {
        values.push(std::rand() % 100);
        values.push(std::rand() % dims[1]);
        values.push(std::rand() % dims[0]);
    }
    // Check set and get elements
    for(size_t i=0; i<base; i++)
    {
        double currElement = (double) values.front();
        values.pop();
        size_t col = values.front();
        values.pop();
        size_t row = values.front();
        values.pop();
        precon.setElement(row, col, currElement);
        double returnedElement = precon.getElement(row, col);
        if(std::abs(currElement) >= thresh || row == col)
        {
            EXPECT_DOUBLE_EQ (returnedElement, currElement);
        }
        else
        {
            EXPECT_DOUBLE_EQ (returnedElement, 0.0);
        }
    }
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
    precon.setReactorStart(randSetGet);
    EXPECT_EQ(precon.getReactorStart(), randSetGet);
    precon.setAbsoluteTolerance(randSetGet);
    EXPECT_EQ(precon.getAbsoluteTolerance(), randSetGet);
}

TEST(PreconditionerBaseTestSet, test_preconditioner_base)
{
    // Required variables
    PreconditionerBase precon;
    IdealGasConstPressureMoleReactor pressureReactor;
    Reactor generalReactor;
    // Test error throwing of preconditioner base
    EXPECT_THROW(precon.solve(nullptr, nullptr, nullptr), NotImplementedError);
    EXPECT_THROW(precon.setup(nullptr, 0.0, nullptr, nullptr, 1), NotImplementedError);
    EXPECT_THROW(generalReactor.acceptPreconditioner(&precon, 0.0, nullptr, nullptr, nullptr), NotImplementedError);
    EXPECT_THROW(pressureReactor.acceptPreconditioner(&precon, 0.0, nullptr, nullptr, nullptr), NotImplementedError);
    EXPECT_EQ(PRECONDITIONER_NOT_SET, precon.getPreconditionerType());
    EXPECT_THROW(precon.initialize(nullptr), NotImplementedError);
    EXPECT_THROW(precon.setElement(0, 0, 0.0), NotImplementedError);
    EXPECT_THROW(precon.getElement(0, 0), NotImplementedError);
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
    EXPECT_THROW(network.FuncEval::preconditionerSetup(0, nullptr, nullptr, 0.0), NotImplementedError);
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
