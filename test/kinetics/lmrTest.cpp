// #include "ext/googletest/googletest/include/gtest/gtest.h"
#include "gtest/gtest.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/base/Solution.h"
#include "cantera/base/global.h"

// COMPILE/EXECUTE COMMANDS
// gcc test/kinetics/lmrTest.cpp -o test/kinetics/lmrTest -lstdc++
// lmrTest.exe

namespace Cantera
{

class lmrTest : public testing::Test
{
public:
    lmrTest() {}

    static void SetUpTestCase() {
        soln_ = newSolution("../data/kineticsfromscratch_LMRtest.yaml"); //import kineticsfromscratch-lmrtest.yaml
    }

    static void TearDownTestCase() {
        soln_.reset();
    }

protected:
    static shared_ptr<Solution> soln_;

    // void set_TP(double T, double P) { //will likely need to set T, P, X (need to set whatever params are needed to initialize the state)
    //     T_ = T;                         //maybe put the T,P,X defn in the TEST_F directly
    //     RT_ = GasConst_cal_mol_K * T;
    //     P_ = P;
    //     soln_->thermo()->setState_TP(T_, P_);
    // }

    // double k(double A, double n, double Ea) {
    //     return A * pow(T_, n) * exp(-Ea/RT_);
    // }
    // double T_, RT_, P_;
};

shared_ptr<Solution> lmrTest::soln_;

TEST_F(lmrTest, reactionCounts) 
{
    EXPECT_EQ((size_t) 14, soln_->kinetics()->nReactions());
}

TEST_F(lmrTest, PlogLowPressure)
{

    string X = "H:1.0, O2:0.0";
    soln_->thermo()->setState_TPX(900.0, 101325 * 8.0, X);

    double ktest1=10;
    double ktest2=10;


    // Test that P-log reactions have the right low-pressure limit
    //set_TP(500.0, 1e-7);
    //vector<double> kf(7); //7 is number of rxns in that input file. But the solution obj can be used to automate this (the test above chaks that numb rxns is 7)
    // soln_->kinetics()->getFwdRateConstants(&kf[0]);

    // Pre-exponential factor decreases by 10^3 for second-order reaction
    // when converting from cm + mol to m + kmol
    // double kf0 = k(1.212400e+13, -0.5779, 10872.7);

    // EXPECT_NEAR(kf0, kf[0], 1e-9 * kf0); //third arg is the relative tolerance (we have floating point numerals which are inexact)
    EXPECT_DOUBLE_EQ(ktest1, ktest2);
    //EXPECT_NEAR(kf3, kf[3], 1e-9 * kf3); //kf[3], e.g. selects rxn 3 from input file
}
} // namespace Cantera

int main(int argc, char** argv)
{
    printf("Running main() from lmrTest.cpp\n");
    Cantera::make_deprecation_warnings_fatal();
    Cantera::printStackTraceOnSegfault();
    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
    Cantera::appdelete();
    return result;
}
