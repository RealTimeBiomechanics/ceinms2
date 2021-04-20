#include "ceinms2/Lloyd2003Muscle.h"
#include "ceinms2/DataTable.h"
#include "ceinms2/RienerActivation.h"
#include <iostream>
#include <vector>
#include <string>
#include <filesystem>
#include <limits>
#include <chrono>
using namespace ceinms;
using std::string;

int test1() {

    const CubicSpline act{ ceinms::getDefaultActiveForceLengthCurve() };
    const CubicSpline pas(ceinms::getDefaultPassiveForceLengthCurve());
    const CubicSpline vel(ceinms::getDefaultForceVelocityCurve());
    const CubicSpline ten(ceinms::getDefaultTendonForceStrainCurve());

    Lloyd2003Muscle::Parameters p;
    p.damping = 0.1;
    p.maxContractionVelocity = 10;
    p.maxIsometricForce = 6194;
    p.optimalFiberLength = 0.04;
    p.pennationAngleAtOptimalFiberLength = 0.38;
    p.percentageChange = 0.15;
    p.strengthCoefficient = 1;
    p.tendonSlackLength = 0.259;
    p.activeForceLengthCurve = act;
    p.passiveForceLengthCurve = pas;
    p.forceVelocityCurve = vel;
    p.tendonForceStrainCurve = ten;

	Lloyd2003Muscle mtu(p);
    std::filesystem::path root(ROOT_DIR);
    std::filesystem::path filename{ root
                                    / std::filesystem::path("data/gait/soleus.csv") };
    DataTable<double> data = fillFromCSV<double>(filename.string());
    vector<double> activation = data.getColumn("activation");
    vector<double> mtuLength = data.getColumn("mtuLength");
    vector<double> mtuStiffness, fiberStiffness, tendonStiffness;
    mtu.equilibrate();
    for (int i{ 0 }; i < activation.size(); ++i) {
        mtu.setActivation(activation.at(i));
        mtu.setMusculotendonLength(mtuLength.at(i));
        mtu.evaluate(0.005);
        mtu.calculateOutput();
        mtuStiffness.push_back(mtu.getOutput().musculotendonStiffness);
        fiberStiffness.push_back(mtu.getOutput().fiberStiffness);
        tendonStiffness.push_back(mtu.getOutput().tendonStiffness);
    }
    data.pushColumn("mtuStiffness", mtuStiffness);
    data.pushColumn("fiberStiffness", fiberStiffness);
    data.pushColumn("tendonStiffness", tendonStiffness);

    data.print("outStiffness.csv");

    double baseVel = mtu.getState().fiberVelocity;
    mtu.integrate(0.1);
    mtu.validateState();
    mtu.getState().fiberVelocity = 0;
    mtu.calculateOutput();
    return 0;
}

auto calculateJacobianWithCentralDifferences(Lloyd2003Muscle mtu) {
    auto param = mtu.getParameters();
    auto input = mtu.getInput();
    std::vector<DoubleT*> w(10);
    w[0] = &input.musculotendonLength;
    w[1] = &input.activation;
    w[2] = &param.damping;
    w[3] = &param.maxContractionVelocity;
    w[4] = &param.maxIsometricForce;
    w[5] = &param.optimalFiberLength;
    w[6] = &param.pennationAngleAtOptimalFiberLength;
    w[7] = &param.percentageChange;
    w[8] = &param.strengthCoefficient;
    w[9] = &param.tendonSlackLength;
  
    auto calculateForce{
        [](Lloyd2003Muscle &mtu, Lloyd2003Muscle::Input i, Lloyd2003Muscle::Parameters p) {
            mtu.setInput(i);
            mtu.setParameters(p);
            mtu.integrate(0.005);
            return mtu.calculateFiberForceProjectedOnTendon();
        }
    };
    Eigen::Matrix<DoubleT, 1, 10> J;
    DoubleT eps = std::numeric_limits<DoubleT>::epsilon()*10.;
    for (int i{ 0 }; i < w.size(); ++i) {
        DoubleT center = *w[i];
        *w[i] = center + eps;
        DoubleT fUp = calculateForce(mtu, input, param);
        *w[i] = center - eps;
        DoubleT fDown = calculateForce(mtu, input, param);
        J(i) = (fUp - fDown) / (2 * eps);
    }
    return J;
}


int test2() {
    const CubicSpline act{ ceinms::getDefaultActiveForceLengthCurve() };
    const CubicSpline pas(ceinms::getDefaultPassiveForceLengthCurve());
    const CubicSpline vel(ceinms::getDefaultForceVelocityCurve());
    const CubicSpline ten(ceinms::getDefaultTendonForceStrainCurve());

    Lloyd2003Muscle::Parameters p;
    p.damping = 0.1;
    p.maxContractionVelocity = 10;
    p.maxIsometricForce = 6194;
    p.optimalFiberLength = 0.04;
    p.pennationAngleAtOptimalFiberLength = 0.38;
    p.percentageChange = 0.15;
    p.strengthCoefficient = 1;
    p.tendonSlackLength = 0.259;
    p.activeForceLengthCurve = act;
    p.passiveForceLengthCurve = pas;
    p.forceVelocityCurve = vel;
    p.tendonForceStrainCurve = ten;

    Lloyd2003Muscle mtu(p);
    std::filesystem::path root(ROOT_DIR);
    std::filesystem::path filename{ root / std::filesystem::path("data/gait/soleus.csv") };
    DataTable<double> data = fillFromCSV<double>(filename.string());
    vector<double> activation = data.getColumn("activation");
    vector<double> mtuLength = data.getColumn("mtuLength");
    vector<double> mtuStiffness, fiberStiffness, tendonStiffness;
    vector<Eigen::Matrix<DoubleT, 1, 10>> Jautos, Jnumerics;

    mtu.equilibrate();
    for (int i{ 0 }; i < activation.size(); ++i) {
        mtu.setActivation(activation.at(i));
        mtu.setMusculotendonLength(mtuLength.at(i));
        mtu.evaluate(0.005);
        mtu.calculateOutput();
        auto J = mtu.calculateJacobian();
        Jautos.push_back(J);
        auto J1 = calculateJacobianWithCentralDifferences(mtu);
        Jnumerics.push_back(J1);
    }
    
    auto print { [](vector<Eigen::Matrix<DoubleT, 1, 10>> &jac, std::string filename) {
        std::ofstream outF(filename);
        for (auto &j : jac) {
            for (int i{ 0 }; i < j.size(); ++i)
                outF << j(i) << ", ";
            outF << '\n';
        }
    }};

    print(Jnumerics, "numericJ.csv");
    print(Jautos, "autoJ.csv");

    return 0;
}

int main() {
    test2();

}

