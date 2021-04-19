#include "ceinms2/Lloyd2003Muscle.h"
#include "ceinms2/DataTable.h"
#include "ceinms2/RienerActivation.h"
#include <iostream>
#include <vector>
#include <string>
#include <filesystem>
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
    mtu.equilibrate();
    for (int i{ 0 }; i < activation.size(); ++i) {
        mtu.setActivation(activation.at(i));
        mtu.setMusculotendonLength(mtuLength.at(i));
        mtu.evaluate(0.005);
        mtu.calculateOutput();
        auto J = mtu.calculateJacobian();
        std::cout << J << std::endl;
    }
}

int main() {
    test2();

}

