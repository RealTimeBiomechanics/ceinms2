#include "ceinms2/Lloyd2003Muscle.h"
#include "ceinms2/DataTable.h"
//#include "ceinms2/RienerActivation.h"
#include <iostream>
#include <vector>
#include <string>
#include <filesystem>
#include <limits>
#include <chrono>
#include <algorithm>
using namespace ceinms;
using std::string;
using std::back_inserter;

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



class Adam {
  public:
    Adam(unsigned N)
        : theta(N, 1)
        , v(N, 1)
        , m(N, 1) {}
    Eigen::MatrixXd step(Eigen::MatrixXd &gradients);

    double alpha = 0.01;// learning rate
    double beta1 = 0.9;// decay1
    double beta2 = 0.99;// decay2
    double epsilon = 1e-8;
    Eigen::MatrixXd theta;// parameters
    Eigen::MatrixXd v;// exponential moving average
    Eigen::MatrixXd m;// exponential moving average
    int t = 0;
};

Eigen::MatrixXd Adam::step(Eigen::MatrixXd &gradients) {
    auto g = gradients.array() * gradients.array();
    t += 1;
    m = beta1 * m + (1 - beta1) * gradients;
    v = beta2 * v + (1 - beta2) * Eigen::MatrixXd(gradients.array() * gradients.array());
    Eigen::MatrixXd m_hat = m.array() / (1 - std::pow(beta1, t));
    Eigen::MatrixXd v_hat = v.array() / (1 - std::pow(beta2, t));
    theta = theta.array() - alpha * (m_hat.array() / (v_hat.array().sqrt() - epsilon).array());
    return theta;
}

int test3() {
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

    auto runMuscle{ [&activation, &mtuLength, &mtu](Lloyd2003Muscle::Parameters &par) {
        vector<double> forces;
        mtu.setParameters(par);
        mtu.equilibrate();
        for (int i{ 0 }; i < activation.size(); ++i) {
            mtu.setActivation(activation.at(i));
            mtu.setMusculotendonLength(mtuLength.at(i));
            mtu.evaluate(0.005);
            mtu.calculateOutput();
            forces.push_back(mtu.getOutput().force);
        }
        return forces;
    } };

    auto reference = runMuscle(p);

    auto J{ [&](const Eigen::MatrixXd &theta) {
        Lloyd2003Muscle::Parameters par = p;
   /*     par.damping = theta(0);
        par.maxContractionVelocity = theta(1);
        par.maxIsometricForce = theta(2);
        par.optimalFiberLength = theta(3);
        par.pennationAngleAtOptimalFiberLength = theta(4);
        par.percentageChange = theta(5);
        par.strengthCoefficient = theta(6);
        par.tendonSlackLength = theta(7);*/
        par.strengthCoefficient = theta(0);
        mtu.setParameters(par);
        std::vector<double> out;
        auto estimates = runMuscle(par);
        std::transform(
            reference.cbegin(),
            reference.cend(),
            estimates.cbegin(),
            back_inserter(out),
            [](auto a, auto b) { return (a - b) * (a - b); });
        return std::accumulate(out.cbegin(), out.cend(), 0.) / (2*out.size());
    } };

    auto grad{ [&](const Eigen::MatrixXd& theta) {
        static int step = 0;
        ++step;
        step = step % reference.size();
        Lloyd2003Muscle::Parameters par = p;
       /* par.damping = theta(0);
        par.maxContractionVelocity = theta(1);
        par.maxIsometricForce = theta(2);
        par.optimalFiberLength = theta(3);
        par.pennationAngleAtOptimalFiberLength = theta(4);
        par.percentageChange = theta(5);
        par.strengthCoefficient = theta(6);
        par.tendonSlackLength = theta(7);
        */
        par.strengthCoefficient = theta(0);
       
        mtu.setParameters(par);

        Eigen::MatrixXd G(1, 1);
        G.setZero();

        mtu.setActivation(activation.at(step));
        mtu.setMusculotendonLength(mtuLength.at(step));
        if (step == 0) mtu.equilibrate();
        mtu.evaluate(0.005);
        mtu.calculateOutput();
        auto force = mtu.getOutput().force;
        auto Jac = mtu.calculateJacobian();
        G(0) = (force - reference.at(step)) * Jac(8);
 //       for (int k{0}; k < 8; ++k)
  //          G(k) = (force - reference.at(step)) * Jac(k + 2);
        return G;
    } };

    int epochs = 1000;
    int N = 1;

    Eigen::MatrixXd theta(1, 1);
  /*  theta(0) = p.damping;
    theta(1) = p.maxContractionVelocity;
    theta(2) = p.maxIsometricForce;
    theta(3) = p.optimalFiberLength;
    theta(4) = p.pennationAngleAtOptimalFiberLength;
    theta(5) = p.percentageChange;
    theta(6) = p.strengthCoefficient *0.5;
    theta(7) = p.tendonSlackLength;
    */
    theta(0) = p.strengthCoefficient * 0.5;
  
    Adam opt(N);
    opt.theta = theta;
    while (--epochs > 0) {
        auto g = grad(theta);
        theta = opt.step(g);
        if (epochs % 1 == 0) std::cout << J(theta) << std::endl;
    };
    std::cout << "Parameters: " << theta;


    return 0;
}

int main() {
    test3();

}

