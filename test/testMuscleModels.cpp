#include "ceinms2/Lloyd2003Muscle.h"
#include "ceinms2/Curve.h"
#include "ceinms2/DataTable.h"
#include <cmath>
#include <vector>
#include <fstream>
#include <filesystem>
#include <map>
#include <algorithm>
#include <numeric>
#include <chrono>
using std::cout;
using std::endl;
using std::vector;
using std::string;

constexpr double pi = 3.14159265358979323846;

template<typename MuscleT>
void printStep(MuscleT &muscle) {
    std::ofstream outF("outStepResponse.csv");
    DoubleT mtuLStart = 2.1, mtuLStop = 2.3, activation = 0.5;
    muscle.setActivation(activation);
    muscle.setMusculotendonLength(mtuLStart);
    muscle.equilibrate();
    for (int i{ 0 }; i < 100; ++i) {
        muscle.setActivation(activation);
        muscle.setMusculotendonLength(mtuLStart);
        muscle.evaluate(0.001);
        outF << muscle.getOutput().fiberForce << ", "
             << muscle.getOutput().normalisedFiberLength << endl;
    }

    for (int i{ 0 }; i < 100; ++i) {
        muscle.setActivation(activation);
        muscle.setMusculotendonLength(mtuLStop);
        muscle.evaluate(0.001);
        outF << muscle.getOutput().fiberForce << ", "
             << muscle.getOutput().normalisedFiberLength << endl;
    }
}

template<typename MuscleT>
double runTrial(ceinms::DataTable<double> &trial) {
   
    const ceinms::CubicSpline act{ ceinms::getDefaultActiveForceLengthCurve() };
    const ceinms::CubicSpline pas( ceinms::getDefaultPassiveForceLengthCurve());
    const ceinms::CubicSpline vel( ceinms::getDefaultForceVelocityCurve() );
    const ceinms::CubicSpline ten( ceinms::getDefaultTendonForceStrainCurve() );

    typename MuscleT::Parameters p;
    p.damping = 0.1;
    p.maxContractionVelocity = 5.3156293513;
    p.maxIsometricForce = 1.2758948506;
    p.optimalFiberLength = 17.1e-3;// m;
    p.pennationAngleAtOptimalFiberLength = 6.0 * (pi / 180.0);// rad
    p.percentageChange = 0.15;
    p.strengthCoefficient = 1;
    p.tendonSlackLength = p.optimalFiberLength;
    p.activeForceLengthCurve = act;
    p.passiveForceLengthCurve = pas;
    p.forceVelocityCurve = vel;
    p.tendonForceStrainCurve = ten;

    DoubleT muscleRestLength =
        p.optimalFiberLength * std::cos(p.pennationAngleAtOptimalFiberLength)
        + p.tendonSlackLength;
    DoubleT posOffset = -2.0e-3;// m
    MuscleT muscle(p);

    auto times = trial.getTimeColumn();
    auto displacements = trial.getColumn("Experimental_displacement_mm");
    vector<DoubleT> outForce, outStiffness;
    auto startTime = std::chrono::system_clock::now();
    for (size_t i{ 0 }; i < times.size(); ++i) {
        DoubleT disp = displacements.at(i);
        DoubleT dt = 0.001;
        if (i > 0) dt = times.at(i) - times.at(i - 1);
        DoubleT mtLength = muscleRestLength + posOffset + disp * 1.0e-3;
        muscle.setActivation(1.0);
        muscle.setMusculotendonLength(mtLength);
        if (i == 0) muscle.equilibrate();
        muscle.integrate(dt);
        muscle.validateState();
        muscle.calculateOutput();
        outForce.push_back(muscle.getOutput().tendonForce);
        outStiffness.push_back(muscle.getOutput().musculotendonStiffness);
    }
    auto stopTime = std::chrono::system_clock::now();
    trial.pushColumn("Predicted_Force_N", outForce);
    trial.pushColumn("Predicted_MTUstiffness_N", outStiffness);
    return std::chrono::duration_cast<std::chrono::microseconds>(stopTime - startTime).count()/1000.;
}

/*This test runs the Maximal-Activation Biological Benchmark test described in
Millard M, Uchida T, Seth A, Delp SL. Flexing computational muscle: modeling and simulation 
of musculotendon dynamics. J Biomech Eng 2013; 135: 021005.

In their paper, Millard et al. report mean absolute errors relative to experimentally measured 
maximum isometric force between 3.3% and 8.9%. This test will fail if the calculated error for this same test
is greater than 10%.
*/
template<typename MuscleT>
bool runMillardBenchmark() {
    // 1.17 is the experimentally measured maximum isometric force. See Millard et al. 2013
    double experimentallyMeasuredFiso = 1.17;
    std::filesystem::path root(ROOT_DIR);
    std::filesystem::path dataPath{
        root / std::filesystem::path("data/millard2013/maximalActivation")
    };

    bool failed(false);
    const std::map<std::string, double> trialToDisplacement{
        { "force_trial1.dat", 0.05 },
        { "force_trial2.dat", 0.1 },
        { "force_trial3.dat", 0.25 },
        { "force_trial4.dat", 0.5 },
        { "force_trial5.dat", 1.0 },
        { "force_trial6.dat", 2.0 },
    };
    auto normalisedDisplacementDataTable =
        fillFromCSV<double>((dataPath / "displacement.dat").string(), '\t', 15);
    auto normalisedDisplacement =
        normalisedDisplacementDataTable.getColumn("displacement_mm");
    auto times = normalisedDisplacementDataTable.getTimeColumn();
    ceinms::Curve<ceinms::CurveMode::Offline> displacementCurve{ times, normalisedDisplacement };

    auto squareError = [](double a, double b) {
        double e = a - b;
        return e * e;
    };

       auto absoluteError = [](double a, double b) {
        return std::abs(a - b);
    };
    std::vector<double> runTimes;
    for (auto &[trialName, disp] : trialToDisplacement) {
        auto currentTrial =
            fillFromCSV<double>((dataPath / trialName).string(), '\t', 15);
        currentTrial.setLabels({ "Experimental_Force_N" });
        std::vector<double> currentDisplacementData;
        for (auto &t : currentTrial.getTimeColumn())
            currentDisplacementData.emplace_back(
                displacementCurve.getValue(t) * disp);
        currentTrial.pushColumn(
            "Experimental_displacement_mm", currentDisplacementData);
        double elapsedTime = runTrial<MuscleT>(currentTrial);
        runTimes.push_back(elapsedTime);
        auto experimentalForce = currentTrial.getColumn("Experimental_Force_N");
        auto predictedForce = currentTrial.getColumn("Predicted_Force_N");
        
        double sum = std::transform_reduce(experimentalForce.begin(),
            experimentalForce.end(),
            predictedForce.begin(),
            0.,
            std::plus<>(),
            absoluteError);
 
        auto absoluteRelativeError =
            sum / experimentalForce.size() * 100. / experimentallyMeasuredFiso;
        failed |= (absoluteRelativeError > 10.);
        cout << "AbsoluteRelativeError @ " << disp << " mm: " <<  absoluteRelativeError << "% calculated in " << elapsedTime << " ms" << endl;
        currentTrial.print("out" + trialName);
    }
    return failed;
}

int main() {
    bool failed = runMillardBenchmark<ceinms::Lloyd2003Muscle>();
    return failed;
}
