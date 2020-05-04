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
using std::cout;
using std::endl;
using std::vector;
using std::string;

constexpr double pi = 3.14159265358979323846;
using DoubleT = double;
vector<string> tokenize(const string &line, char sep = ',') {
    vector<string> result;
    size_t first(0), last(0);
    if (line.find(sep) == string::npos)
        result.emplace_back(line);
    else {
        while (last < line.size()) {
            if (line[last] == '\"') {
                last++;
                last =
                    line.substr(last, line.size() - last).find('\"') + last + 1;
            } else if (line[last] == sep) {
                result.emplace_back(line.substr(first, last - first));
                first = ++last;
            } else
                ++last;
        }
        if (last != first)
            result.emplace_back(line.substr(first, last - first));
    }
    return result;
}

template<typename T>
static ceinms::DataTable<T> fillFromCSV(const string filename,
    const char sep = ',',
    unsigned skiprows = 0) {
    ceinms::DataTable<T> data;
    std::ifstream inF(filename);
    if (!inF.is_open()) {
        std::cout << "Cannot open " << filename << std::endl;
        exit(EXIT_FAILURE);
    }

    string line;
    for (unsigned i = 0; i < skiprows; ++i) getline(inF, line, '\n');
    getline(inF, line, '\n');
    auto labels(tokenize(line, sep));
    data.setLabels(vector<string>(labels.begin() + 1, labels.end()));

    while (!inF.eof()) {
        getline(inF, line);
        if (!line.empty()) {
            vector<string> result(tokenize(line, sep));
            vector<T> values;
            std::transform(std::begin(result),
                std::end(result),
                std::back_inserter(values),
                [](std::string v) { return std::stod(v); });
            data.pushRow(
                values.at(0), std::vector<T>(values.begin() + 1, values.end()));
        }
    }
    return data;
}

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
void runTrial(ceinms::DataTable<double> &trial) {
    const vector<DoubleT> act_x{ 0.4241,
        0.4641,
        0.5441,
        0.6241,
        0.7041,
        0.7441,
        0.8641,
        0.9441,
        0.9841,
        1.0641,
        1.1041,
        1.2241,
        1.7841,
        1.8241 };
    const vector<DoubleT> act_y{ 0,
        0.0036,
        0.2531,
        0.5362,
        0.7739,
        0.8204,
        0.926,
        0.9919,
        1,
        0.9907,
        0.9524,
        0.7902,
        0.0117,
        0 };

    vector<DoubleT> pas_x{ 1, 1.116, 1.22, 1.324, 1.428, 1.61, 2 };
    vector<DoubleT> pas_y{ 0, 0.0208, 0.0604, 0.1536, 0.3117, 0.7448, 1.8571 };

    vector<DoubleT> vel_x{
        -1, -0.624, -0.312, -0.104, 0, 0.104, 0.312, 0.624, 0.832, 0.988
    };
    vector<DoubleT> vel_y{
        0, 0.0939, 0.3174, 0.6162, 1, 1.2278, 1.2972, 1.3507, 1.3823, 1.4
    };

    vector<DoubleT> ten_x{ 0,
        0.01,
        0.02,
        0.06,
        0.1,
        0.14,
        0.18,
        0.24,
        0.32,
        0.38,
        0.42,
        0.44,
        0.46,
        0.48,
        0.5 };
    vector<DoubleT> ten_y{ 0,
        0.094,
        0.2609,
        1.3087,
        2.4311,
        3.5536,
        4.676,
        6.3597,
        8.6046,
        10.2883,
        11.4107,
        11.9719,
        12.5332,
        13.0944,
        13.6556 };

    ceinms::CubicSpline act( act_x, act_y );
    ceinms::CubicSpline pas( pas_x, pas_y );
    ceinms::CubicSpline vel( vel_x, vel_y );
    ceinms::CubicSpline ten( ten_x, ten_y );

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
    vector<DoubleT> outForce;
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
    }
    trial.pushColumn("Predicted_Force_N", outForce);
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
        runTrial<MuscleT>(currentTrial);
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
        cout << "AbsoluteRelativeError @ " << disp << " mm: " <<  absoluteRelativeError << "%" << endl;
        currentTrial.print("out" + trialName);
    }
    return failed;
}

int main() {
    bool failed = runMillardBenchmark<ceinms::Lloyd2003Muscle>();
    return failed;
}
