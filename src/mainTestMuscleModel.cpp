#include "ceinms2/Lloyd2019Activation.h"
#include "ceinms2/Lloyd2019Muscle.h"
#include "ceinms2/Curve.h"
#include "ceinms2/DataTable.h"
#include <cmath>
#include <vector>
#include <fstream>
#include <filesystem>
#include <map>
#include <algorithm>
using std::cout;
using std::endl;
using std::vector;
using std::string;

constexpr double pi = 3.14159265358979323846;

vector<string> tokenize(const string &line, char sep = ',') {
    vector<string> result;
    size_t first(0), last(0);
    if (line.find(sep) == string::npos)
        result.emplace_back(line);
    else {
        while (last < line.size()) {
            if (line[last] == '\"') {
                last++;
                last = line.substr(last, line.size() - last).find('\"') + last + 1;
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
static ceinms::DataTable<T> fillFromCSV(const string filename, const char sep = ',', unsigned skiprows = 0) {
    ceinms::DataTable<T> data;
    std::ifstream inF(filename);
    if (!inF.is_open()) {
        std::cout << "Cannot open " << filename << std::endl;
        exit(EXIT_FAILURE);
    }

    string line;
    for (unsigned i = 0; i < skiprows; ++i)
        getline(inF, line, '\n');
    getline(inF, line, '\n');
    auto labels(tokenize(line, sep));
    data.setLabels(vector<string>(labels.begin() + 1, labels.end()));

    while (!inF.eof()) {
        getline(inF, line);
        if (!line.empty()) {
            vector<string> result(tokenize(line, sep));
            vector<T> values;
            std::transform(std::begin(result), std::end(result), std::back_inserter(values), [](std::string v) {
                return std::stod(v);
            });
            data.pushRow(values.at(0), std::vector<T>(values.begin() + 1, values.end()));
        }
    }
    return data;
}

/*
void printStates(ceinms::Lloyd2019Muscle& muscle) {

	ofstream outF("out.csv");
	outF << "calcFibreForce, activation, activeForce, dampingForce, fibreForce,fibreLength, musculotendonLength,"
		"normalisedFibreLength, normalisedFibreLengthAtT, normalisedFibreVelocity, optimalFibreLengthAtT,"
		"passiveForce, pennationAngle, tendonLength, tendonStrain, fa, fv, fp" << endl;
	DoubleT flStart = 0.1, flStop = 1.9;
	const unsigned N = 100;
	DoubleT dx = (flStop - flStart) / N;
	for (int i{ 0 }; i < N; ++i) {
		outF << muscle.calculateFibreForce(1., (i+1)*dx, 0.) << ", ";
		outF << muscle.updStates().activation << ", ";
		outF << muscle.updStates().activeForce << ", ";
		outF << muscle.updStates().dampingForce << ", ";
		outF << muscle.updStates().fibreForce << ", ";
		outF << muscle.updStates().fibreLength << ", ";
		outF << muscle.updStates().musculotendonLength << ", ";
		outF << muscle.updStates().normalisedFibreLength << ", ";
		outF << muscle.updStates().normalisedFibreLengthAtT << ", ";
		outF << muscle.updStates().normalisedFibreVelocity << ", ";
		outF << muscle.updStates().optimalFibreLengthAtT << ", ";
		outF << muscle.updStates().passiveForce << ", ";
		outF << muscle.updStates().pennationAngle << ", ";
		outF << muscle.updStates().tendonLength << ", ";
		outF << muscle.updStates().tendonStrain << ", ";
		outF << muscle.updStates().fa << ", ";
		outF << muscle.updStates().fv << ", ";
		outF << muscle.updStates().fp << endl;

	}
	outF.close();
}
*/


void printSolutionSpace(ceinms::Lloyd2019Muscle &muscle) {
    std::ofstream outF("outSolutionSpace.csv");
    DoubleT flStart = 0.5, flStop = 1.5, fvStart = -1, fvStop = 1;
    const unsigned N = 100;
    DoubleT dFl = (flStop - flStart) / N;
    DoubleT dFv = (fvStop - fvStart) / N;
    for (int iFl(0); iFl < N; ++iFl)
        for (int iFv(0); iFv < N; ++iFv) {
            DoubleT fl = flStart + dFl * iFl;
            DoubleT fv = fvStart + dFv * iFv;
            DoubleT u = muscle.calculateFibreForce(1., fl, fv, muscle.updParameters());
            DoubleT t = muscle.calculateTendonForce(2, fl, muscle.updParameters());
            outF << fl << ", ";
            outF << fv << ", ";
            outF << u << ", ";
            outF << t << ", ";
            outF << DoubleT(u - t) << endl;
        }
}

void printStep(ceinms::Lloyd2019Muscle &muscle) {
    std::ofstream outF("outStepResponse.csv");
    DoubleT mtuLStart = 2.1, mtuLStop = 2.3, activation = 0.5;
    muscle.setInput(activation, mtuLStart);
    muscle.equilibrate();
    for (int i{ 0 }; i < 100; ++i) {
        muscle.setInput(activation, mtuLStart);
        muscle.integrate(0.001);
        muscle.validateState();
        muscle.calculateOutput();
        outF << muscle.updOutput().fibreForce << ", " << muscle.updOutput().normalisedFibreLength << endl;
    }

    for (int i{ 0 }; i < 100; ++i) {
        muscle.setInput(activation, mtuLStop);
        muscle.integrate(0.001);
        muscle.validateState();
        muscle.calculateOutput();
        outF << muscle.updOutput().fibreForce << ", " << muscle.updOutput().normalisedFibreLength << endl;
    }
}

/*
void testAutoDiff(ceinms::Lloyd2019Muscle& muscle) {
	dual a = 1.;
	dual fl = 1.;
	dual fv = 0;
	dual u = muscle.calculateFibreForce(1, 1, 0);
	cout << u << endl;
	dual dudfl = derivative([&](auto x1, auto x2, auto x3) {return muscle.calculateFibreForce(x1, x2, x3); }, wrt(fl), a, fl, fv);
	cout << "automatic derivative: " << dudfl << endl;

	dual dx = 0.0001;
	dual d1 = muscle.calculateFibreForce(1, 1 + dx, 0);
	dual d2 = muscle.calculateFibreForce(1, 1 - dx, 0);
	dual diff = (d1 - d2) / (2 * dx);
	cout << "Numeric differentiation: " << diff << endl;
}
*/

int main_() {
    vector<DoubleT> act_x{ 0.4241, 0.4641, 0.5441, 0.6241, 0.7041, 0.7441, 0.8641, 0.9441, 0.9841, 1.0641, 1.1041, 1.2241, 1.7841, 1.8241 };
    vector<DoubleT> act_y{ 0, 0.0036, 0.2531, 0.5362, 0.7739, 0.8204, 0.926, 0.9919, 1, 0.9907, 0.9524, 0.7902, 0.0117, 0 };

    vector<DoubleT> pas_x{ 1, 1.116, 1.22, 1.324, 1.428, 1.61, 2 };
    vector<DoubleT> pas_y{ 0, 0.0208, 0.0604, 0.1536, 0.3117, 0.7448, 1.8571 };

    vector<DoubleT> vel_x{ -1, -0.624, -0.312, -0.104, 0, 0.104, 0.312, 0.624, 0.832, 0.988 };
    vector<DoubleT> vel_y{ 0, 0.0939, 0.3174, 0.6162, 1, 1.2278, 1.2972, 1.3507, 1.3823, 1.4 };

    vector<DoubleT> ten_x{ 0, 0.01, 0.02, 0.06, 0.1, 0.14, 0.18, 0.24, 0.32, 0.38, 0.42, 0.44, 0.46, 0.48, 0.5 };
    vector<DoubleT> ten_y{ 0, 0.094, 0.2609, 1.3087, 2.4311, 3.5536, 4.676, 6.3597, 8.6046, 10.2883, 11.4107, 11.9719, 12.5332, 13.0944, 13.6556 };

    CurveOffline act(act_x, act_y);
    CurveOffline pas(pas_x, pas_y);
    CurveOffline vel(vel_x, vel_y);
    CurveOffline ten(ten_x, ten_y);

    ceinms::Lloyd2019Muscle::Parameters p;
    p.damping = 0.1;
    p.maxContractionVelocity = 1;
    p.maxIsometricForce = 1;
    p.optimalFibreLength = 1;
    p.pennationAngleAtOptimalFibreLength = 0;
    p.percentageChange = 0.15;
    p.strengthCoefficient = 1;
    p.tendonSlackLength = 1;
    p.activeForceLengthCurve = act;
    p.passiveForceLengthCurve = pas;
    p.forceVelocityCurve = vel;
    p.tendonForceStrainCurve = ten;

    ceinms::Lloyd2019Muscle muscle(p);
    printSolutionSpace(muscle);
    printStep(muscle);

    return 0;
}


int main_1() {
    /*
	using DoubleT = dual;
	DoubleT a = 1;
	a = a - 2*a;
	cout << "it's " << a << " should be -1\n";

	a = 1;
	a = a + 2 * a;
	cout << "it's " << a << " should be 3\n";

	a = 1;
	a = a - a;
	cout << "it's " << a << " should be 0\n";

	a = 1;
	a = a + a;
	cout << "it's " << a << " should be 2\n";

	a = 1;
	a = 2*a - a;
	cout << "it's " << a << " should be 1\n";

	a = 1;
	a = 2*a + a;
	cout << "it's " << a << " should be 3\n";

	a = 1;
	a = DoubleT{ 2 * a - a };
	cout << "it's " << a << " should be 1\n";
	 */
    return 0;
}


void runTrial(ceinms::DataTable<double> &trial) {
    vector<DoubleT> act_x{ 0.4241, 0.4641, 0.5441, 0.6241, 0.7041, 0.7441, 0.8641, 0.9441, 0.9841, 1.0641, 1.1041, 1.2241, 1.7841, 1.8241 };
    vector<DoubleT> act_y{ 0, 0.0036, 0.2531, 0.5362, 0.7739, 0.8204, 0.926, 0.9919, 1, 0.9907, 0.9524, 0.7902, 0.0117, 0 };

    vector<DoubleT> pas_x{ 1, 1.116, 1.22, 1.324, 1.428, 1.61, 2 };
    vector<DoubleT> pas_y{ 0, 0.0208, 0.0604, 0.1536, 0.3117, 0.7448, 1.8571 };

    vector<DoubleT> vel_x{ -1, -0.624, -0.312, -0.104, 0, 0.104, 0.312, 0.624, 0.832, 0.988 };
    vector<DoubleT> vel_y{ 0, 0.0939, 0.3174, 0.6162, 1, 1.2278, 1.2972, 1.3507, 1.3823, 1.4 };

    vector<DoubleT> ten_x{ 0, 0.01, 0.02, 0.06, 0.1, 0.14, 0.18, 0.24, 0.32, 0.38, 0.42, 0.44, 0.46, 0.48, 0.5 };
    vector<DoubleT> ten_y{ 0, 0.094, 0.2609, 1.3087, 2.4311, 3.5536, 4.676, 6.3597, 8.6046, 10.2883, 11.4107, 11.9719, 12.5332, 13.0944, 13.6556 };

    CurveOffline act(act_x, act_y);
    CurveOffline pas(pas_x, pas_y);
    CurveOffline vel(vel_x, vel_y);
    CurveOffline ten(ten_x, ten_y);

    ceinms::Lloyd2019Muscle::Parameters p;
    p.damping = 0.1;
    p.maxContractionVelocity = 5.3156293513;
    p.maxIsometricForce = 1.2758948506;
    p.optimalFibreLength = 17.1e-3;//m;
    p.pennationAngleAtOptimalFibreLength = 6.0 * (pi / 180.0);//rad
    p.percentageChange = 0.15;
    p.strengthCoefficient = 1;
    p.tendonSlackLength = p.optimalFibreLength;
    p.activeForceLengthCurve = act;
    p.passiveForceLengthCurve = pas;
    p.forceVelocityCurve = vel;
    p.tendonForceStrainCurve = ten;

    DoubleT muscleRestLength = p.optimalFibreLength * std::cos(p.pennationAngleAtOptimalFibreLength) + p.tendonSlackLength;
    DoubleT posOffset = -2.0e-3;//m
    ceinms::Lloyd2019Muscle muscle(p);

    auto times = trial.getTimeColumn();
    auto displacements = trial.getColumn("Experimental_displacement_mm");
    vector<DoubleT> outForce;
    for (size_t i{ 0 }; i < times.size(); ++i) {
        DoubleT disp = displacements.at(i);
        DoubleT dt = 0.001;
        if (i > 0)
            dt = times.at(i) - times.at(i - 1);
        DoubleT mtLength = muscleRestLength + posOffset + disp * 1.0e-3;
        muscle.setInput(1.0, mtLength);
        if (i == 0)
            muscle.equilibrate();
        muscle.integrate(dt);
        muscle.validateState();
        muscle.calculateOutput();
        outForce.push_back(muscle.updOutput().tendonForce);
    }
    trial.pushColumn("Predicted_Force_N", outForce);
}

void runMillardBenchmark() {
    std::filesystem::path root(ROOT_DIR);
    std::filesystem::path dataPath{ root / std::filesystem::path("data/millard2013/maximalActivation") };
    const std::map<std::string, double> trialToDisplacement{
        { "force_trial1.dat", 0.05 },
        { "force_trial2.dat", 0.1 },
        { "force_trial3.dat", 0.25 },
        { "force_trial4.dat", 0.5 },
        { "force_trial5.dat", 1.0 },
        { "force_trial6.dat", 2.0 },
    };
    auto normalisedDisplacementDataTable = fillFromCSV<double>((dataPath / "displacement.dat").string(), '\t', 15);
    auto normalisedDisplacement = normalisedDisplacementDataTable.getColumn("displacement_mm");
    auto times = normalisedDisplacementDataTable.getTimeColumn();
    CurveOffline displacementCurve{ times, normalisedDisplacement };

    for (auto &[trialName, disp] : trialToDisplacement) {
        auto currentTrial = fillFromCSV<double>((dataPath / trialName).string(), '\t', 15);
        currentTrial.setLabels({ "Experimental_Force_N" });
        std::vector<double> currentDisplacementData;
        for (auto &t : currentTrial.getTimeColumn())
            currentDisplacementData.emplace_back(displacementCurve.getValue(t) * disp);
        currentTrial.pushColumn("Experimental_displacement_mm", currentDisplacementData);
        runTrial(currentTrial);
        currentTrial.print("out" + trialName);
    }
}

int main() {
    runMillardBenchmark();
    return 0;
}