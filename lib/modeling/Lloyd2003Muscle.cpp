#include "ceinms2/Lloyd2003Muscle.h"

namespace ceinms {

CubicSpline getDefaultTendonForceStrainCurve() {
    const std::vector<DoubleT> x{
        0, 0.01, 0.02, 0.06, 0.1, 0.14, 0.18, 0.24, 0.32, 0.38, 0.42, 0.44, 0.46, 0.48, 0.5
    };
    const std::vector<DoubleT> y{ 0,
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

    return CubicSpline(x, y);
}

CubicSpline getDefaultActiveForceLengthCurve() {
    const std::vector<DoubleT> x{ 0.4241,
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
    const std::vector<DoubleT> y{ 0,
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
    return CubicSpline(x, y);
}

CubicSpline getDefaultPassiveForceLengthCurve() {
    const std::vector<DoubleT> x{ 1, 1.116, 1.22, 1.324, 1.428, 1.61, 2 };
    const std::vector<DoubleT> y{ 0, 0.0208, 0.0604, 0.1536, 0.3117, 0.7448, 1.8571 };
    return CubicSpline(x, y);
}

CubicSpline getDefaultForceVelocityCurve() {
    const std::vector<DoubleT> x{
        -1, -0.624, -0.312, -0.104, 0, 0.104, 0.312, 0.624, 0.832, 0.988
    };
    const std::vector<DoubleT> y{
        0, 0.0939, 0.3174, 0.6162, 1, 1.2278, 1.2972, 1.3507, 1.3823, 1.4
    };
    return CubicSpline(x, y);
}



void Lloyd2003Muscle::setActivation(DoubleT activation) {
    i_.activation = activation;
}

void Lloyd2003Muscle::setMusculotendonLength(DoubleT musculotendonLength) {
    i_.musculotendonLength = musculotendonLength;
}

void Lloyd2003Muscle::setInput(Activation value) {
    setActivation(value.get());
}

void Lloyd2003Muscle::setInput(MusculotendonLength value) {
    setMusculotendonLength(value.get());
}

void Lloyd2003Muscle::setState(State state) {
    sNew_ = state;
}



DoubleT Lloyd2003Muscle::integrateFiberLength(DoubleT dt) {
    DoubleT minFiberLength = 0.2 * p_.optimalFiberLength;
    DoubleT maxFiberLength = 2 * p_.optimalFiberLength;

    auto f = [&](double fiberLength) -> double {
        const auto fiberVelocity = Lloyd2003Muscle::calculateFiberVelocityFromFiberLength(
            this->s_.fiberLength, fiberLength, dt);
        const auto tendonForce = Lloyd2003Muscle::calculateTendonForce(this->i_.musculotendonLength,
            fiberLength,
            this->p_.optimalFiberLength,
            this->p_.pennationAngleAtOptimalFiberLength,
            this->p_.tendonSlackLength,
            this->p_.strengthCoefficient,
            this->p_.maxIsometricForce,
            this->p_.tendonForceStrainCurve);
        const auto fiberForceProjectedOnTendon =
            Lloyd2003Muscle::calculateFiberForceProjectedOnTendon(this->i_.activation,
                fiberLength,
                fiberVelocity,
                this->p_.optimalFiberLength,
                this->p_.pennationAngleAtOptimalFiberLength,
                this->p_.maxIsometricForce,
                this->p_.strengthCoefficient,
                this->p_.damping,
                this->p_.percentageChange,
                this->p_.maxContractionVelocity,
                this->p_.activeForceLengthCurve,
                this->p_.passiveForceLengthCurve,
                this->p_.forceVelocityCurve);
        const auto v = tendonForce - fiberForceProjectedOnTendon;
        return v;
    };

    const auto currentFiberLength = wdbSolve(f, minFiberLength, maxFiberLength, 1e-9);
    return currentFiberLength;
}

void Lloyd2003Muscle::calculateOutput() {
    o_.force = calculateFiberForceProjectedOnTendon();
    if (properties_.doCalculateSecondaryOutput) {
        o_.normalizedFiberLengthAtT = calculateNormalizedFiberLengthAtT();
        o_.normalizedFiberLength = calculateNormalizedFiberLength();
        o_.normalizedFiberVelocity = calculateNormalizedFiberVelocity();
        o_.pennationAngle = calculatePennationAngle();
        o_.fa = p_.activeForceLengthCurve.get(o_.normalizedFiberLengthAtT);
        o_.fp = p_.passiveForceLengthCurve.get(o_.normalizedFiberLength);
        o_.fv = p_.forceVelocityCurve.get(o_.normalizedFiberVelocity);
        o_.activeForce =
            p_.maxIsometricForce * p_.strengthCoefficient * o_.fa * o_.fv * i_.activation;
        o_.passiveForce = p_.maxIsometricForce * p_.strengthCoefficient * o_.fp;
        o_.dampingForce =
            p_.maxIsometricForce * p_.strengthCoefficient * p_.damping * o_.normalizedFiberVelocity;
        o_.fiberForce = calculateFiberForce();
        o_.tendonStrain = calculateTendonStrain();
        o_.tendonLength = calculateTendonLength();
        o_.tendonForce = calculateTendonForce();
        o_.musculotendonStiffness = calculateMusculotendonStiffness();
        o_.fiberStiffness = calculateFiberStiffness();
        o_.tendonStiffness = calculateTendonStiffness();
    }
}

DoubleT Lloyd2003Muscle::calculateNormalizedFiberLengthAtT() const {
    return calculateNormalizedFiberLengthAtT(
        i_.activation, sNew_.fiberLength, p_.optimalFiberLength, p_.percentageChange);
}

DoubleT Lloyd2003Muscle::calculateNormalizedFiberLength() const {
    return calculateNormalizedFiberLength(sNew_.fiberLength, p_.optimalFiberLength);
}

DoubleT Lloyd2003Muscle::calculateNormalizedFiberVelocity() const {
    return calculateNormalizedFiberVelocity(
        sNew_.fiberVelocity, p_.optimalFiberLength, p_.maxContractionVelocity);
}

DoubleT Lloyd2003Muscle::calculatePennationAngle() const {
    return calculatePennationAngle(
        sNew_.fiberLength, p_.optimalFiberLength, p_.pennationAngleAtOptimalFiberLength);
}

DoubleT Lloyd2003Muscle::calculateTendonLength() const {
    return calculateTendonLength(i_.musculotendonLength,
        sNew_.fiberLength,
        p_.optimalFiberLength,
        p_.pennationAngleAtOptimalFiberLength);
}

DoubleT Lloyd2003Muscle::calculateTendonForce() const {
    return calculateTendonForce(i_.musculotendonLength,
        sNew_.fiberLength,
        p_.optimalFiberLength,
        p_.pennationAngleAtOptimalFiberLength,
        p_.tendonSlackLength,
        p_.strengthCoefficient,
        p_.maxIsometricForce,
        p_.tendonForceStrainCurve);
}

DoubleT Lloyd2003Muscle::calculateFiberForceProjectedOnTendon() const {
    return calculateFiberForceProjectedOnTendon(i_.activation,
        sNew_.fiberLength,
        sNew_.fiberVelocity,
        p_.optimalFiberLength,
        p_.pennationAngleAtOptimalFiberLength,
        p_.maxIsometricForce,
        p_.strengthCoefficient,
        p_.damping,
        p_.percentageChange,
        p_.maxContractionVelocity,
        p_.activeForceLengthCurve,
        p_.passiveForceLengthCurve,
        p_.forceVelocityCurve);
}

DoubleT Lloyd2003Muscle::calculateFiberStiffness() const {
    using boost::math::differentiation::make_fvar;
    const auto flen = make_fvar<DoubleT, 1>(sNew_.fiberLength);
    const auto force = calculateFiberForce(i_.activation,
        flen,
        sNew_.fiberVelocity,
        p_.maxIsometricForce,
        p_.strengthCoefficient,
        p_.damping,
        p_.optimalFiberLength,
        p_.percentageChange,
        p_.maxContractionVelocity,
        p_.activeForceLengthCurve,
        p_.passiveForceLengthCurve,
        p_.forceVelocityCurve);
    return force.derivative(1);
}

DoubleT Lloyd2003Muscle::calculateFiberStiffnessProjectedOnTendon() const {
    using boost::math::differentiation::make_fvar;
    const auto flen = make_fvar<DoubleT, 1>(sNew_.fiberLength);
    const auto force = calculateFiberForceProjectedOnTendon(i_.activation,
        flen,
        sNew_.fiberVelocity,
        p_.optimalFiberLength,
        p_.pennationAngleAtOptimalFiberLength,
        p_.maxIsometricForce,
        p_.strengthCoefficient,
        p_.damping,
        p_.percentageChange,
        p_.maxContractionVelocity,
        p_.activeForceLengthCurve,
        p_.passiveForceLengthCurve,
        p_.forceVelocityCurve);
    DoubleT pennationAngle = calculatePennationAngle(
        sNew_.fiberLength, p_.optimalFiberLength, p_.pennationAngleAtOptimalFiberLength);
    // To check
    return force.derivative(1) * cos(pennationAngle);
}

DoubleT Lloyd2003Muscle::calculateFiberForce() const {
    return calculateFiberForce(i_.activation,
        sNew_.fiberLength,
        sNew_.fiberVelocity,
        p_.maxIsometricForce,
        p_.strengthCoefficient,
        p_.damping,
        p_.optimalFiberLength,
        p_.percentageChange,
        p_.maxContractionVelocity,
        p_.activeForceLengthCurve,
        p_.passiveForceLengthCurve,
        p_.forceVelocityCurve);
}

DoubleT Lloyd2003Muscle::calculateTendonStrain() const {
    return calculateTendonStrain(i_.musculotendonLength,
        sNew_.fiberLength,
        p_.optimalFiberLength,
        p_.pennationAngleAtOptimalFiberLength,
        p_.tendonSlackLength);
}

DoubleT Lloyd2003Muscle::calculateMusculotendonStiffness() const {
    const DoubleT tendonStiffness = calculateTendonStiffness();
    const DoubleT fiberStiffnessProjectedOnTendon = calculateFiberStiffnessProjectedOnTendon();
    return 1. / (1. / tendonStiffness + 1. / fiberStiffnessProjectedOnTendon);
}

DoubleT Lloyd2003Muscle::calculateTendonStiffness() const {
    const auto tendonStrain{ calculateTendonStrain() };
    const auto tendonStiffness{ p_.strengthCoefficient * p_.maxIsometricForce / p_.tendonSlackLength
                                * p_.tendonForceStrainCurve.getFirstDerivative(tendonStrain) };
    return tendonStiffness;
}

auto Lloyd2003Muscle::calculateJacobian() const {
    constexpr int N = 1;
    using boost::math::differentiation::make_ftuple;
    const auto variables =
        make_ftuple<DoubleT, N, N, N, N, N, N, N, N, N, N>(i_.musculotendonLength,
            i_.activation,
            p_.damping,
            p_.maxContractionVelocity,
            p_.maxIsometricForce,
            p_.optimalFiberLength,
            p_.pennationAngleAtOptimalFiberLength,
            p_.percentageChange,
            p_.strengthCoefficient,
            p_.tendonSlackLength);

    const auto &fMusculotendonLength = std::get<0>(variables);
    const auto &fActivation = std::get<1>(variables);
    const auto &fDamping = std::get<2>(variables);
    const auto &fMaxContractionVelocity = std::get<3>(variables);
    const auto &fMaxIsometricForce = std::get<4>(variables);
    const auto &fOptimalFiberLength = std::get<5>(variables);
    const auto &fPennationAngleAtOptimalFiberLength = std::get<6>(variables);
    const auto &fPercentageChange = std::get<7>(variables);
    const auto &fStrengthCoefficient = std::get<8>(variables);
    const auto &fTendonSlackLength = std::get<9>(variables);

    auto force =
        [](const auto &aMusculotendonLength,
            const auto &aActivation,
            const auto &aFiberVelocity,
            const auto &aDamping,
            const auto &aMaxContractionVelocity,
            const auto &aMaxIsometricForce,
            const auto &aOptimalFiberLength,
            const auto &aPennationAngleAtOptimalFiberLength,
            const auto &aPercentageChange,
            const auto &aStrengthCoefficient,
            const auto &aTendonSlackLength,
            const CubicSpline &aActiveForceLengthCurve,
            const CubicSpline &aPassiveForceLengthCurve,
            const CubicSpline &aForceVelocityCurve) {
            const auto first = aOptimalFiberLength * sin(aPennationAngleAtOptimalFiberLength);
            const auto second = aMusculotendonLength - aTendonSlackLength;
            const auto aNewFiberLength = sqrt(first * first + second * second);

            // using fNewFiberLength we now have the calculation of the muscle force as function of
            // the musculotendon length, without having to do integration, which cannot work with
            // automatic differentiation
            const auto aForce = calculateFiberForceProjectedOnTendon(aActivation,
                aNewFiberLength,
                aFiberVelocity,
                aOptimalFiberLength,
                aPennationAngleAtOptimalFiberLength,
                aMaxIsometricForce,
                aStrengthCoefficient,
                aDamping,
                aPercentageChange,
                aMaxContractionVelocity,
                aActiveForceLengthCurve,
                aPassiveForceLengthCurve,
                aForceVelocityCurve);

            return aForce;
        }(fMusculotendonLength,
            fActivation,
            sNew_.fiberVelocity,
            fDamping,
            fMaxContractionVelocity,
            fMaxIsometricForce,
            fOptimalFiberLength,
            fPennationAngleAtOptimalFiberLength,
            fPercentageChange,
            fStrengthCoefficient,
            fTendonSlackLength,
            p_.activeForceLengthCurve,
            p_.passiveForceLengthCurve,
            p_.forceVelocityCurve);

    Eigen::Matrix<DoubleT, 1, 10> J;
    J(0) = force.at(1, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    J(1) = force.at(0, 1, 0, 0, 0, 0, 0, 0, 0, 0);
    J(2) = force.at(0, 0, 1, 0, 0, 0, 0, 0, 0, 0);
    J(3) = force.at(0, 0, 0, 1, 0, 0, 0, 0, 0, 0);
    J(4) = force.at(0, 0, 0, 0, 1, 0, 0, 0, 0, 0);
    J(5) = force.at(0, 0, 0, 0, 0, 1, 0, 0, 0, 0);
    J(6) = force.at(0, 0, 0, 0, 0, 0, 1, 0, 0, 0);
    J(7) = force.at(0, 0, 0, 0, 0, 0, 0, 1, 0, 0);
    J(8) = force.at(0, 0, 0, 0, 0, 0, 0, 0, 1, 0);
    J(9) = force.at(0, 0, 0, 0, 0, 0, 0, 0, 0, 1);
    return J;
}

void Lloyd2003Muscle::equilibrate() {
    const double tol = 1e-9;
    do {
        integrate(0.001);
        validateState();
    } while (sNew_.fiberVelocity > tol);
}

void Lloyd2003Muscle::integrate(DoubleT dt) {
    sNew_.fiberLength = integrateFiberLength(dt);
    sNew_.fiberVelocity =
        calculateFiberVelocityFromFiberLength(s_.fiberLength, sNew_.fiberLength, dt);
}

void Lloyd2003Muscle::validateState() {
    s_ = sNew_;
}

void Lloyd2003Muscle::evaluate(DoubleT dt) {
    integrate(dt);
    calculateOutput();
    validateState();
}

void connectSocket(const Activation &parent, Lloyd2003Muscle &child) {
    child.setInput(parent);
}

void connectSocket(const MusculotendonLength &parent, Lloyd2003Muscle &child) {
    child.setInput(parent);
}

}// namespace ceinms