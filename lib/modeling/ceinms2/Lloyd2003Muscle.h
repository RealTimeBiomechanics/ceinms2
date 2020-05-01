#ifndef ceinms2_Lloyd2003Muscle_h
#define ceinms2_Lloyd2003Muscle_h
#include <cmath>
#include <string_view>
#include "ceinms2/Types.h"
#include "ceinms2/Curve.h"
#include "ceinms2/WDBsolver.h"
#include "ceinms2/ExponentialActivation.h"

namespace ceinms {

using CurveOffline = ceinms::Curve<ceinms::CurveMode::Offline,
    ceinms::CurveMode::Interpolation::Cubic>;
class Lloyd2003Muscle {
  public:
    using type = nullptr_t;
    using concept_t = component_t;
    static constexpr std::string_view class_name = "Lloyd2003Muscle";
    struct Parameters {
        Parameters() = default;
        DoubleT optimalFiberLength{ 1. };
        DoubleT pennationAngleAtOptimalFiberLength{ 0. };
        DoubleT tendonSlackLength{ 0.8 };
        DoubleT maxContractionVelocity{ 5. };
        DoubleT damping{ 0.1 };
        DoubleT maxIsometricForce{ 100. };
        DoubleT strengthCoefficient{ 1. };
        DoubleT percentageChange{ 0.15 };
        CurveOffline forceVelocityCurve{};
        CurveOffline activeForceLengthCurve{};
        CurveOffline passiveForceLengthCurve{};
        CurveOffline tendonForceStrainCurve{};
    };

    struct Input {
        Input() = default;
        DoubleT activation{ 0. };
        DoubleT musculotendonLength{ 0. };
    };

    struct State {
        State() = default;
        DoubleT fiberLength{ 0. };
        DoubleT fiberVelocity{ 0. };
    };

    struct Output {
        Output() = default;
        DoubleT optimalFiberLengthAtT{ 0. };
        DoubleT pennationAngle{ 0. };
        DoubleT fiberForce{ 0. };
        DoubleT activeForce{ 0. };
        DoubleT passiveForce{ 0. };
        DoubleT dampingForce{ 0. };
        DoubleT normalisedFiberLength{ 0. };
        DoubleT normalisedFiberLengthAtT{ 0. };
        DoubleT normalisedFiberVelocity{ 0. };
        DoubleT tendonLength{ 0. };
        DoubleT tendonStrain{ 0. };
        DoubleT tendonForce{ 0. };
        DoubleT fa{ 0. };
        DoubleT fv{ 0. };
        DoubleT fp{ 0. };
    };

    Lloyd2003Muscle(Parameters parameters)
        : p_(parameters) {
        s_.fiberLength = p_.optimalFiberLength;
        s_.fiberVelocity = 0.;
    }

    void setActivation(DoubleT activation);
    void setMusculotendonLength(DoubleT musculotendonLength);
    void setInput(Activation value);
    void setInput(MusculotendonLength value);

    void equilibrate();
    void setState(State state);
    // From input and current state calculate the new state of the system
    void integrate(DoubleT dt);
    // The temporary state calculated via `integrate` becomes the new state
    void validateState();
    // from the internal state of the system and the input, calculate all the
    // output;
    void calculateOutput();
    [[nodiscard]] std::string getName() const { return name_; }
    void setName(std::string name) { name_ = name; }
    // Convenience function that, from input and current state, calculate the
    // new state and all the output
    void evaluate(DoubleT dt);

    [[nodiscard]] Parameters &getParameters() { return p_; }
    [[nodiscard]] const Parameters& getParameters() const { return p_; }
    [[nodiscard]] const State &getState() const { return s_; }
    [[nodiscard]] State &getState() { return s_; }
    [[nodiscard]] const Output &getOutput() const { return o_; }
    template<typename T, std::enable_if_t<std::is_same<T, Force>::value, int> = 0>
    [[nodiscard]] Force getOutput() const {
        return Force{ o_.fiberForce };
    }
    /* calculateXYZ functions perform calculations and return the calculated
     * value*/
    static DoubleT calculateFiberVelocityFromFiberLength(DoubleT previousFiberLength,
        DoubleT currentFiberLength,
        DoubleT dt);
    static DoubleT calculateFiberForce(DoubleT activation,
        DoubleT fiberLength,
        DoubleT fiberVelocity,
        const Parameters &p);
    static DoubleT
        calculateTendonForce(DoubleT musculotendonLength, DoubleT fiberLength, const Parameters &p);
    static DoubleT calculatePennationAngle(DoubleT fiberLength, const Parameters &p);
    static DoubleT calculateOptimalFiberLengthAtT(DoubleT activation, const Parameters &p);
    static DoubleT calculateNormalisedFiberLengthAtT(DoubleT activation,
        DoubleT fiberLength,
        const Parameters &p);
    static DoubleT calculateNormalisedFiberLength(DoubleT fiberLength, const Parameters &p);
    static DoubleT calculateNormalisedFiberVelocity(DoubleT fiberVelocity, const Parameters &p);
    static DoubleT calculateTendonLength(DoubleT musculotendonLength,
        DoubleT fiberLength,
        const Parameters &p);
    static DoubleT calculateTendonStrain(DoubleT musculotendonLength,
        DoubleT fiberLength,
        const Parameters &p);

  private:
    DoubleT integrateFiberLength(DoubleT dt);
    std::string name_;
    Parameters p_;
    State s_, sNew_;
    Input i_;
    Output o_;
};

void Lloyd2003Muscle::setActivation(DoubleT activation) { i_.activation = activation; }

void Lloyd2003Muscle::setMusculotendonLength(DoubleT musculotendonLength) {
    i_.musculotendonLength = musculotendonLength;
}

void Lloyd2003Muscle::setInput(Activation value) { setActivation(value.get()); }

void Lloyd2003Muscle::setInput(MusculotendonLength value) { setMusculotendonLength(value.get()); }

void Lloyd2003Muscle::setState(State state) { s_ = state; }

void Lloyd2003Muscle::equilibrate() {
    double diff = 1;
    double fLength = s_.fiberLength;
    double tol = 1e-9;
    while (diff > tol) {
        integrate(0.001);
        diff = std::abs(sNew_.fiberLength - fLength);
        fLength = sNew_.fiberLength;
    }
    validateState();
}

void Lloyd2003Muscle::integrate(DoubleT dt) {
    sNew_.fiberLength = integrateFiberLength(dt);
    sNew_.fiberVelocity =
        calculateFiberVelocityFromFiberLength(s_.fiberLength, sNew_.fiberLength, dt);
}

void Lloyd2003Muscle::validateState() { s_ = sNew_; }

void Lloyd2003Muscle::calculateOutput() {
    o_.normalisedFiberLengthAtT =
        calculateNormalisedFiberLengthAtT(i_.activation, s_.fiberLength, p_);
    o_.normalisedFiberLength = calculateNormalisedFiberLength(s_.fiberLength, p_);
    o_.normalisedFiberVelocity = calculateNormalisedFiberVelocity(s_.fiberVelocity, p_);
    o_.pennationAngle = calculatePennationAngle(s_.fiberLength, p_);
    o_.fa = p_.activeForceLengthCurve.getValue(o_.normalisedFiberLengthAtT);
    o_.fp = p_.passiveForceLengthCurve.getValue(o_.normalisedFiberLength);
    o_.fv = p_.forceVelocityCurve.getValue(o_.normalisedFiberVelocity);
    o_.activeForce = p_.maxIsometricForce * p_.strengthCoefficient * o_.fa * o_.fv;
    o_.passiveForce = p_.maxIsometricForce * p_.strengthCoefficient * o_.fp;
    o_.dampingForce =
        p_.maxIsometricForce * p_.strengthCoefficient * p_.damping * o_.normalisedFiberVelocity;
    o_.fiberForce = calculateFiberForce(i_.activation, s_.fiberLength, s_.fiberVelocity, p_);
    o_.tendonStrain = calculateTendonStrain(i_.musculotendonLength, s_.fiberLength, p_);
    o_.tendonLength = calculateTendonLength(i_.musculotendonLength, s_.fiberLength, p_);
    o_.tendonForce = calculateTendonForce(i_.musculotendonLength, s_.fiberLength, p_);
}

void Lloyd2003Muscle::evaluate(DoubleT dt) {
    integrate(dt);
    validateState();
    calculateOutput();
}

DoubleT Lloyd2003Muscle::calculateFiberVelocityFromFiberLength(DoubleT previousFiberLength,
    DoubleT currentFiberLength,
    DoubleT dt) {
    return (currentFiberLength - previousFiberLength) / dt;
}

DoubleT Lloyd2003Muscle::integrateFiberLength(DoubleT dt) {
    DoubleT minFiberLength = 0.2 * p_.optimalFiberLength;
    DoubleT maxFiberLength = 2 * p_.optimalFiberLength;

    auto f = [&](double fiberLength) -> double {
        DoubleT fiberVelocity = Lloyd2003Muscle::calculateFiberVelocityFromFiberLength(
            this->s_.fiberLength, fiberLength, dt);
        DoubleT tendonForce = Lloyd2003Muscle::calculateTendonForce(
            this->i_.musculotendonLength, fiberLength, this->p_);
        DoubleT fiberForce = Lloyd2003Muscle::calculateFiberForce(
            this->i_.activation, fiberLength, fiberVelocity, this->p_);
        DoubleT v = tendonForce - fiberForce;
        return v;
    };

    double currentFiberLength = wdbSolve(f, minFiberLength, maxFiberLength, 1e-9);
    return currentFiberLength;
}

DoubleT Lloyd2003Muscle::calculateOptimalFiberLengthAtT(DoubleT activation,
    const Lloyd2003Muscle::Parameters &p) {
    return p.optimalFiberLength * (p.percentageChange * (1.0 - activation) + 1);
}

DoubleT Lloyd2003Muscle::calculateNormalisedFiberLengthAtT(DoubleT activation,
    DoubleT fiberLength,
    const Lloyd2003Muscle::Parameters &p) {
    return fiberLength / calculateOptimalFiberLengthAtT(activation, p);
}

DoubleT Lloyd2003Muscle::calculateNormalisedFiberLength(DoubleT fiberLength,
    const Lloyd2003Muscle::Parameters &p) {
    return fiberLength / p.optimalFiberLength;
}

DoubleT Lloyd2003Muscle::calculateNormalisedFiberVelocity(DoubleT fiberVelocity,
    const Lloyd2003Muscle::Parameters &p) {
    DoubleT normalisedFiberVelocity =
        fiberVelocity / (p.optimalFiberLength * p.maxContractionVelocity);
    if (normalisedFiberVelocity > 1) { normalisedFiberVelocity = 1; }
    if (normalisedFiberVelocity < -1) { normalisedFiberVelocity = -1; }
    return normalisedFiberVelocity;
}

DoubleT Lloyd2003Muscle::calculatePennationAngle(DoubleT fiberLength,
    const Lloyd2003Muscle::Parameters &p) {
    DoubleT value{ p.optimalFiberLength * sin(p.pennationAngleAtOptimalFiberLength) / fiberLength };
    if (value < 0) { value = 0; }
    if (value > 0.99) { value = 0.99; }
    return asin(value);
}

DoubleT Lloyd2003Muscle::calculateTendonLength(DoubleT musculotendonLength,
    DoubleT fiberLength,
    const Lloyd2003Muscle::Parameters &p) {
    return musculotendonLength - fiberLength * cos(calculatePennationAngle(fiberLength, p));
}

DoubleT Lloyd2003Muscle::calculateTendonStrain(DoubleT musculotendonLength,
    DoubleT fiberLength,
    const Lloyd2003Muscle::Parameters &p) {
    DoubleT tendonStrain =
        (calculateTendonLength(musculotendonLength, fiberLength, p) - p.tendonSlackLength)
        / p.tendonSlackLength;
    if (tendonStrain < 0.) { tendonStrain = 0; }
    return tendonStrain;
}

DoubleT Lloyd2003Muscle::calculateFiberForce(DoubleT activation,
    DoubleT fiberLength,
    DoubleT fiberVelocity,
    const Lloyd2003Muscle::Parameters &p) {
    DoubleT normalisedFiberLengthAtT{ calculateNormalisedFiberLengthAtT(
        activation, fiberLength, p) };
    DoubleT normalisedFiberLength{ calculateNormalisedFiberLength(fiberLength, p) };
    DoubleT normalisedFiberVelocity{ calculateNormalisedFiberVelocity(fiberVelocity, p) };
    DoubleT pennationAngle{ calculatePennationAngle(fiberLength, p) };
    DoubleT fa{ p.activeForceLengthCurve.getValue(normalisedFiberLengthAtT) };
    DoubleT fp{ p.passiveForceLengthCurve.getValue(normalisedFiberLength) };
    DoubleT fv{ p.forceVelocityCurve.getValue(normalisedFiberVelocity) };

    DoubleT fiberForce{ p.maxIsometricForce * p.strengthCoefficient
                        * (fa * fv * activation + fp + p.damping * normalisedFiberVelocity)
                        * cos(pennationAngle) };
    // clamp muscle force
    if (fiberForce < 0.) { fiberForce = 0; }
    return fiberForce;
}


DoubleT Lloyd2003Muscle::calculateTendonForce(DoubleT musculotendonLength,
    DoubleT fiberLength,
    const Lloyd2003Muscle::Parameters &p) {
    DoubleT tendonStrain{ calculateTendonStrain(musculotendonLength, fiberLength, p) };
    DoubleT tendonForce{ p.strengthCoefficient * p.maxIsometricForce
                         * p.tendonForceStrainCurve.getValue(tendonStrain) };
    return tendonForce;
}


//This can be made generic when Concepts are available
void connectSocket(const ExponentialActivation &parent, Lloyd2003Muscle &child) {
    child.setInput(parent.getOutput<Activation>());
}

void connectSocket(const Activation &parent, Lloyd2003Muscle &child) {
    child.setInput(parent);
}

void connectSocket(const MusculotendonLength &parent, Lloyd2003Muscle &child) {
    child.setInput(parent);
}

}// namespace ceinms

#endif