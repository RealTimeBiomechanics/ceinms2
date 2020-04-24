#ifndef ceinms2_Lloyd2019Muscle_h
#define ceinms2_Lloyd2019Muscle_h

#include "ceinms2/Curve.h"
#include "ceinms2/WDBsolver.h"
#include <cmath>
using DoubleT = double;
using CurveOffline = ceinms::Curve<ceinms::CurveMode::Offline,
    ceinms::CurveMode::Interpolation::Cubic>;
namespace ceinms {

class Lloyd2019Muscle {
  public:
    using type = nullptr_t;
    static constexpr std::string_view class_name = "Lloyd2019Muscle";
    struct Parameters {
        Parameters()
            : optimalFiberLength(1.)
            , pennationAngleAtOptimalFiberLength(0.)
            , tendonSlackLength(0.8)
            , maxContractionVelocity(5.)
            , damping(0.1)
            , maxIsometricForce(100.)
            , strengthCoefficient(1.)
            , percentageChange(0.15) {}


        DoubleT optimalFiberLength;
        DoubleT pennationAngleAtOptimalFiberLength;
        DoubleT tendonSlackLength;
        DoubleT maxContractionVelocity;
        DoubleT damping;
        DoubleT maxIsometricForce;
        DoubleT strengthCoefficient;
        DoubleT percentageChange;
        CurveOffline forceVelocityCurve;
        CurveOffline activeForceLengthCurve;
        CurveOffline passiveForceLengthCurve;
        CurveOffline tendonForceStrainCurve;
    };

    struct Input {
        Input()
            : activation(0)
            , musculotendonLength(0) {}

        DoubleT activation;
        DoubleT musculotendonLength;
    };

    struct State {
        State()
            : fiberLength(0)
            , fiberVelocity(0) {}

        DoubleT fiberLength;
        DoubleT fiberVelocity;
    };

    struct Output {
        Output()
            : optimalFiberLengthAtT(0)
            , pennationAngle(0)
            , fiberForce(0)
            , activeForce(0)
            , passiveForce(0)
            , dampingForce(0)
            , normalisedFiberLength(0)
            , normalisedFiberLengthAtT(0)
            , normalisedFiberVelocity(0)
            , tendonLength(0)
            , tendonStrain(0)
            , tendonForce(0)
            , fa(0)
            , fv(0)
            , fp(0) {}

        DoubleT optimalFiberLengthAtT;
        DoubleT pennationAngle;
        DoubleT fiberForce;
        DoubleT activeForce;
        DoubleT passiveForce;
        DoubleT dampingForce;
        DoubleT normalisedFiberLength;
        DoubleT normalisedFiberLengthAtT;
        DoubleT normalisedFiberVelocity;
        DoubleT tendonLength;
        DoubleT tendonStrain;
        DoubleT tendonForce;
        DoubleT fa;
        DoubleT fv;
        DoubleT fp;
        DoubleT getPrimary() const { return fiberForce; }
    };

    Lloyd2019Muscle(Parameters parameters)
        : p_(parameters) {
        s_.fiberLength = p_.optimalFiberLength;
        s_.fiberVelocity = 0.;
    }

    void setActivation(DoubleT activation);
    void setMusculotendonLength(DoubleT musculotendonLength);
    void setInput(Input input);
    void equilibrate();
    void setState(State state);
    // From input and current state calculate the new state of the system
    void integrate(DoubleT dt);
    // The temporary state calculated via `integrate` becomes the new state
    void validateState();
    // from the internal state of the system and the input, calculate all the
    // output;
    void calculateOutput();
    std::string getName() const { return name_; }
    void setName(std::string name) { name_ = name; }
    // Convenience function that, from input and current state, calculate the
    // new state and all the output
    void evaluate(DoubleT dt);

    Parameters &updParameters() { return p_; }
    Parameters getParameters() { return p_; }
    State &updState() { return s_; }
    State getState() { return s_; }
    Output getOutput() { return o_; }

    /* calculateXYZ functions perform calculations and return the calculated
     * value*/
    static DoubleT calculateFiberVelocityFromFiberLength(
        DoubleT previousFiberLength,
        DoubleT currentFiberLength,
        DoubleT dt);
    static DoubleT calculateFiberForce(DoubleT activation,
        DoubleT fiberLength,
        DoubleT fiberVelocity,
        const Parameters &p);
    static DoubleT calculateTendonForce(DoubleT musculotendonLength,
        DoubleT fiberLength,
        const Parameters &p);
    static DoubleT calculatePennationAngle(DoubleT fiberLength,
        const Parameters &p);
    static DoubleT calculateOptimalFiberLengthAtT(DoubleT activation,
        const Parameters &p);
    static DoubleT calculateNormalisedFiberLengthAtT(DoubleT activation,
        DoubleT fiberLength,
        const Parameters &p);
    static DoubleT calculateNormalisedFiberLength(DoubleT fiberLength,
        const Parameters &p);
    static DoubleT calculateNormalisedFiberVelocity(DoubleT fiberVelocity,
        const Parameters &p);
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

void Lloyd2019Muscle::setActivation(DoubleT activation) {
    i_.activation = activation;
}

void Lloyd2019Muscle::setMusculotendonLength(DoubleT musculotendonLength) {
    i_.musculotendonLength = musculotendonLength;
}


void Lloyd2019Muscle::setInput(Input input) { i_ = input; }


void Lloyd2019Muscle::equilibrate() {
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

void Lloyd2019Muscle::integrate(DoubleT dt) {
    sNew_.fiberLength = integrateFiberLength(dt);
    sNew_.fiberVelocity = calculateFiberVelocityFromFiberLength(
        s_.fiberLength, sNew_.fiberLength, dt);
}

void Lloyd2019Muscle::validateState() { s_ = sNew_; }

void Lloyd2019Muscle::calculateOutput() {
    o_.normalisedFiberLengthAtT =
        calculateNormalisedFiberLengthAtT(i_.activation, s_.fiberLength, p_);
    o_.normalisedFiberLength =
        calculateNormalisedFiberLength(s_.fiberLength, p_);
    o_.normalisedFiberVelocity =
        calculateNormalisedFiberVelocity(s_.fiberVelocity, p_);
    o_.pennationAngle = calculatePennationAngle(s_.fiberLength, p_);
    o_.fa = p_.activeForceLengthCurve.getValue(o_.normalisedFiberLengthAtT);
    o_.fp = p_.passiveForceLengthCurve.getValue(o_.normalisedFiberLength);
    o_.fv = p_.forceVelocityCurve.getValue(o_.normalisedFiberVelocity);
    o_.activeForce =
        p_.maxIsometricForce * p_.strengthCoefficient * o_.fa * o_.fv;
    o_.passiveForce = p_.maxIsometricForce * p_.strengthCoefficient * o_.fp;
    o_.dampingForce = p_.maxIsometricForce * p_.strengthCoefficient * p_.damping
                      * o_.normalisedFiberVelocity;
    o_.fiberForce = calculateFiberForce(
        i_.activation, s_.fiberLength, s_.fiberVelocity, p_);
    o_.tendonStrain =
        calculateTendonStrain(i_.musculotendonLength, s_.fiberLength, p_);
    o_.tendonLength =
        calculateTendonLength(i_.musculotendonLength, s_.fiberLength, p_);
    o_.tendonForce =
        calculateTendonForce(i_.musculotendonLength, s_.fiberLength, p_);
}

void Lloyd2019Muscle::evaluate(DoubleT dt) {
    integrate(dt);
    validateState();
    calculateOutput();
}

DoubleT Lloyd2019Muscle::calculateFiberVelocityFromFiberLength(
    DoubleT previousFiberLength,
    DoubleT currentFiberLength,
    DoubleT dt) {
    return (currentFiberLength - previousFiberLength) / dt;
}

DoubleT Lloyd2019Muscle::integrateFiberLength(DoubleT dt) {
    DoubleT minFiberLength = 0.2 * p_.optimalFiberLength;
    DoubleT maxFiberLength = 2 * p_.optimalFiberLength;

    auto f = [&](double fiberLength) -> double {
        DoubleT fiberVelocity =
            Lloyd2019Muscle::calculateFiberVelocityFromFiberLength(
                this->s_.fiberLength, fiberLength, dt);
        DoubleT tendonForce = Lloyd2019Muscle::calculateTendonForce(
            this->i_.musculotendonLength, fiberLength, this->p_);
        DoubleT fiberForce = Lloyd2019Muscle::calculateFiberForce(
            this->i_.activation, fiberLength, fiberVelocity, this->p_);
        DoubleT v = tendonForce - fiberForce;
        return v;
    };

    double currentFiberLength =
        wdbSolve(f, minFiberLength, maxFiberLength, 1e-9);
    return currentFiberLength;
}

DoubleT Lloyd2019Muscle::calculateOptimalFiberLengthAtT(DoubleT activation,
    const Lloyd2019Muscle::Parameters &p) {
    return p.optimalFiberLength * (p.percentageChange * (1.0 - activation) + 1);
}

DoubleT Lloyd2019Muscle::calculateNormalisedFiberLengthAtT(DoubleT activation,
    DoubleT fiberLength,
    const Lloyd2019Muscle::Parameters &p) {
    return fiberLength / calculateOptimalFiberLengthAtT(activation, p);
}

DoubleT Lloyd2019Muscle::calculateNormalisedFiberLength(DoubleT fiberLength,
    const Lloyd2019Muscle::Parameters &p) {
    return fiberLength / p.optimalFiberLength;
}

DoubleT Lloyd2019Muscle::calculateNormalisedFiberVelocity(DoubleT fiberVelocity,
    const Lloyd2019Muscle::Parameters &p) {
    DoubleT normalisedFiberVelocity =
        fiberVelocity / (p.optimalFiberLength * p.maxContractionVelocity);
    if (normalisedFiberVelocity > 1) normalisedFiberVelocity = 1;
    if (normalisedFiberVelocity < -1) normalisedFiberVelocity = -1;
    return normalisedFiberVelocity;
}

DoubleT Lloyd2019Muscle::calculatePennationAngle(DoubleT fiberLength,
    const Lloyd2019Muscle::Parameters &p) {
    DoubleT value{ p.optimalFiberLength
                   * sin(p.pennationAngleAtOptimalFiberLength) / fiberLength };
    if (value < 0) value = 0;
    if (value > 0.99) value = 0.99;
    return asin(value);
}

DoubleT Lloyd2019Muscle::calculateTendonLength(DoubleT musculotendonLength,
    DoubleT fiberLength,
    const Lloyd2019Muscle::Parameters &p) {
    return musculotendonLength
           - fiberLength * cos(calculatePennationAngle(fiberLength, p));
}

DoubleT Lloyd2019Muscle::calculateTendonStrain(DoubleT musculotendonLength,
    DoubleT fiberLength,
    const Lloyd2019Muscle::Parameters &p) {
    DoubleT tendonStrain =
        (calculateTendonLength(musculotendonLength, fiberLength, p)
            - p.tendonSlackLength)
        / p.tendonSlackLength;
    if (tendonStrain < 0.) tendonStrain = 0;
    return tendonStrain;
}

DoubleT Lloyd2019Muscle::calculateFiberForce(DoubleT activation,
    DoubleT fiberLength,
    DoubleT fiberVelocity,
    const Lloyd2019Muscle::Parameters &p) {
    DoubleT normalisedFiberLengthAtT{ calculateNormalisedFiberLengthAtT(
        activation, fiberLength, p) };
    DoubleT normalisedFiberLength{ calculateNormalisedFiberLength(
        fiberLength, p) };
    DoubleT normalisedFiberVelocity{ calculateNormalisedFiberVelocity(
        fiberVelocity, p) };
    DoubleT pennationAngle{ calculatePennationAngle(fiberLength, p) };
    DoubleT fa{ p.activeForceLengthCurve.getValue(normalisedFiberLengthAtT) };
    DoubleT fp{ p.passiveForceLengthCurve.getValue(normalisedFiberLength) };
    DoubleT fv{ p.forceVelocityCurve.getValue(normalisedFiberVelocity) };

    DoubleT fiberForce{ p.maxIsometricForce * p.strengthCoefficient
                        * (fa * fv * activation + fp
                            + p.damping * normalisedFiberVelocity)
                        * cos(pennationAngle) };
    // clamp muscle force
    if (fiberForce < 0.) fiberForce = 0;
    return fiberForce;
}


DoubleT Lloyd2019Muscle::calculateTendonForce(DoubleT musculotendonLength,
    DoubleT fiberLength,
    const Lloyd2019Muscle::Parameters &p) {
    DoubleT tendonStrain{ calculateTendonStrain(
        musculotendonLength, fiberLength, p) };
    DoubleT tendonForce{ p.strengthCoefficient * p.maxIsometricForce
                         * p.tendonForceStrainCurve.getValue(tendonStrain) };
    return tendonForce;
}
}// namespace ceinms

#endif