#ifndef ceinms2_Lloyd2003Muscle_h
#define ceinms2_Lloyd2003Muscle_h
#include <cmath>
#include <string_view>
#include <boost/math/differentiation/autodiff.hpp>
#include "ceinms2/Types.h"
#include "ceinms2/Curve.h"
#include "ceinms2/WDBsolver.h"
#include "ceinms2/ExponentialActivation.h"

namespace ceinms {

class CubicSpline {
  public:
    struct Parameters {
        std::vector<DoubleT> x;
        std::vector<DoubleT> y;
    };
    CubicSpline() = default;
    CubicSpline(const Parameters &p) { set(p); }
    CubicSpline(const std::vector<DoubleT> &x, const std::vector<DoubleT> &y) { set({ x, y }); }

    DoubleT get(DoubleT x) const {
        x = std::min(x, spline_.getMaxX());
        x = std::max(x, spline_.getMinX());
        return spline_.getValue(x);
    }

    DoubleT getFirstDerivative(DoubleT x) const {
        x = std::min(x, spline_.getMaxX());
        x = std::max(x, spline_.getMinX());
        return spline_.getFirstDerivative(x);
    }

    template<typename T>
    auto get(const T &cr) const {
        using root_type = typename T::root_type;
        constexpr size_t order = T::order_sum;
        root_type const d0 = get(static_cast<root_type>(cr));
        if constexpr (order == 0)
            return fvar(d0);
        else {
            const root_type d1 = getFirstDerivative(static_cast<root_type>(cr));
            const root_type derivatives[2]{ d0, d1 };
            return cr.apply_derivatives(
                order, [&derivatives](size_t i) { return derivatives[i & 1]; });
        }
    }

    void set(const Parameters &p) { spline_.resetPointsWith(p.x, p.y); }

  private:
    Curve<CurveMode::Offline> spline_;
};


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


class Lloyd2003Muscle {
  public:
    using type = nullptr_t;
    using concept_t = component_t;
    static constexpr std::string_view class_name = "Lloyd2003Muscle";
    struct Parameters {
        DoubleT optimalFiberLength{ 1. };
        DoubleT pennationAngleAtOptimalFiberLength{ 0. };
        DoubleT tendonSlackLength{ 0.8 };
        DoubleT maxContractionVelocity{ 5. };
        DoubleT damping{ 0.1 };
        DoubleT maxIsometricForce{ 100. };
        DoubleT strengthCoefficient{ 1. };
        DoubleT percentageChange{ 0.15 };
        CubicSpline activeForceLengthCurve;
        CubicSpline passiveForceLengthCurve;
        CubicSpline forceVelocityCurve;
        CubicSpline tendonForceStrainCurve;

        Parameters() {
            activeForceLengthCurve = getDefaultActiveForceLengthCurve();
            passiveForceLengthCurve = getDefaultPassiveForceLengthCurve();
            forceVelocityCurve = getDefaultForceVelocityCurve();
            tendonForceStrainCurve = getDefaultTendonForceStrainCurve();
        }
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
    template<typename W, typename X, typename T>
    static auto calculateFiberVelocityFromFiberLength(const W& previousFiberLength,
        const X& currentFiberLength,
        const T& dt);

    template<typename W, typename X, typename Y>
    static auto calculateFiberForce(const W& activation,
        const X& fiberLength,
        const Y& fiberVelocity,
        const Parameters &p);

    template<typename W, typename X, typename Y>
    static auto calculateFiberStiffness(const W& activation,
        const X& fiberLength,
        const Y& fiberVelocity,
        const Parameters &p);

  
    template<typename W>
    static auto calculatePennationAngle(const W& fiberLength, const Parameters &p);

    template<typename W>
    static auto calculateOptimalFiberLengthAtT(const W& activation, const Parameters &p);

    template<typename W, typename X>
    static auto calculateNormalisedFiberLengthAtT(const W& activation,
        const X& fiberLength,
        const Parameters &p);

    template<typename W>
    static auto calculateNormalisedFiberLength(const W& fiberLength, const Parameters &p);

    template<typename W>
    static auto calculateNormalisedFiberVelocity(const W& fiberVelocity, const Parameters &p);

    template<typename W, typename X>
    static auto calculateTendonLength(const W& musculotendonLength,
        const X& fiberLength,
        const Parameters &p);

    template<typename W, typename X>
    static auto calculateTendonStrain(const W& musculotendonLength,
        const X& fiberLength,
        const Parameters &p);

    template<typename W, typename X>
    static auto calculateTendonForce(const W &musculotendonLength,
        const X &fiberLength,
        const Parameters &p);

    static DoubleT calculateTendonStiffness(DoubleT musculotendonLength,
        const DoubleT fiberLength,
        const Parameters &p);

  private:

    DoubleT integrateFiberLength(DoubleT dt);
    std::string name_;
    Parameters p_;
    State s_, sNew_;
    Input i_;
    Output o_;
};

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
    s_ = state;
}


template<typename W, typename X, typename T>
auto Lloyd2003Muscle::calculateFiberVelocityFromFiberLength(
    const W& previousFiberLength,
    const X& currentFiberLength,
    const T& dt) {
    return (currentFiberLength - previousFiberLength) / dt;
}

template<typename W>
auto Lloyd2003Muscle::calculateOptimalFiberLengthAtT(const W& activation,
    const Parameters &p) {
    return p.optimalFiberLength * (p.percentageChange * (1.0 - activation) + 1);
}

template<typename W, typename X>
auto Lloyd2003Muscle::calculateNormalisedFiberLengthAtT(const W& activation,
    const X& fiberLength,
    const Parameters &p) {
    return fiberLength / calculateOptimalFiberLengthAtT(activation, p);
}

template<typename W>
auto Lloyd2003Muscle::calculateNormalisedFiberLength(const W& fiberLength,
    const Parameters &p) {
    return fiberLength / p.optimalFiberLength;
}

template<typename W>
auto Lloyd2003Muscle::calculateNormalisedFiberVelocity(const W& fiberVelocity,
    const Parameters &p) {
    auto normalisedFiberVelocity =
        fiberVelocity / (p.optimalFiberLength * p.maxContractionVelocity);
    if (normalisedFiberVelocity > 1) { normalisedFiberVelocity = 1; }
    if (normalisedFiberVelocity < -1) { normalisedFiberVelocity = -1; }
    return normalisedFiberVelocity;
}

template<typename W>
auto Lloyd2003Muscle::calculatePennationAngle(const W& fiberLength,
    const Parameters &p) {
    auto value{ p.optimalFiberLength * sin(p.pennationAngleAtOptimalFiberLength) / fiberLength };
    if (value < 0) { value = 0; }
    if (value > 0.99) { value = 0.99; }
    return asin(value);
}

template<typename W, typename X>
auto Lloyd2003Muscle::calculateTendonLength(const W& musculotendonLength,
    const X& fiberLength,
    const Parameters &p) {
    return musculotendonLength - fiberLength * cos(calculatePennationAngle(fiberLength, p));
}

template<typename W, typename X>
auto Lloyd2003Muscle::calculateTendonStrain(const W& musculotendonLength,
    const X& fiberLength,
    const Parameters &p) {
    auto tendonStrain =
        (calculateTendonLength(musculotendonLength, fiberLength, p) - p.tendonSlackLength)
        / p.tendonSlackLength;
    if (tendonStrain < 0.) { tendonStrain = 0; }
    return tendonStrain;
}

template<typename W, typename X, typename Y>
auto Lloyd2003Muscle::calculateFiberForce(const W& activation,
    const X& fiberLength,
    const Y& fiberVelocity,
    const Parameters &p) {
    const auto normalisedFiberLengthAtT{ calculateNormalisedFiberLengthAtT(
        activation, fiberLength, p) };
    const auto normalisedFiberLength{ calculateNormalisedFiberLength(fiberLength, p) };
    const auto normalisedFiberVelocity{ calculateNormalisedFiberVelocity(fiberVelocity, p) };
    const auto pennationAngle{ calculatePennationAngle(fiberLength, p) };
    const auto fa{ p.activeForceLengthCurve.get(normalisedFiberLengthAtT) };
    const auto fp{ p.passiveForceLengthCurve.get(normalisedFiberLength) };
    const auto fv{ p.forceVelocityCurve.get(normalisedFiberVelocity) };

    auto fiberForce{ p.maxIsometricForce * p.strengthCoefficient
                        * (fa * fv * activation + fp + p.damping * normalisedFiberVelocity)
                        * cos(pennationAngle) };
    // clamp muscle force
    if (fiberForce < 0.) { fiberForce = 0; }
    return fiberForce;
}

template<typename W, typename X, typename Y>
static auto Lloyd2003Muscle::calculateFiberStiffness(const W& activation,
    const X& fiberLength,
    const Y& fiberVelocity,
    const Parameters &p) {
    return 0.
}

template<typename W, typename X>
auto Lloyd2003Muscle::calculateTendonForce(const W& musculotendonLength,
    const X& fiberLength,
    const Parameters &p) {
    const auto tendonStrain{ calculateTendonStrain(musculotendonLength, fiberLength, p) };
    const auto tendonForce{ p.strengthCoefficient * p.maxIsometricForce
                         * p.tendonForceStrainCurve.get(tendonStrain) };
    return tendonForce;
}

DoubleT Lloyd2003Muscle::calculateTendonStiffness(DoubleT musculotendonLength,
    const DoubleT fiberLength,
    const Parameters &p) {
    const auto tendonStrain{ calculateTendonStrain(musculotendonLength, fiberLength, p) };
    const auto tendonStiffness{ p.tendonForceStrainCurve.getFirstDerivative(tendonStrain) };
    return tendonStiffness;
}


DoubleT Lloyd2003Muscle::integrateFiberLength(DoubleT dt) {
    DoubleT minFiberLength = 0.2 * p_.optimalFiberLength;
    DoubleT maxFiberLength = 2 * p_.optimalFiberLength;

    auto f = [&](double fiberLength) -> double {
        const auto fiberVelocity = Lloyd2003Muscle::calculateFiberVelocityFromFiberLength(
            this->s_.fiberLength, fiberLength, dt);
        const auto tendonForce = Lloyd2003Muscle::calculateTendonForce(
            this->i_.musculotendonLength, fiberLength, this->p_);
        const auto fiberForce = Lloyd2003Muscle::calculateFiberForce(
            this->i_.activation, fiberLength, fiberVelocity, this->p_);
        const auto v = tendonForce - fiberForce;
        return v;
    };

    const auto currentFiberLength = wdbSolve(f, minFiberLength, maxFiberLength, 1e-9);
    return currentFiberLength;
}


void Lloyd2003Muscle::calculateOutput() {
    o_.normalisedFiberLengthAtT =
        calculateNormalisedFiberLengthAtT(i_.activation, s_.fiberLength, p_);
    o_.normalisedFiberLength = calculateNormalisedFiberLength(s_.fiberLength, p_);
    o_.normalisedFiberVelocity = calculateNormalisedFiberVelocity(s_.fiberVelocity, p_);
    o_.pennationAngle = calculatePennationAngle(s_.fiberLength, p_);
    o_.fa = p_.activeForceLengthCurve.get(o_.normalisedFiberLengthAtT);
    o_.fp = p_.passiveForceLengthCurve.get(o_.normalisedFiberLength);
    o_.fv = p_.forceVelocityCurve.get(o_.normalisedFiberVelocity);
    o_.activeForce = p_.maxIsometricForce * p_.strengthCoefficient * o_.fa * o_.fv;
    o_.passiveForce = p_.maxIsometricForce * p_.strengthCoefficient * o_.fp;
    o_.dampingForce =
        p_.maxIsometricForce * p_.strengthCoefficient * p_.damping * o_.normalisedFiberVelocity;
    o_.fiberForce = calculateFiberForce(i_.activation, s_.fiberLength, s_.fiberVelocity, p_);
    o_.tendonStrain = calculateTendonStrain(i_.musculotendonLength, s_.fiberLength, p_);
    o_.tendonLength = calculateTendonLength(i_.musculotendonLength, s_.fiberLength, p_);
    o_.tendonForce = calculateTendonForce(i_.musculotendonLength, s_.fiberLength, p_);
}
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

void Lloyd2003Muscle::validateState() {
    s_ = sNew_;
}

void Lloyd2003Muscle::evaluate(DoubleT dt) {
    integrate(dt);
    validateState();
    calculateOutput();
}
// This can be made generic when Concepts are available
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