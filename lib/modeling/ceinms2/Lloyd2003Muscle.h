#ifndef ceinms2_Lloyd2003Muscle_h
#define ceinms2_Lloyd2003Muscle_h
#include <cmath>
#include <string_view>
#include <boost/math/differentiation/autodiff.hpp>
#include "ceinms2/Types.h"
#include "ceinms2/Curve.h"
#include "ceinms2/WDBsolver.h"
#include <Eigen/Dense>

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

    constexpr DoubleT get(DoubleT x) const {
        x = std::min(x, spline_.getMaxX());
        x = std::max(x, spline_.getMinX());
        return spline_.getValue(x);
    }

    constexpr DoubleT getFirstDerivative(DoubleT x) const {
        x = std::min(x, spline_.getMaxX());
        x = std::max(x, spline_.getMinX());
        return spline_.getFirstDerivative(x);
    }

     constexpr DoubleT getSecondDerivative(DoubleT x) const {
        x = std::min(x, spline_.getMaxX());
        x = std::max(x, spline_.getMinX());
        return spline_.getSecondDerivative(x);
    }

    template<typename T>
    constexpr auto get(const T &cr) const {
        using root_type = typename T::root_type;
        constexpr size_t order = T::order_sum;
        root_type const d0 = get(static_cast<root_type>(cr));
        if constexpr (order == 0)
            return fvar(d0);
        else {
            const root_type d1 = getFirstDerivative(static_cast<root_type>(cr));
            const root_type d2 = getSecondDerivative(static_cast<root_type>(cr));
            const root_type derivatives[3]{ d0, d1, d2 };
            return cr.apply_derivatives(
                order, [&derivatives](size_t i) { return derivatives[i]; });
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
        DoubleT force{ 0. };
        DoubleT activeForce{ 0. };
        DoubleT passiveForce{ 0. };
        DoubleT dampingForce{ 0. };
        DoubleT normalizedFiberLength{ 0. };
        DoubleT normalizedFiberLengthAtT{ 0. };
        DoubleT normalizedFiberVelocity{ 0. };
        DoubleT tendonLength{ 0. };
        DoubleT tendonStrain{ 0. };
        DoubleT tendonForce{ 0. };
        DoubleT musculotendonStiffness{ 0. };
        DoubleT fiberStiffness{ 0. };
        DoubleT tendonStiffness{ 0. };
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
    void setInput(Input input) { i_ = input; }
    void setParameters(Parameters parameters) { p_ = parameters; }
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
    void setName(std::string_view name) { name_ = name; }
    // Convenience function that, from input and current state, calculate the
    // new state and all the output
    void evaluate(DoubleT dt);

    [[nodiscard]] Parameters &getParameters() { return p_; }
    [[nodiscard]] const Parameters& getParameters() const { return p_; }
    [[nodiscard]] Input &getInput() { return i_; }
    [[nodiscard]] const Input &getInput() const { return i_; }
    [[nodiscard]] const State &getState() const { return s_; }
    [[nodiscard]] State &getState() { return s_; }
    [[nodiscard]] const Output &getOutput() const { return o_; }
    template<typename T, std::enable_if_t<std::is_same<T, Force>::value, int> = 0>
    [[nodiscard]] Force getOutput() const {
        return Force{ o_.force };
    }
   
    [[nodiscard]] DoubleT calculateFiberForceProjectedOnTendon() const;
    [[nodiscard]] DoubleT calculatePennationAngle() const;
    [[nodiscard]] DoubleT calculateNormalizedFiberVelocity() const;
    [[nodiscard]] DoubleT calculateNormalizedFiberLengthAtT() const;
    [[nodiscard]] DoubleT calculateNormalizedFiberLength() const;
    [[nodiscard]] DoubleT calculateFiberForce() const;
    [[nodiscard]] DoubleT calculateTendonStiffness() const;
    [[nodiscard]] DoubleT calculateFiberStiffness() const;
    [[nodiscard]] DoubleT calculateFiberStiffnessProjectedOnTendon() const;
    [[nodiscard]] DoubleT calculateMusculotendonStiffness() const;
    [[nodiscard]] DoubleT calculateTendonStrain() const;
    [[nodiscard]] DoubleT calculateTendonLength() const;
    [[nodiscard]] DoubleT calculateTendonForce() const;
    [[nodiscard]] auto calculateJacobian() const;
  private:

     template<typename I0,
        typename S0,
        typename S1,
        typename P0,
        typename P1,
        typename P2,
        typename P3,
        typename P4,
        typename P5,
        typename P6>
    constexpr static auto calculateFiberForceProjectedOnTendon(const I0 &activation,
        const S0 &fiberLength,
        const S1 &fiberVelocity,
        const P0 &optimalFiberLength,
        const P1 &pennationAngleAtOptimalFiberLength,
        const P2 &maxIsometricForce,
        const P3 &strengthCoefficient,
        const P4 &damping,
        const P5 &percentageChange,
        const P6 &maxContractionVelocity,
        const CubicSpline &activeForceLengthCurve,
        const CubicSpline &passiveForceLengthCurve,
        const CubicSpline &forceVelocityCurve);

   template<typename I0,
         typename S0,
         typename S1,
         typename P0,
         typename P1,
         typename P2,
         typename P3,
         typename P4,
         typename P5>
    static constexpr auto calculateFiberForce(
        const I0 &activation,
        const S0 &fiberLength,
        const S1 &fiberVelocity,
        const P0 &maxIsometricForce,
        const P1 &strengthCoefficient,
        const P2 &damping,
        const P3 &optimalFiberLength,
        const P4 &percentageChange,
        const P5 &maxContractionVelocity,
        const CubicSpline &activeForceLengthCurve,
        const CubicSpline &passiveForceLengthCurve,
        const CubicSpline &forceVelocityCurve);
    
    template<typename S0, typename T>
    static constexpr auto calculateFiberVelocityFromFiberLength(const S0 &previousFiberLength,
        const S0 &currentFiberLength,
        const T &dt);

    template<typename I0, typename P0, typename P1>
    static constexpr auto calculateOptimalFiberLengthAtT(const I0 &activation,
        const P0 &optimalFiberLength,
        const P1 &percentageChange);

    template<typename S0, typename P0, typename P1>
    static constexpr auto calculatePennationAngle(const S0 &fiberLength,
        const P0 &optimalFiberLength,
        const P1 &pennationAngleAtOptimalFiberLength);

    template<typename S0, typename P0, typename P1>
    static constexpr auto calculateNormalizedFiberVelocity(
        const S0 &fiberVelocity,
        P0 &optimalFiberLength,
        P1 &maxContractionVelocity);

    template<typename I0, typename S0, typename P0, typename P1>
    static constexpr auto calculateTendonLength(const I0 &musculotendonLength,
        const S0 &fiberLength,
        const P0 &optimalFiberLength,
        const P1 &pennationAngleAtOptimalFiberLength);

    template<typename I0, typename S0, typename P0, typename P1>
    static constexpr auto calculateNormalizedFiberLengthAtT(const I0 &activation,
        const S0 &fiberLength,
        const P0 &optimalFiberLength,
        const P1 &percentageChange);

    template<typename S0, typename P0>
    static constexpr auto calculateNormalizedFiberLength(const S0 &fiberLength,
        const P0 &optimalFiberLength);
    
    template<typename I0, typename S0, typename P0, typename P1, typename P2>
    static constexpr auto calculateTendonStrain(const I0 &musculotendonLength,
        const S0 &fiberLength,
        const P0 &optimalFiberLength,
        const P1 &pennationAngleAtOptimalFiberLength,
        const P2 &tendonSlackLength);

    template<typename I0, typename S0, typename P0, typename P1, typename P2, typename P3, typename P4>
    static constexpr auto calculateTendonForce(const I0 &musculotendonLength,
    const S0 &fiberLength,
    const P0 &optimalFiberLength,
    const P1 &pennationAngleAtOptimalFiberLength,
    const P2 &tendonSlackLength,
    const P3 &strengthCoefficient,
    const P4 &maxIsometricForce,
    const CubicSpline &tendonForceStrainCurve);

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
    sNew_ = state;
}


  template<typename I0,
    typename S0,
    typename S1,
    typename P0,
    typename P1,
    typename P2,
    typename P3,
    typename P4,
    typename P5,
    typename P6>
constexpr auto Lloyd2003Muscle::calculateFiberForceProjectedOnTendon(const I0 &activation,
    const S0 &fiberLength,
    const S1 &fiberVelocity,
    const P0 &optimalFiberLength,
    const P1 &pennationAngleAtOptimalFiberLength,
    const P2 &maxIsometricForce,
    const P3 &strengthCoefficient,
    const P4 &damping,
    const P5 &percentageChange,
    const P6 &maxContractionVelocity,
    const CubicSpline &activeForceLengthCurve,
    const CubicSpline &passiveForceLengthCurve,
    const CubicSpline &forceVelocityCurve) {

    const  auto pennationAngle{ calculatePennationAngle(
        fiberLength, optimalFiberLength, pennationAngleAtOptimalFiberLength) };
      return calculateFiberForce(activation,
                 fiberLength,
                 fiberVelocity,
                 maxIsometricForce,
                 strengthCoefficient,
                 damping,
                 optimalFiberLength,
                 percentageChange,
                 maxContractionVelocity,
                 activeForceLengthCurve,
                 passiveForceLengthCurve,
                 forceVelocityCurve)
             * cos(pennationAngle);
  }

template<typename I0,
    typename S0,
    typename S1,
    typename P0,
    typename P1,
    typename P2,
    typename P3,
    typename P4,
    typename P5>
constexpr auto Lloyd2003Muscle::calculateFiberForce(const I0 &activation,
    const S0 &fiberLength,
    const S1 &fiberVelocity,
    const P0 &maxIsometricForce,
    const P1 &strengthCoefficient,
    const P2 &damping,
    const P3 &optimalFiberLength,
    const P4 &percentageChange,
    const P5 &maxContractionVelocity,
    const CubicSpline &activeForceLengthCurve,
    const CubicSpline &passiveForceLengthCurve,
    const CubicSpline &forceVelocityCurve) {
    const auto normalizedFiberLengthAtT{ calculateNormalizedFiberLengthAtT(
        activation, fiberLength, optimalFiberLength, percentageChange) };
    const auto normalizedFiberLength{ calculateNormalizedFiberLength(
        fiberLength, optimalFiberLength) };
    const auto normalizedFiberVelocity{ calculateNormalizedFiberVelocity(fiberVelocity,optimalFiberLength, maxContractionVelocity) };
    const auto fa{ activeForceLengthCurve.get(normalizedFiberLengthAtT) };
    const auto fp{ passiveForceLengthCurve.get(normalizedFiberLength) };
    const auto fv{ forceVelocityCurve.get(normalizedFiberVelocity) };

    const auto fiberForce{ maxIsometricForce * strengthCoefficient
                           * (fa * fv * activation + fp + damping * normalizedFiberVelocity) };
    // clamp muscle force
    if (fiberForce < 0.) return .0 * fiberForce;
    return fiberForce;
}

template<typename S0, typename P0, typename P1>
constexpr auto Lloyd2003Muscle::calculatePennationAngle(const S0 &fiberLength,
    const P0 &optimalFiberLength,
    const P1 &pennationAngleAtOptimalFiberLength) {
    const auto value{ optimalFiberLength * sin(pennationAngleAtOptimalFiberLength) / fiberLength };
    if (value < 0) { return  asin(.0 *value); }
    if (value > 0.99) { return asin(.0*value + 0.99); }
    return asin(value);
}

template<typename S0, typename T>
constexpr auto Lloyd2003Muscle::calculateFiberVelocityFromFiberLength(
    const S0& previousFiberLength,
    const S0& currentFiberLength,
    const T& dt) {
    return (currentFiberLength - previousFiberLength) / dt;
}

template<typename I0, typename P0, typename P1>
constexpr auto Lloyd2003Muscle::calculateOptimalFiberLengthAtT(const I0 &activation,
    const P0 &optimalFiberLength,
    const P1& percentageChange) {
    return optimalFiberLength * (percentageChange * (1.0 - activation) + 1);
}

template<typename I0, typename S0, typename P0, typename P1>
constexpr auto Lloyd2003Muscle::calculateNormalizedFiberLengthAtT(const I0 &activation,
    const S0 &fiberLength,
    const P0 &optimalFiberLength,
    const P1 &percentageChange) {
    return fiberLength
           / calculateOptimalFiberLengthAtT(activation, optimalFiberLength, percentageChange);
}

template<typename S0, typename P0>
constexpr auto Lloyd2003Muscle::calculateNormalizedFiberLength(const S0 &fiberLength,
    const P0 &optimalFiberLength) {
    return fiberLength / optimalFiberLength;
}

template<typename S0, typename P0, typename P1>
constexpr auto Lloyd2003Muscle::calculateNormalizedFiberVelocity(const S0 &fiberVelocity,
    P0 &optimalFiberLength,
    P1 &maxContractionVelocity) {
    auto normalizedFiberVelocity =
        fiberVelocity / (optimalFiberLength * maxContractionVelocity);
    if (normalizedFiberVelocity > 1) { return 0. * normalizedFiberVelocity + 1; }
    if (normalizedFiberVelocity < -1) { return 0. * normalizedFiberVelocity - 1; }
    return normalizedFiberVelocity;
}

template<typename I0, typename S0, typename P0, typename P1>
constexpr auto Lloyd2003Muscle::calculateTendonLength(const I0 &musculotendonLength,
    const S0 &fiberLength,
    const P0 &optimalFiberLength,
    const P1 &pennationAngleAtOptimalFiberLength) {
    return musculotendonLength
           - fiberLength
                 * cos(calculatePennationAngle(
                     fiberLength, optimalFiberLength, pennationAngleAtOptimalFiberLength));
}


template<typename I0, typename S0, typename P0, typename P1, typename P2>
constexpr auto Lloyd2003Muscle::calculateTendonStrain(
    const I0 &musculotendonLength,
    const S0 &fiberLength,
    const P0 &optimalFiberLength,
    const P1 &pennationAngleAtOptimalFiberLength,
    const P2 &tendonSlackLength) {
    auto tendonStrain = (calculateTendonLength(musculotendonLength,
                             fiberLength,
                             optimalFiberLength,
                             pennationAngleAtOptimalFiberLength)
                            - tendonSlackLength)
                        / tendonSlackLength;
    if (tendonStrain < 0.) { return 0. * tendonStrain; }
    return tendonStrain;
}

template<typename I0, typename S0, typename P0, typename P1, typename P2, typename P3, typename P4>
constexpr auto Lloyd2003Muscle::calculateTendonForce(
    const I0 &musculotendonLength,
    const S0 &fiberLength,
    const P0 &optimalFiberLength,
    const P1 &pennationAngleAtOptimalFiberLength,
    const P2 &tendonSlackLength,
    const P3 &strengthCoefficient,
    const P4 &maxIsometricForce,
    const CubicSpline &tendonForceStrainCurve) {
    const auto tendonStrain{ calculateTendonStrain(musculotendonLength,
        fiberLength,
        optimalFiberLength,
        pennationAngleAtOptimalFiberLength,
        tendonSlackLength) };
    const auto tendonForce{ strengthCoefficient * maxIsometricForce
                            * tendonForceStrainCurve.get(tendonStrain) };
    return tendonForce;
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
    o_.normalizedFiberLengthAtT = calculateNormalizedFiberLengthAtT();
 
    o_.normalizedFiberLength = calculateNormalizedFiberLength();
    o_.normalizedFiberVelocity = calculateNormalizedFiberVelocity();
    o_.pennationAngle = calculatePennationAngle();
    o_.fa = p_.activeForceLengthCurve.get(o_.normalizedFiberLengthAtT);
    o_.fp = p_.passiveForceLengthCurve.get(o_.normalizedFiberLength);
    o_.fv = p_.forceVelocityCurve.get(o_.normalizedFiberVelocity);
    o_.activeForce = p_.maxIsometricForce * p_.strengthCoefficient * o_.fa * o_.fv * i_.activation;
    o_.passiveForce = p_.maxIsometricForce * p_.strengthCoefficient * o_.fp;
    o_.dampingForce =
        p_.maxIsometricForce * p_.strengthCoefficient * p_.damping * o_.normalizedFiberVelocity;
    o_.fiberForce = calculateFiberForce();
    o_.force = calculateFiberForceProjectedOnTendon();
    o_.tendonStrain = calculateTendonStrain();
    o_.tendonLength = calculateTendonLength();
    o_.tendonForce = calculateTendonForce();
    o_.musculotendonStiffness = calculateMusculotendonStiffness();
    o_.fiberStiffness = calculateFiberStiffness();
    o_.tendonStiffness = calculateTendonStiffness();
}

DoubleT Lloyd2003Muscle::calculateNormalizedFiberLengthAtT() const {
    return calculateNormalizedFiberLengthAtT(
        i_.activation, sNew_.fiberLength, p_.optimalFiberLength, p_.percentageChange);
}

DoubleT Lloyd2003Muscle::calculateNormalizedFiberLength() const {
    return calculateNormalizedFiberLength(sNew_.fiberLength, p_.optimalFiberLength);
}

DoubleT Lloyd2003Muscle::calculateNormalizedFiberVelocity() const {
    return calculateNormalizedFiberVelocity(sNew_.fiberVelocity, p_.optimalFiberLength, p_.maxContractionVelocity);
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
    return calculateTendonForce(
        i_.musculotendonLength,
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
    const auto force = calculateFiberForceProjectedOnTendon(
        i_.activation,
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
        make_ftuple<DoubleT, N, N, N, N, N, N, N, N, N, N>(
            i_.musculotendonLength,
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
            const auto aForce = calculateFiberForceProjectedOnTendon(
                aActivation,
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
            p_.forceVelocityCurve
            );

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

template<ActivationGenerator T>
void connectSocket(const T &parent, Lloyd2003Muscle &child) {
    child.setActivation(parent.getOutput().activation);
}

void connectSocket(const Activation &parent, Lloyd2003Muscle &child) {
    child.setInput(parent);
}

void connectSocket(const MusculotendonLength &parent, Lloyd2003Muscle &child) {
    child.setInput(parent);
}

}// namespace ceinms

#endif
