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

CubicSpline getDefaultTendonForceStrainCurve();
CubicSpline getDefaultActiveForceLengthCurve();
CubicSpline getDefaultPassiveForceLengthCurve();
CubicSpline getDefaultForceVelocityCurve();


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


template<ActivationGenerator T>
void connectSocket(const T &parent, Lloyd2003Muscle &child) {
    child.setActivation(parent.getOutput().activation);
}

void connectSocket(const Activation &parent, Lloyd2003Muscle &child);
void connectSocket(const MusculotendonLength &parent, Lloyd2003Muscle &child);

}// namespace ceinms

#endif
