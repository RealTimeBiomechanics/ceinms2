#ifndef ceinms_Mileusnic2006MuscleSpindle_h
#define ceinms_Mileusnic2006MuscleSpindle_h
#include "ceinms2/Types.h"
#include <algorithm>

namespace ceinms {

/*This class implements a modified version of the model of a muscle spindle as 
presented in Mileusnic, M. P., I. E. Brown, N. Lan and G. E. Loeb (2006). 
"Mathematical models of proprioceptors. I. Control and transduction in the 
muscle spindle." J Neurophysiol 96(4): 1772-1788.

In this implementation we removed the mass of the intrafusal fiber and solved for 
the length of polar region using by solving the differential equation implicitly.
*/

class Mileusnic2006MuscleSpindle;

class Mileusnic2006IntrafusalFiber {
  public:
    using type = nullptr_t;
    using concept_t = component_t;
    static constexpr std::string_view class_name = "Mileusnic2006IntrafusalFiber";
    friend Mileusnic2006MuscleSpindle;
    struct State {
        DoubleT polarRegionLength = 1.;// normalized
        DoubleT polarRegionVelocity = 0.;// calculated from the normalized lengthPolarRegion
    };

    struct Parameters {
        DoubleT coefficientOfAsymmetryShortening;// `C_L`
        DoubleT coefficientOfAsymmetryLengthening;// `C_S`
        DoubleT
            muscleFascicleSlackLength;// `R` fascicle length below which force production is zero
        DoubleT sensoryRegionRestLength;// `Lsr`
        DoubleT sensoryRegionStiffness;// `Ksr
        DoubleT sensoryRegionThresholdLength;// `L_Nsr`
        DoubleT polarRegionRestLength;// `Lpr`
        DoubleT polarRegionStiffness;// `Kpr`
        DoubleT polarRegionThresholdLength;// `L_Npr`
        DoubleT velocityPowerTerm;// `a` in Mileusnic et al. 2006
        DoubleT beta0;// Passive damping coefficient
        DoubleT beta1;// Coef. of damping due to dyn. fusimotor input
        DoubleT beta2;// Coef. of damping due to stat. fusimotor input
        DoubleT gamma1;// Coef. of force generation due to dyn. fusimotor input
        DoubleT gamma2;// Coef. of force generation due to stat. fusimotor input
        DoubleT sensoryRegionStretchToPrimaryAfferentFiring;// `G`
        DoubleT sensoryRegionStretchToSecondaryAfferentFiring;// `G`
        DoubleT percentageSecondaryAfferentOnSensoryRegion;// `X`
        DoubleT secondaryAfferentRestLength;// `Lsecondary`
    };

    struct Input {
        DoubleT normalisedMuscleFiberLength = 1.;// `L` length of the muscle fiber
        DoubleT fStatic = 0.;// Static fusimotor activation, between 0 and 1
        DoubleT fDynamic = 0.;// Dynamic fusimotor activation, between 0 and 1
    };

    struct Output {
        DoubleT fiberTension = 0.;
        DoubleT beta = 0.;// polar region’s damping term
        DoubleT gamma = 0.;// active-state force generator term
        DoubleT primaryAfferent = 0.;
        DoubleT secondaryAfferent = 0.;
        DoubleT fspr = 0.;
        DoubleT fvpr = 0.;
    };

    struct Properties {};
    void setActivation(DoubleT fStatic, DoubleT fDynamic);
    void setNormalizedMuscleFiberLength(DoubleT normalizedMuscleFiberLength);
    void setInput(Activation fStatic, Activation fDynamic);
    void setInput(Input input);
    // From input and current state calculate the new state of the system
    void integrate(DoubleT dt);
    // The temporary state calculated via `integrate` becomes the new state
    void validateState();
    // from the internal state of the system and the input, calculate all the
    // output;
    void calculateOutput();
    void evaluate(DoubleT dt);
    [[nodiscard]] std::string getName() const { return name_; }
    void setName(std::string_view name) { name_ = name; }
    void equilibrate();

    [[nodiscard]] Parameters &getParameters() { return p_; }
    [[nodiscard]] const Parameters &getParameters() const { return p_; }
    [[nodiscard]] Input &getInput() { return i_; }
    [[nodiscard]] const Input &getInput() const { return i_; }
    [[nodiscard]] const State &getState() const { return s_; }
    [[nodiscard]] State &getState() { return s_; }
    [[nodiscard]] const Properties &getProperties() const { return properties_; }
    [[nodiscard]] Properties &getProperties() { return properties_; }
    [[nodiscard]] const Output &getOutput() const { return o_; }

    // Output is in pulses per second
    [[nodiscard]] DoubleT calculatePrimaryAfferent() const;

    // Output is in pulses per second
    [[nodiscard]] DoubleT calculateSecondaryAfferent() const;

    [[nodiscard]] DoubleT calculateFiberTension() const;

    [[nodiscard]] DoubleT calculateBeta() const;
    [[nodiscard]] DoubleT calculateGamma() const;
    [[nodiscard]] DoubleT calculateStiffnessContributionToPolarRegionForce() const;
    [[nodiscard]] DoubleT calculateVelocityContributionToPolarRegionForce() const;

  private:
    DoubleT integratePolarRegionLength(DoubleT dt);

    [[nodiscard]] static DoubleT calculatePolarRegionVelocityFromPolarRegionLength(
        DoubleT previousLengthPolarRegion,
        DoubleT currentLengthPolarRegion,
        DoubleT dt) {
        DoubleT polarRegionVelocity = (currentLengthPolarRegion - previousLengthPolarRegion) / dt;
        polarRegionVelocity = std::clamp(polarRegionVelocity, -10., 10.);
        return polarRegionVelocity;
    }

    [[nodiscard]] static DoubleT calculateVelocityContributionToPolarRegionForce(
        DoubleT polarRegionLength,
        DoubleT polarRegionVelocity,
        DoubleT beta,
        DoubleT coefficientOfAsymmetry,
        DoubleT muscleFascicleSlackLength,
        DoubleT velocityPowerTerm);

    [[nodiscard]] static DoubleT calculateStiffnessContributionToPolarRegionForce(
        DoubleT polarRegionLength,
        DoubleT polarRegionStiffness,
        DoubleT polarRegionRestLength);

    [[nodiscard]] static DoubleT calculateSensoryRegionForce(DoubleT normalisedMuscleFiberLength,
        DoubleT lengthPolarRegion,
        DoubleT stiffnessSensoryRegion,
        DoubleT restLengthSensoryRegion);

    [[nodiscard]] static DoubleT calculateBeta(DoubleT fDynamic,
        DoubleT fStatic,
        DoubleT beta0,
        DoubleT beta1,
        DoubleT beta2);

    [[nodiscard]] static DoubleT
        calculateGamma(DoubleT fDynamic, DoubleT fStatic, DoubleT gamma1, DoubleT gamma2);

    std::string name_;
    Parameters p_;
    State s_, sNew_;
    Input i_;
    Output o_;
    Properties properties_;
};

class Mileusnic2006IntrafusalFiberActivationDynamics {
  public:
    using type = nullptr_t;
    using concept_t = component_t;
    static constexpr std::string_view class_name = "Mileusnic2006IntrafusalFiberActivationDynamics";

    struct State {
        DoubleT activation = 0.;// `f`
    };

    struct Parameters {
        DoubleT tau = 1.;// Low-pass filter time constant
        DoubleT freq = 1.;// Constant relating the fusimotor frequency to activation
    };

    struct Input {
        DoubleT fusimotorFrequency = 0.;// pulses per second
    };

    struct Properties {
        bool isChain = false;
    };

    using Output = State;
    [[nodiscard]] Parameters &getParameters() { return p_; }
    [[nodiscard]] const Parameters &getParameters() const { return p_; }
    [[nodiscard]] Input &getInput() { return i_; }
    [[nodiscard]] const Input &getInput() const { return i_; }
    [[nodiscard]] const State &getState() const { return s_; }
    [[nodiscard]] State &getState() { return s_; }
    [[nodiscard]] const Output &getOutput() const { return s_; }
    [[nodiscard]] Properties &getProperties() { return properties_; }
    [[nodiscard]] const Properties &getProperties() const { return properties_; }

    void setInput(SpikeRate fusimotorFrequency) {
        i_.fusimotorFrequency = fusimotorFrequency.get();
    }
    void setFusimotorFrequency(DoubleT fusimotorFrequency) {
        i_.fusimotorFrequency = fusimotorFrequency;
    }
    void integrate(DoubleT dt);
    void validateState() { s_ = sNew_; }
    void calculateOutput() { /*intentionally left blank*/
    }
    void evaluate(DoubleT dt) {
        integrate(dt);
        calculateOutput();
        validateState();
    }

  private:
    Input i_;
    Parameters p_;
    State s_, sNew_;
    Properties properties_;
};

class Mileusnic2006MuscleSpindle {
  public:
    using type = nullptr_t;
    using concept_t = component_t;
    static constexpr std::string_view class_name = "Mileusnic2006MuscleSpindle";

    struct State {};

    struct Parameters {
        DoubleT primaryAfferentPartialOcclusion{ 0.156 };// `S`
    };

    struct Output {
        DoubleT primaryAfferent{ 0. };
        DoubleT secondaryAfferent{ 0. };
    };

    struct Input {
        DoubleT normalizedMuscleFiberLength = 1.;// `L` length of the muscle fiber
        DoubleT gStatic = 0.;// Static fusimotor frequency in pulses per second
        DoubleT gDynamic = 0.;// Dynamic fusimotor frequency in pulses per second
    };
    struct Properties {};

    Mileusnic2006MuscleSpindle();

    // void setActivation(DoubleT fStatic, DoubleT fDynamic);
    void setNormalizedMuscleFiberLength(DoubleT normalizedMuscleFiberLength);
    void setFusimotorFrequency(DoubleT gStatic, DoubleT gDynamic);
    void setInput(SpikeRate gStatic, SpikeRate gDynamic);
    void setInput(NormalizedFiberLength normalizedMuscleFiberLength);
    void setInput(Input input);
    // From input and current state calculate the new state of the system
    void integrate(DoubleT dt);
    // The temporary state calculated via `integrate` becomes the new state
    void validateState();
    // from the internal state of the system and the input, calculate all the
    // output;
    void calculateOutput();
    void evaluate(DoubleT dt);
    [[nodiscard]] std::string getName() const { return name_; }
    void setName(std::string_view name) { name_ = name; }
    void equilibrate();

    const auto &getBag1() const { return bag1_; }
    auto &getBag1() { return bag1_; }
    const auto &getBag2() const { return bag2_; }
    auto &getBag2() { return bag2_; }
    const auto &getChain() const { return chain_; }
    auto &getChain() { return chain_; }

    const auto &getDynamicsBag1() const { return dynamicsBag1_; }
    auto &getDynamicsBag1() { return dynamicsBag1_; }
    const auto &getDynamicsBag2() const { return dynamicsBag2_; }
    auto &getDynamicsBag2() { return dynamicsBag2_; }
    const auto &getDynamicsChain() const { return dynamicsChain_; }
    auto &getDynamicsChain() { return dynamicsChain_; }

    [[nodiscard]] const Output &getOutput() const { return o_; }
  private:
    Mileusnic2006IntrafusalFiber bag1_, bag2_, chain_;
    Mileusnic2006IntrafusalFiberActivationDynamics dynamicsBag1_, dynamicsBag2_, dynamicsChain_;
    Parameters p_;
    Input i_;
    Output o_;
    std::string name_;
};

}// namespace ceinms
#endif
