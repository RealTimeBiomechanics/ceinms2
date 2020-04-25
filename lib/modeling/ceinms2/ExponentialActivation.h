#ifndef ceinms2_ExponentialActivation_h
#define ceinms2_ExponentialActivation_h
#include <cmath>
#include <algorithm>
#include <type_traits>
#include <string_view>
#include "ceinms2/Types.h"
namespace ceinms {



class ExponentialActivation {
  public:
    // this is required to allow some of the template magic in the `NMSmodel` class to work
    // it might be removed once we move to Concepts
    using type = nullptr_t;
    using concept_t = component_t;
    static constexpr std::string_view class_name = "ExponentialActivation";
    struct Parameters {
        Parameters() = default;
        DoubleT c1{ -0.5 };
        DoubleT c2{ -0.5 };
        DoubleT shapefactor{ -1 };
        DoubleT scalefactor{ 1. };
    };

    struct Input {
        Input() = default;
        DoubleT excitation{ 0. };
    };

    struct State {
        State() = default;
        DoubleT neuralActivationT1{ 0. }, neuralActivationT2{ 0. };
    };

    struct Output {
        Output() = default;
        DoubleT activation{ 0. };
        [[nodiscard]] DoubleT getPrimary() const { return activation; }
    };

    ExponentialActivation() { updateCoefficients(); };
    ExponentialActivation(Parameters parameters)
        : p_(parameters) {
        updateCoefficients();
    }

    void setInput(DoubleT excitation);
    void setInput(Input input);

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
    void evaluate(DoubleT dt);
    Parameters &updParameters() { return p_; }
    [[nodiscard]] Parameters getParameters() const { return p_; }
    State &updState() { return s_; }
    [[nodiscard]] State getState() const { return s_; }
    [[nodiscard]] Output getOutput() const { return o_; }

  private:
    /*Calculates internal coefficients based on the current `Parameters`*/
    void updateCoefficients();

    std::string name_;
    Parameters p_;
    State s_, sNew_;
    Input i_;
    Output o_;
    DoubleT alpha_, beta1_, beta2_, expShapefactor_;
};


void ExponentialActivation::setInput(DoubleT excitation) {
    excitation = std::max(0., excitation);
    excitation = std::min(1., excitation);
    i_.excitation = excitation;
}

void ExponentialActivation::setInput(Input input) { setInput(input.excitation); }

void ExponentialActivation::setState(State state) { s_ = state; }

void ExponentialActivation::integrate(DoubleT) {
    sNew_.neuralActivationT1 = (alpha_ * i_.excitation) - (beta1_ * s_.neuralActivationT1)
                               - (beta2_ * s_.neuralActivationT2);
}

void ExponentialActivation::validateState() {
    sNew_.neuralActivationT2 = s_.neuralActivationT1;
    s_ = sNew_;
}

void ExponentialActivation::calculateOutput() {
    o_.activation = p_.scalefactor * (std::exp(p_.shapefactor * s_.neuralActivationT1) - 1)
                    / (expShapefactor_ - 1);
}

void ExponentialActivation::evaluate(DoubleT dt) {
    integrate(dt);
    validateState();
    calculateOutput();
}

void ExponentialActivation::updateCoefficients() {
    beta1_ = p_.c1 + p_.c2;
    beta2_ = p_.c1 * p_.c2;
    alpha_ = 1 + beta1_ + beta2_;
    expShapefactor_ = std::exp(p_.shapefactor);
}


void connectSocket(const Excitation &parent, ExponentialActivation &child) {
    child.setInput(parent.value);
}


}// namespace ceinms

#endif