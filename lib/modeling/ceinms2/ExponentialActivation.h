#ifndef ceinms2_ExponentialActivation_h
#define ceinms2_ExponentialActivation_h
#include <cmath>
#include <algorithm>
#include <type_traits>
#include <string>
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
    };

    ExponentialActivation() { updateCoefficients(); };
    ExponentialActivation(Parameters parameters)
        : p_(parameters) {
        updateCoefficients();
    }

    void setInput(Excitation value);
    void setExcitation(DoubleT value);

    void setState(State state);
    // From input and current state calculate the new state of the system
    void integrate(DoubleT dt);
    // from the temporary internal state of the system and the input, calculate all the
    // output; This is done on the temporary state so we can probe different input wihtout affecting the state
    void calculateOutput();
    // The temporary state calculated via `integrate` becomes the new state
    void validateState();

    [[nodiscard]] std::string getName() const { return name_; }
    void setName(std::string name) { name_ = name; }
    void evaluate(DoubleT dt);
  //  Parameters &getParameters() { return p_; }
    [[nodiscard]] Parameters getParameters() const { return p_; }
    void setParameters(const Parameters &p) {
        p_ = p;
        updateCoefficients();
    }
    
    State &getState() { return s_; }
    [[nodiscard]] State getState() const { return s_; }
    [[nodiscard]] Output getOutput() const { return o_; }
    [[nodiscard]] Input &getInput() { return i_; }
    [[nodiscard]] const Input &getInput() const { return i_; }
    template<typename T, std::enable_if_t<std::is_same<T, Activation>::value, int> = 0>
    [[nodiscard]] Activation getOutput() const {
        return Activation{ o_.activation };
    }

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

template<ExcitationGenerator T>
void connectSocket(const T &parent, ExponentialActivation &child) {
    child.setExcitation(parent.getOutput().excitation);
}

void connectSocket(const Excitation &parent, ExponentialActivation &child);

}// namespace ceinms

#endif
