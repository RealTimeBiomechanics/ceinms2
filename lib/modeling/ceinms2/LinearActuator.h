#ifndef ceinms2_LinearActuator_h
#define ceinms2_LinearActuator_h
#include <ceinms2/Types.h>
#include <string_view>
#include <string>

namespace ceinms {

class LinearActuator {
  public:
    using type = nullptr_t;
    using concept_t = component_t;
    static constexpr std::string_view class_name = "LinearActuator";

    struct Parameters {
        DoubleT maximumForce{ 1. };
    };

    struct Input {
        DoubleT activation{ 0. };
    };

    struct Output {
        DoubleT force{ 0. };
    };

    void setActivation(DoubleT activation) { i_.activation = activation;
    }
    void setInput(Activation value) { i_.activation = value.get();
    }

    // From input and current state calculate the new state of the system
    void integrate(DoubleT){};

    // from the internal state of the system and the input, calculate all the
    // output;
    void calculateOutput() { o_.force = i_.activation * p_.maximumForce;
    }
    [[nodiscard]] std::string getName() const { return name_; }
    void setName(std::string_view name) { name_ = name; }
    // Convenience function that, from input and current state, calculate the
    // new state and all the output
    void evaluate(DoubleT) { calculateOutput();
    }

    [[nodiscard]] Parameters &getParameters() { return p_; }
    [[nodiscard]] const Parameters &getParameters() const { return p_; }
    [[nodiscard]] const Output &getOutput() const { return o_; }
    template<typename T, std::enable_if_t<std::is_same<T, Force>::value, int> = 0>
    [[nodiscard]] Force getOutput() const {
        return Force{ o_.force };
    }

  private:
    std::string name_;
    Parameters p_;
    Input i_;
    Output o_;
};

template<ActivationGenerator T>
void connectSocket(const T &parent, LinearActuator &child) {
    child.setActivation(parent.getOutput().activation);
}

void connectSocket(const Activation &parent, LinearActuator &child) {
    child.setInput(parent);
}

}


#endif
