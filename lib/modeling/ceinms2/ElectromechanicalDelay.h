#ifndef ceinms2_ElectromechanicalDelay_h
#define ceinms2_ElectromechanicalDelay_h
#include <ceinms2/Types.h>
#include <ceinms2/Curve.h>

namespace ceinms {

class ElectromechanicalDelay {
  public:
    using type = nullptr_t;
    using concept_t = component_t;
    static constexpr std::string_view class_name = "ElectromechanicalDelay";
    struct Parameters {
        DoubleT delay{ 0. };
    };

    struct Output {
        DoubleT excitation;
    };

    struct State {
        Curve<CurveMode::Online, CurveMode::Linear, 2000> excitationNodes;
        DoubleT elapsedTime{ 0. };
    };

    struct Input {
        DoubleT excitation;
    };

    ElectromechanicalDelay() = default;
    ElectromechanicalDelay(const Parameters &p)
        : p_(p) {
        validateState();
    }


    void setInput(Excitation value) { i_.excitation = value;
    }
    void setExcitation(DoubleT value) { i_.excitation = value; }
    void setState(const State& state) { s_ = state;
    }
    // From input and current state calculate the new state of the system
    void integrate(DoubleT dt) { 
        tempElapsedTime_ = s_.elapsedTime;
        tempElapsedTime_ += dt;
        s_.excitationNodes.removeLastPointNoUpdate();
        s_.excitationNodes.addPoint(tempElapsedTime_, i_.excitation);
    }
    // from the internal state of the system and the input, calculate all the
    // output;
    void calculateOutput() { 
        o_.excitation = s_.excitationNodes.getValue(tempElapsedTime_ - p_.delay);
    }

    // The temporary state calculated via `integrate` becomes the new state
    void validateState() { 
        s_.elapsedTime = tempElapsedTime_;
        //Add buffer point at the end. This is going to be removed the next time `integrate` is called
        s_.excitationNodes.addPointNoUpdate(std::numeric_limits<DoubleT>::max(), 0);
    }
    void evaluate(DoubleT dt) {
        integrate(dt);
        calculateOutput();
        validateState();
    }

    [[nodiscard]] std::string getName() const { return name_; }
    void setName(std::string name) { name_ = name; }

    Parameters &updParameters() { return p_; }
    [[nodiscard]] Parameters getParameters() const { return p_; }
    State &updState() { return s_; }
    [[nodiscard]] State getState() const { return s_; }
    [[nodiscard]] Output getOutput() const { return o_; }

    template<typename T, std::enable_if_t<std::is_same<T, Excitation>::value, int> = 0>
    [[nodiscard]] Excitation getOutput() const {
        return o_.excitation;
    }


  private:
    Parameters p_;
    Output o_;
    Input i_;
    State s_;
    DoubleT tempElapsedTime_{ 0. };
    std::string name_;
};

void connectSocket(const Excitation &parent, ElectromechanicalDelay &child);
}// namespace ceinms

#endif