#ifndef ceinms2_Dof_h
#define ceinms2_Dof_h

#include <vector>
#include <string>
#include <numeric>
#include "ceinms2/Types.h"

namespace ceinms {
class Dof {
  public:
    // this is required to allow some of the template magic in the `NMSmodel` class to work
    // it might be removed once we move to Concepts
    using type = nullptr_t;
    using concept_t = component_t;
    static constexpr std::string_view class_name = "Dof";

    struct Input {
        Input() = default;
        std::vector<DoubleT> forces;
        std::vector<DoubleT> momentArms;
    };

    struct Parameters {
        unsigned n{ 1 };
    };

    struct Output {
        Output() = default;
        DoubleT torque{ 0. };
    };

    Dof(Parameters parameters);

    [[nodiscard]] std::string getName() const { return name_; }
    void setName(std::string name) { name_ = name; }

    void setInput(const std::vector<Force> &values);
    void setInput(const Force &value, size_t slot);
    void setInput(const std::vector<MomentArm> &values);
    void setInput(const MomentArmsOnDof &values);

    void setForces(const std::vector<DoubleT> &values);
    void setMomentArms(const std::vector<DoubleT> &values);

    void integrate(DoubleT) {
    //intentionally blank
    }
    void validateState() {
    // intentionally blank
    }

    void evaluate(DoubleT) { calculateOutput(); }
    void calculateOutput();

    [[nodiscard]] const Output &getOutput() const { return o_; }

    template<typename T, std::enable_if_t<std::is_same<T, Torque>::value, int> = 0>
    [[nodiscard]] Torque getOutput() const { return Torque{ o_.torque }; }

  private:
    std::string name_;
    Output o_;
    Input i_;
    Parameters p_;
};


Dof::Dof(Parameters parameters)
    : p_(parameters) {
    i_.forces.resize(p_.n);
    i_.momentArms.resize(p_.n);
}


void Dof::setInput(const std::vector<Force>& values) {
    for (unsigned k{ 0 }; k < p_.n; ++k)
        i_.forces.at(k) = values.at(k).get();
 }

void Dof::setInput(const Force &value, size_t slot) {
     i_.forces.at(slot) = value.get();
 }


void Dof::setInput(const std::vector<MomentArm>& values) {
     for (unsigned k{ 0 }; k < p_.n; ++k)
         i_.momentArms.at(k) = values.at(k).get();
 }

void Dof::setInput(const MomentArmsOnDof &values) {
     i_.momentArms = values.get();
 }



void Dof::setForces(const std::vector<DoubleT> &values) {
    i_.forces = values;
 }

void Dof::setMomentArms(const std::vector<DoubleT> &values) {
     i_.momentArms = values;
 }


void Dof::calculateOutput() {
     o_.torque =
         std::inner_product(i_.forces.cbegin(), i_.forces.cend(), i_.momentArms.cbegin(), 0.);
 }


template<ForceGenerator F>
 void connectSocket(const F &parent, Dof &child, size_t childSlot) {
    child.setInput(parent.getOutput<Force>(), childSlot);
}

 void connectSocket(const MomentArmsOnDof &parent, Dof &child) {
    child.setInput(parent);
}

}
#endif