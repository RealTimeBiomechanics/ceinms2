#include "ceinms2/Dof.h"

namespace ceinms {

Dof::Dof(Parameters parameters)
    : p_(parameters) {
    i_.forces.resize(p_.n);
    i_.momentArms.resize(p_.n);
}

void Dof::setInput(const std::vector<Force> &values) {
    for (unsigned k{ 0 }; k < p_.n; ++k)
        i_.forces.at(k) = values.at(k).get();
}

void Dof::setInput(const Force &value, size_t slot) {
    i_.forces.at(slot) = value.get();
}

void Dof::setInput(const std::vector<MomentArm> &values) {
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

void connectSocket(const MomentArmsOnDof &parent, Dof &child) {
    child.setInput(parent);
}
}// namespace ceinms