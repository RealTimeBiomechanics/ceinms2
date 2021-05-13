#include "ceinms2/ExponentialActivation.h"

namespace ceinms {

void ExponentialActivation::setExcitation(DoubleT value) {
    i_.excitation = std::clamp(value, 0., 1.);
}

void ExponentialActivation::setInput(Excitation value) {
    setExcitation(value.get());
}

void ExponentialActivation::setState(State state) {
    s_ = state;
}

void ExponentialActivation::integrate(DoubleT) {
    sNew_.neuralActivationT1 = (alpha_ * i_.excitation) - (beta1_ * s_.neuralActivationT1)
                               - (beta2_ * s_.neuralActivationT2);
}

void ExponentialActivation::validateState() {
    sNew_.neuralActivationT2 = s_.neuralActivationT1;
    s_ = sNew_;
}

void ExponentialActivation::calculateOutput() {
    o_.activation = p_.scalefactor * (std::exp(p_.shapefactor * sNew_.neuralActivationT1) - 1)
                    / (expShapefactor_ - 1);
}

void ExponentialActivation::evaluate(DoubleT dt) {
    integrate(dt);
    calculateOutput();
    validateState();
}

void ExponentialActivation::updateCoefficients() {
    beta1_ = p_.c1 + p_.c2;
    beta2_ = p_.c1 * p_.c2;
    alpha_ = 1 + beta1_ + beta2_;
    expShapefactor_ = std::exp(p_.shapefactor);
}

void connectSocket(const Excitation &parent, ExponentialActivation &child) {
    child.setInput(parent);
}


}