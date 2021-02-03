#include "ceinms2/LinearActuator.h"
#include "ceinms2/NMSmodel.h"
#include "ceinms2/testingUtilities.h"
#include <iostream>

int test1() {
    ceinms::LinearActuator actuator;
    actuator.getParameters().maximumForce = 100;
    actuator.setActivation(0.5);
    actuator.evaluate(0.1);
    return !(actuator.getOutput<ceinms::Force>() == 50);
}

int test2() {
    ceinms::NMSmodel<ceinms::LinearActuator> model;
    ceinms::LinearActuator actuator;
    actuator.getParameters().maximumForce = 100;
    actuator.setName("mtu1");
    model.addComponent<ceinms::LinearActuator>(actuator);
    model.addInput<ceinms::Activation>("mtu1");
    model.connect<ceinms::Activation, ceinms::LinearActuator>();
    model.setInput(std::vector<ceinms::Activation>{ 1});
    model.evaluate(0.1);
    return !(model.getOutput<ceinms::LinearActuator>().at(0).force == 100.);
}

int main() {
    bool failed = false;
    failed |= ceinms::runTest(&test1, "Correctness test");
    failed |= ceinms::runTest(&test2, "Correctness test within NMS model");
    return failed;
}