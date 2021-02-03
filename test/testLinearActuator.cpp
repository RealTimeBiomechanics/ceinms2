#include "ceinms2/LinearActuator.h"
#include "ceinms2/NMSmodel.h"
#include <iostream>

int main() {
    ceinms::LinearActuator actuator;
    actuator.getParameters().maximumForce = 100;
    actuator.setActivation(0.5);
    actuator.evaluate(0.1);
    
    std::cout << actuator.getOutput<ceinms::Force>() << " should be 50" << std::endl;

    ceinms::NMSmodel<ceinms::LinearActuator> model;
    actuator.setName("mtu1");
    model.addComponent<ceinms::LinearActuator>(actuator);
    model.addInput<ceinms::Activation>("mtu1");
    model.connect<ceinms::Activation, ceinms::LinearActuator>();
    model.setInput(std::vector<ceinms::Activation>{ 1});
    model.evaluate(0.1);
    std::cout << model.getOutput<ceinms::LinearActuator>().at(0).force << std::endl;

    return 0;
}