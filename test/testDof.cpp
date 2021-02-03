#include <ceinms2/NMSmodel.h>
#include <ceinms2/Dof.h>
#include <ceinms2/Lloyd2003Muscle.h>
#include <ceinms2/LinearActuator.h>
#include <ceinms2/TestingUtilities.h>
#include <iostream>

using namespace ceinms;

auto getDefaultMuscle() {
    ceinms::Lloyd2003Muscle::Parameters p;
    p.damping = 0.1;
    p.maxContractionVelocity = 1;
    p.maxIsometricForce = 1;
    p.optimalFiberLength = 1;
    p.pennationAngleAtOptimalFiberLength = 0;
    p.percentageChange = 0.15;
    p.strengthCoefficient = 1;
    p.tendonSlackLength = 1;
    ceinms::Lloyd2003Muscle muscle(p);
    muscle.setName("muscle");
    return muscle;
}

int testDofConnectionWithNMSmodel() {
    using MyNMSmodel = NMSmodel<LinearActuator, Dof>;
    MyNMSmodel model;
    ceinms::LinearActuator muscle;
    muscle.setName("mtu1");
    model.addComponent(muscle);
    muscle.setName("mtu2");
    model.addComponent(muscle);
    muscle.setName("mtu3");
    model.addComponent(muscle);

    model.addInput<Activation>("mtu1");
    model.addInput<Activation>("mtu2");
    model.addInput<Activation>("mtu3");

    ceinms::Dof dof({ 3 });
    dof.setName("dof1");
    model.addComponent(dof);

    model.addInput<MomentArmsOnDof>("dof1");

    model.connect<Activation, LinearActuator>();
    model.connect<LinearActuator, Dof>(Socket("mtu1"), Socket{ "dof1", 0 });
    model.connect<LinearActuator, Dof>(Socket("mtu2"), Socket{ "dof1", 1 });
    model.connect<LinearActuator, Dof>(Socket("mtu3"), Socket{ "dof1", 2 });

    model.connect<MomentArmsOnDof, Dof>();

    model.setInput(vector<Activation>{ 1.,  1., 1. });
    model.setInput(vector{ MomentArmsOnDof{ { 1., 1., 1. } } });
    model.evaluate(0.01);

    return !(model.getComponent<Dof>("dof1").getOutput<Torque>() == 3.);
 
}

int testDof() {
    ceinms::Dof d({ 3 });
    std::vector<DoubleT> forces{ 1., 2., 3. };
    std::vector<DoubleT> momentArms{ 1., 1. / 2., 1. / 3. };

    d.setForces(forces);
    d.setMomentArms(momentArms);
    d.calculateOutput();
    std::cout << d.getOutput<Torque>() << std::endl;

    d.setInput(std::vector{ Force{ 1. }, Force{ 2. }, Force{ 3. } });
    d.setInput(std::vector{ MomentArm{ 2. }, MomentArm{ 2. / 2. }, MomentArm{ 2. / 3. } });
    d.calculateOutput();
    std::cout << d.getOutput<Torque>() << std::endl;

    return 0;
}

int main() {
    bool failed = false;
    failed |= runTest(&testDof, "Compile test");
    failed |= runTest(&testDofConnectionWithNMSmodel, "Works within NMS model");
    return failed;
}