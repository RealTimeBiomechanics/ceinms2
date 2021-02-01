#include <ceinms2/NMSmodel.h>
#include <ceinms2/Dof.h>
#include <ceinms2/Lloyd2003Muscle.h>
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

int testDofConnection() {
    using MyNMSmodel = NMSmodel<Lloyd2003Muscle, Dof>;
    MyNMSmodel model;
    ceinms::Lloyd2003Muscle muscle(getDefaultMuscle());
    muscle.setName("mtu1");
    model.addComponent(muscle);
    muscle.setName("mtu2");
    model.addComponent(muscle);
    muscle.setName("mtu3");
    model.addComponent(muscle);

    model.addInput<Activation>("mtu1");
    model.addInput<Activation>("mtu2");
    model.addInput<Activation>("mtu3");
    model.addInput<MusculotendonLength>("mtu1");
    model.addInput<MusculotendonLength>("mtu2");
    model.addInput<MusculotendonLength>("mtu3");

    ceinms::Dof dof({ 3 });
    dof.setName("dof1");
    model.addComponent(dof);

    model.addInput<MomentArmsOnDof>("dof1");

    model.connect<Activation, Lloyd2003Muscle>();
    model.connect<MusculotendonLength, Lloyd2003Muscle>();
    model.connect<Lloyd2003Muscle, Dof>(Socket("mtu1"), Socket{ "dof1", 0 });
    model.connect<Lloyd2003Muscle, Dof>(Socket("mtu2"), Socket{ "dof1", 1 });
    model.connect<Lloyd2003Muscle, Dof>(Socket("mtu3"), Socket{ "dof1", 2 });

    model.connect<MomentArmsOnDof, Dof>();


    model.setInput(vector{ Activation{ 1. }, Activation{ 1. }, Activation{ 1. } });
    model.setInput(
        vector{ MusculotendonLength{ 2. }, MusculotendonLength{ 2. }, MusculotendonLength{ 2. } });
    model.setInput(vector{ MomentArmsOnDof{ { 1., 1., 1. } } });
    model.evaluate(0.01);
    std::cout << model.getComponent<Dof>("dof1").getOutput<Torque>().get() << std::endl;
    return 0;
}

int testDof() {
    ceinms::Dof d({ 3 });
    std::vector<DoubleT> forces{ 1., 2., 3. };
    std::vector<DoubleT> momentArms{ 1., 1. / 2., 1. / 3. };

    d.setForces(forces);
    d.setMomentArms(momentArms);
    d.calculateOutput();
    std::cout << d.getOutput<Torque>().get() << std::endl;

    d.setInput(std::vector{ Force{ 1. }, Force{ 2. }, Force{ 3. } });
    d.setInput(std::vector{ MomentArm{ 2. }, MomentArm{ 2. / 2. }, MomentArm{ 2. / 3. } });
    d.calculateOutput();
    std::cout << d.getOutput<Torque>().get() << std::endl;

    return 0;
}

int main() {
    testDof();
    testDofConnection();

    return 0;
}