#include "ceinms2/Mileusnic2006MuscleSpindle.h"
#include "ceinms2/Lloyd2003Muscle.h"
#include "ceinms2/testingUtilities.h"
#include <iostream>
#include <fstream>



int test0() {
    std::ofstream outF("test1.csv");
    ceinms::Lloyd2003Muscle::Parameters mtuParameters;
    mtuParameters.optimalFiberLength = 0.038;
    mtuParameters.tendonSlackLength = 0.071;
    mtuParameters.maxIsometricForce = 1;
    double maxMtuLength = 0.11275;
    double m = 0.005;
    double dt = 0.001;
    ceinms::Lloyd2003Muscle mtu{ mtuParameters };
    ceinms::Mileusnic2006MuscleSpindle spindle;
//    spindle.getBag1().getParameters().beta0 = 0.6;

    for (int i = 0; i < 4000; ++i) {
        double t = i * dt;
        double l = 0.96 * maxMtuLength;
        if (t > 1.) l = 0.96 * maxMtuLength + m * (t - 1.);
        if (l > maxMtuLength) l = maxMtuLength;
        mtu.setActivation(0);
        mtu.setMusculotendonLength(l);
        if (i == 0) { 
            mtu.equilibrate();
            spindle.equilibrate();
        }
        mtu.evaluate(dt);
        spindle.setNormalizedMuscleFiberLength(mtu.getOutput().normalizedFiberLength);
        spindle.setFusimotorFrequency(0, 0);
        spindle.evaluate(dt);
        outF << t << "," << l << "," 
             << spindle.getOutput().primaryAfferent<< ","
             << spindle.getBag1().getState().polarRegionVelocity << ","
            << spindle.getBag1().getOutput().fspr << ","
             << spindle.getBag1().getOutput().fvpr << ","
             << spindle.getBag1().getOutput().fiberTension              << std::endl;
    }

    return 1;
}

int test1() {
    std::ofstream outF("test1.csv");
    double dt = 0.001;
    double m = 0.11;
    double maxLength = 1.08;
    ceinms::Mileusnic2006MuscleSpindle spindle;
        
    for (int i = 0; i < 4000; ++i) {
        double t = i * dt;
        double l = 0.95;
        if (t > 1.) l = 0.95 + m * (t - 1.);
        if (l > maxLength) l = maxLength;
        spindle.setNormalizedMuscleFiberLength(l);
        spindle.setFusimotorFrequency(0, 0);
        spindle.evaluate(dt);
        outF << t << "," << l << "," << spindle.getBag1().getOutput().primaryAfferent << ","
             << spindle.getBag1().getState().polarRegionVelocity << ","

             << spindle.getOutput().primaryAfferent << "," << spindle.getOutput().secondaryAfferent
             << std::endl;
        }
        
    return 1;
}


int test2() {
    ceinms::Mileusnic2006IntrafusalFiberActivationDynamics dynamics;
    dynamics.getParameters().freq = 60;
    dynamics.getParameters().tau = 0.149;
    dynamics.setFusimotorFrequency(200);
    std::ofstream outF("mileusnic2006IntrafusalFiberActivationDynamics.csv");
    for (int i = 0; i < 1000; ++i) {
        dynamics.evaluate(0.01);
        outF << dynamics.getOutput().activation << std::endl;
    }
    return 1;
}

int main() {
    bool failed = false;
    failed |= ceinms::runTest(&test0, "Correctness test");
   // failed |= ceinms::runTest(&test0, "Mileusnic2006IntrafusalFiberActivationDynamics test");
    return failed;
}