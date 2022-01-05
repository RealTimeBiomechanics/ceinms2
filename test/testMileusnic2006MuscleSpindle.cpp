#include "ceinms2/Mileusnic2006MuscleSpindle.h"
#include "ceinms2/testingUtilities.h"
#include <iostream>
#include <fstream>

int test1() {
    ceinms::Mileusnic2006MuscleSpindle spindle;
    
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
    failed |= ceinms::runTest(&test1, "Correctness test");
    failed |= ceinms::runTest(&test2, "Mileusnic2006IntrafusalFiberActivationDynamics test");
    return failed;
}