#include "ceinms2/PoissonGenerator.h"
#include "ceinms2/testingUtilities.h"
#include <iostream>
#include <fstream>

int test1() {
    ceinms::PoissonGenerator gen;
    std::ofstream outF("poissonGenerator.csv");
    double dt = 0.001;
    for (int i{ 0 }; i < 5000; ++i) {
        gen.setInput(25+25*std::sin(i*dt/5.));
        gen.evaluate(dt);
        outF << gen.getOutput().spike << '\n';
    }
    return 0;
}

int test2() {
    return 0;
}

int main() {
    bool failed = false;
    failed |= ceinms::runTest(&test1, "Correctness test");
    failed |= ceinms::runTest(&test2, "Correctness test within NMS model");
    return failed;
}