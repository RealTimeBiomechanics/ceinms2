#include "ceinms2/LIFneuron.h"
#include "ceinms2/testingUtilities.h"
#include <iostream>
#include <fstream>

int test1() {
    ceinms::LIFneuron neuron;
    std::ofstream outF("LIFneuron.csv");
    double dt = 0.0001;
    long double imin = 3e-9;
    long double imax = 5e-9;
    int N = 1000;
    for (int i{ 0 }; i < N; ++i) {
        const double curr = imin + (imax - imin) / N * i;
        neuron.setCurrent(curr);
        neuron.evaluate(dt);
        outF << curr<< "," <<neuron.getOutput().spike << ',' << neuron.getState().membranePotential << "\n";
    }
    return 0;
}


int main() {
    bool failed = false;
    failed |= ceinms::runTest(&test1, "Correctness test");
    return failed;
}