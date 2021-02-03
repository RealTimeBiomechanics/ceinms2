#include <ceinms2/ElectromechanicalDelay.h>
#include <ceinms2/testingUtilities.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <numeric>
#include <fstream>
using namespace ceinms;

int test1() {
    DoubleT delay{ 0.05 };
    ElectromechanicalDelay emDelay({ delay });
    DoubleT dt{ 0.001 };
    std::vector<DoubleT> out, in, truth;
    for (int i{ 0 }; i < 2000; ++i) {
        in.push_back(std::sin(i * dt));
        if (i * dt < delay)
            truth.push_back(0.);
        else
            truth.push_back(std::sin(i * dt - delay));
        emDelay.setExcitation(in.back());
        emDelay.evaluate(dt);
        out.push_back(emDelay.getOutput<Excitation>());
    }

    double rmse = std::sqrt(std::transform_reduce(out.begin(),
                                out.end(),
                                truth.begin(),
                                0.,
                                std::plus<>(),
                                [](auto a, auto b) { return (a - b) * (a - b); })
                            / out.size());

    return !(rmse < 1e-8);
}


int main() {
    bool failed = false;
    failed |= runTest(&test1, "Test correctness");
    return failed;
}