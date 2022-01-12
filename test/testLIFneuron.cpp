#include "ceinms2/LIFneuron.h"
#include "ceinms2/testingUtilities.h"
#include <iostream>
#include <fstream>
#include <random>
#include <algorithm>

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

namespace ceinms {
template<typename T>
class Population {
  public:
    using iterator = typename std::vector<T>::iterator;
    using const_iterator = typename std::vector<T>::const_iterator;
    Population(unsigned int size)
        : population_(size)
        , gen_(0) {}

    Population(unsigned int size, const T &defaultItem)
        : population_(size, defaultItem)
        , gen_(0) {}

    template<typename D, typename L>
    void distributeParameters(D &distribution, const L &function, const DoubleT offset = 0.) {
        for (auto &p : population_) {
            function(p) = distribution(gen_) + offset;
        }
        std::sort(population_.begin(), population_.end(), [&function](auto &lhs, auto &rhs) {
            return function(lhs) < function(rhs);
        });
    }


    template<typename T>
    void setInput(T value) {
        for (auto &p : population_) {
            p.setInput(value);
        }
    }

    void integrate(DoubleT dt) {
        for (auto &p : population_) {
            p.integrate(dt);
        }
    }

    void calculateOutput() {
        for (auto &p : population_) {
            p.calculateOutput();
        }
    }

    void validateState() {
        for (auto &p : population_) {
            p.validateState();
        }
    }

    void evaluate(DoubleT dt) {
        for (auto &p : population_) {
            p.evaluate(dt);
        }
    }

    iterator begin() { return population_.begin(); }
    const_iterator begin() const { return population_.begin(); }
    iterator end() { return population_.end(); }
    const_iterator end() const { return population_.end(); }

  private:
    std::vector<T> population_;
    std::mt19937 gen_;
};
}// namespace ceinms
int test2() {
    std::ofstream outF("LIFneuronPool.csv");
    ceinms::Population<ceinms::LIFneuron> neuronPool(50);
    auto resistance(
        [](ceinms::LIFneuron &neuron) -> ceinms::DoubleT& { return neuron.getParameters().membraneResistance; });
    //std::uniform_real_distribution dis(1e6, 20e6);
    std::normal_distribution dis(10e6, 3e6);
    neuronPool.distributeParameters(dis, resistance);
   
    double dt = 0.0001;
    long double imin = 1e-9;
    long double imax = 5e-9;
    int N = 10000;
    for (int i{ 0 }; i < N; ++i) {
        const double curr = imin + (imax - imin) / N * i;
        neuronPool.setInput(ceinms::Current{ curr });
        neuronPool.evaluate(dt);
        outF <<  curr << ",";
        for (const auto &neuron : neuronPool) {
            if (neuron.getOutput().spike >= 1.)
                outF << i * dt << ',';
            else
                outF << "NaN,";
        }
        outF << '\n';
    }
    return 0;
}



int main() {
    bool failed = false;
    failed |= ceinms::runTest(&test1, "Correctness test");
    failed |= ceinms::runTest(&test2, "Pool test");

    return failed;
}