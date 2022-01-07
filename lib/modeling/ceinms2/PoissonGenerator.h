#ifndef ceinms_PoissonGenerator_h
#define ceinms_PoissonGenerator_h
#include <random>
#include "ceinms2/Types.h"

namespace ceinms {
class PoissonGenerator {
  public:
    PoissonGenerator()
        : dis_(.0, 1.0)
    {
        std::random_device rd;
        gen_.seed(rd());
        validateState();
    }

    struct Input {
        DoubleT spikeRate = 0.;
    };

    struct Parameters {};

    struct State {
        DoubleT occurrence = 1.;
        DoubleT probability = 0.;
    };

    struct Output {
        DoubleT spike = 0.;
    };

    void setInput(SpikeRate rate) { i_.spikeRate = rate; }
    void setSpikeRate(DoubleT rate) { i_.spikeRate = rate; }
    [[nodiscard]] const Output &getOutput() const { return o_; }
    void integrate(DoubleT dt) { 

        sNew_.probability = i_.spikeRate * dt;
    }

    // Update the occurrence only when validating the state.
    // If the generator is used in some large optimisation, 
    // we need it to produce the same results at every call of integrate
    void validateState() {
        sNew_.occurrence = dis_(gen_);
    }

    void calculateOutput() { 
        if (sNew_.occurrence < sNew_.probability)
            o_.spike = 1.;
        else
            o_.spike = 0.;
    }

    void evaluate(DoubleT dt) { 
        integrate(dt);
        calculateOutput();
        validateState();
    }

  private:
    Input i_;
    State sNew_;
    Output o_;
    std::mt19937 gen_;
    std::uniform_real_distribution<> dis_;
};
}


#endif