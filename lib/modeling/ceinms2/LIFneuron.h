#ifndef ceinms2_LIFneuron_h
#define ceinms2_LIFneuron_h
#include "ceinms2/Types.h"
#include <boost/numeric/odeint.hpp>
#include <random>

namespace ceinms {
/* This class represents a leaky integrate-and-fire neuron*/
class LIFneuron {
  public:
      LIFneuron() :
      dis_(0., p_.stochasticThresholdStd){
        std::random_device rd;
        gen_.seed(rd());
        validateState();
    }
    struct Parameters {
        DoubleT tau = 8e-3;// s
        DoubleT restingPotential = -70e-3;// V
        //DoubleT resetVoltage = -65e-3;// V
        DoubleT firingThreshold = -50e-3;// V
        //DoubleT absoluteRefactoryPeriod = 2e-3;// s
        DoubleT membraneResistance = 10e6; //ohm

        // standard deviation of the stochastic noise added to the firing threshold. This is used to
        // simulate the variability observed in real neurons
        DoubleT stochasticThresholdStd = 1e-3; };

    struct State {
        DoubleT membranePotential = -70e-3;
        DoubleT spike = 0.;
        DoubleT thresholdNoise = 0.;
    };

    struct Input {
        DoubleT current = 0.;
    };

    struct Output {
        DoubleT spike = 0.;
    };

    void setCurrent(DoubleT current) { i_.current = current; }
    void setInput(Current current) { i_.current = current.get(); }
    [[nodiscard]] const Output &getOutput() const { return o_; }
    [[nodiscard]] const State &getState() const { return s_; }
    [[nodiscard]] Parameters &getParameters() { return p_; }

    void integrate(DoubleT dt) {
        boost::numeric::odeint::runge_kutta4<DoubleT> stepper;
        stepper.do_step(
            [&](const DoubleT &x, DoubleT &dxdt, const DoubleT /* t */) { 
                dxdt = (-(x - p_.restingPotential) + p_.membraneResistance * i_.current) / p_.tau;
            },
            s_.membranePotential,
            0.,
            sNew_.membranePotential,
            dt);
        
        if (sNew_.membranePotential > p_.firingThreshold + s_.thresholdNoise) {
            sNew_.spike = 1.;
            sNew_.membranePotential = p_.restingPotential;
        }
        else
            sNew_.spike = 0;
            
    }

    void calculateOutput() {
            o_.spike = sNew_.spike;
       }

    void validateState() { 
        s_ = sNew_;
        sNew_.thresholdNoise = dis_(gen_);
    }

    void evaluate(DoubleT dt) { 
        integrate(dt);
        calculateOutput();
        validateState();
    }


  private:
    State s_, sNew_;
    Parameters p_;
    Input i_;
    Output o_;
    std::mt19937 gen_;
    std::normal_distribution<> dis_;
};
}// namespace ceinms

#endif