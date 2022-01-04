#ifndef ceinms2_RienerActivation_h
#define ceinms2_RienerActivation_h
#include <ceinms2/Types.h>
#include <boost/numeric/odeint.hpp>
#include <array>
#include <cmath>
#include <numbers>

namespace ceinms {

/*Implementation of spike train to motor unit activation model.
* Riener, R. and J. Quintern (1997). "A physiologically based model of
muscle activation verified by electrical stimulation." Bioelectrochemistry
and Bioenergetics 43(2): 257-264.
*/

class RienerActivation {

  public:

    struct Parameters {
        // Coefficient for T-tubule membrane model
        DoubleT c1{ 1.84e4 }, c2{ 1.2e7 }, c3{ 0.93 };

        // Coefficient for calcium dynamics model
        DoubleT c4{ 2.24e4 }, c5{ 3.8e5 }, c6{ 4.3e9 };

        // Fatigue coefficients
   //     DoubleT tFat{ 7. }, tRec{ 11. }, fitMin{ 0. };

        //Calcium binding constant
   //     DoubleT rho{ 6.62e4 };
    };

    struct Input {
        bool pulse{ 0 };
    };

    struct Output {
        DoubleT activation{ 0. };
    };

    struct State {
        using type = std::array<DoubleT, 2>;
        DoubleT alpha{ 0. };
        //Depolarisation of T-tubule membrane
        // beta[0] is beta
        // beta[1] is beta_dot
        type beta{ 0., 0. };

        // Release of calcium ion from sarcoplasmic reticumum
        // gamma[0] is gamma
        // gamma[1] is gamma_dot
        type gamma{ 0., 0. };
    };
    static void calculateDepolarization(
        const typename State::type &x,
        typename State::type &dxdt,
        DoubleT pulse,
        const Parameters& p) {

        dxdt[0] = x[1];
        dxdt[1] = p.c3 * pulse - p.c1 * x[1] - p.c2 * x[0];
    }

    static void calculateCalciumRelease(const typename State::type &x,
        typename State::type &dxdt,
        DoubleT depolarization,
        const Parameters &p) {
        dxdt[0] = x[1];
        dxdt[1] = p.c6 * depolarization - p.c4 * x[1] - p.c5 * x[0];
    }

    static DoubleT calculateNervePotential(DoubleT time) {
        const DoubleT halfPeriod = 0.001;
        if (time < halfPeriod) return std::sin(1000. * std::numbers::pi * time);
        return 0;
    }


    void integrate(DoubleT dt) {

        rk1.do_step(
            [&](const typename State::type &x, typename State::type &dxdt, DoubleT) {
                RienerActivation::calculateDepolarization(x, dxdt, i_.pulse, p_);
            },
            s_.beta,
            0.,
            sNew_.beta,
            dt);
        rk2.do_step(
            [&](const typename State::type &x, typename State::type &dxdt, DoubleT) {
                RienerActivation::calculateCalciumRelease(x, dxdt, sNew_.beta[0], p_);
            },
            s_.gamma,
            0.,
            sNew_.gamma,
            dt);
        s_ = sNew_;
        
    }
    boost::numeric::odeint::runge_kutta4<State::type> rk1, rk2;
    Parameters p_;
    State s_, sNew_;
    Input i_;

};

}// namespace ceinms

#endif