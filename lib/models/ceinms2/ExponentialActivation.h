#ifndef ceinms2_ExponentialActivation_h
#define ceinms2_ExponentialActivation_h

namespace ceinms {

class Lloyd2019Activation {
  public:
    using DoubleT = double;
    struct Parameters {
        Parameters()
            : c1(-0.5)
            , c2(-0.5)
            , shapefactor(-1) {}
        DoubleT c1;
        DoubleT c2;
        DoubleT shapefactor;
    };

    struct Input {
        Input()
            : excitation(0) {}

        DoubleT excitation;
    };

    struct State {
        State()
            : fiberVelocity(0) {}

        DoubleT fiberLength;
        DoubleT fiberVelocity;
    };

    struct Output {
        DoubleT activation;
    };
};
}// namespace ceinms

#endif