#include <ceinms2/Lloyd2019Muscle.h>
#include <ceinms2/ExponentialActivation.h>
#include <ceinms2/NMSmodel.h>

using namespace std;
using namespace ceinms;

auto getDefaultMuscle() {
    const vector<DoubleT> act_x{ 0.4241,
        0.4641,
        0.5441,
        0.6241,
        0.7041,
        0.7441,
        0.8641,
        0.9441,
        0.9841,
        1.0641,
        1.1041,
        1.2241,
        1.7841,
        1.8241 };
    const vector<DoubleT> act_y{ 0,
        0.0036,
        0.2531,
        0.5362,
        0.7739,
        0.8204,
        0.926,
        0.9919,
        1,
        0.9907,
        0.9524,
        0.7902,
        0.0117,
        0 };

    const vector<DoubleT> pas_x{ 1, 1.116, 1.22, 1.324, 1.428, 1.61, 2 };
    const vector<DoubleT> pas_y{ 0, 0.0208, 0.0604, 0.1536, 0.3117, 0.7448, 1.8571 };

    const vector<DoubleT> vel_x{
        -1, -0.624, -0.312, -0.104, 0, 0.104, 0.312, 0.624, 0.832, 0.988
    };
    const vector<DoubleT> vel_y{
        0, 0.0939, 0.3174, 0.6162, 1, 1.2278, 1.2972, 1.3507, 1.3823, 1.4
    };

    const vector<DoubleT> ten_x{ 0,
        0.01,
        0.02,
        0.06,
        0.1,
        0.14,
        0.18,
        0.24,
        0.32,
        0.38,
        0.42,
        0.44,
        0.46,
        0.48,
        0.5 };
    const vector<DoubleT> ten_y{ 0,
        0.094,
        0.2609,
        1.3087,
        2.4311,
        3.5536,
        4.676,
        6.3597,
        8.6046,
        10.2883,
        11.4107,
        11.9719,
        12.5332,
        13.0944,
        13.6556 };

    const CurveOffline act(act_x, act_y);
    const CurveOffline pas(pas_x, pas_y);
    const CurveOffline vel(vel_x, vel_y);
    const CurveOffline ten(ten_x, ten_y);

    ceinms::Lloyd2019Muscle::Parameters p;
    p.damping = 0.1;
    p.maxContractionVelocity = 1;
    p.maxIsometricForce = 1;
    p.optimalFiberLength = 1;
    p.pennationAngleAtOptimalFiberLength = 0;
    p.percentageChange = 0.15;
    p.strengthCoefficient = 1;
    p.tendonSlackLength = 1;
    p.activeForceLengthCurve = act;
    p.passiveForceLengthCurve = pas;
    p.forceVelocityCurve = vel;
    p.tendonForceStrainCurve = ten;
    ceinms::Lloyd2019Muscle muscle(p);
    muscle.setName("muscle");
    return muscle;
}

auto getDefaultActivation() {
    ceinms::ExponentialActivation act;
    act.setName("activation");
    return act;
}

int testConnections() {
    try {
        /*Test if automatic connection between input and component and between component and
         * component works well. The matching is based on naming.
         */
        cout << "------------TEST 1------------\n";
        NMSmodel model;
        ceinms::Lloyd2019Muscle muscle(getDefaultMuscle());
        muscle.setName("mtu1");
        ceinms::ExponentialActivation act(getDefaultActivation());
        act.setName("mtu1");
        model.addComponent(act);
        model.addComponent(muscle);
        model.addInput<MusculotendonLength>("mtu1");
        model.addInput<Excitation>("mtu1");

        model.connect<Excitation, ceinms::ExponentialActivation>();
        model.connect<MusculotendonLength, ceinms::Lloyd2019Muscle>();
        model.connect<ceinms::ExponentialActivation, ceinms::Lloyd2019Muscle>();
    } catch (...) { return 1; }
    try {
        /*Test if automatic connection between input and component and between component and
         * component works well. The matching is based on naming. Create a model with multiple input
         * and components
         */
        cout << "------------TEST 2------------\n";
        NMSmodel model;
        ceinms::Lloyd2019Muscle muscle(getDefaultMuscle());
        ceinms::ExponentialActivation act(getDefaultActivation());
        for (int i(0); i < 10; ++i) {
            string name = "mtu" + to_string(i);
            muscle.setName(name);
            act.setName(name);
            model.addComponent(muscle);
            model.addComponent(act);
            model.addInput<MusculotendonLength>(name);
            model.addInput<Excitation>(name);
        }
        cout << "# Connect Excitation input to ExponentialActivation component" << endl;
        model.connect<Excitation, ceinms::ExponentialActivation>();
        cout << "# Connect MusculotendonLength input to Lloyd2019Muscle component" << endl;
        model.connect<MusculotendonLength, ceinms::Lloyd2019Muscle>();
        cout << "# Connect ExponentialActivation component to Lloyd2019Muscle component" << endl;
        model.connect<ceinms::ExponentialActivation, ceinms::Lloyd2019Muscle>();
    } catch (...) { return 1; }
    try {
        /*Test if automatic connection between input and component and between component and
         * component works well. In this test we are going to create a model with a different number of input and components.
         */
        cout << "------------TEST 3------------\n";
        NMSmodel model;
        ceinms::Lloyd2019Muscle muscle(getDefaultMuscle());
        ceinms::ExponentialActivation act(getDefaultActivation());
        for (int i(0); i < 10; ++i) {
            string name = "mtu" + to_string(i);
            muscle.setName(name);
            act.setName(name);
            model.addComponent(muscle);
            model.addComponent(act);
            model.addInput<MusculotendonLength>(name);
            model.addInput<Excitation>(name);
        }
        for (int i(0); i < 3; ++i) {
            string name = "to_not_connect" + to_string(i);
            model.addInput<Excitation>(name);
        }
        cout << "# Connect Excitation input to ExponentialActivation component" << endl;
        model.connect<Excitation, ceinms::ExponentialActivation>();
        cout << "# Connect MusculotendonLength input to Lloyd2019Muscle component" << endl;
        model.connect<MusculotendonLength, ceinms::Lloyd2019Muscle>();
        cout << "# Connect ExponentialActivation component to Lloyd2019Muscle component" << endl;
        model.connect<ceinms::ExponentialActivation, ceinms::Lloyd2019Muscle>();
    } catch (...) { return 1; }
    try {
        /*Test if conncetion with multi input multi output component
         */
        cout << "------------TEST 4------------\n";
        NMSmodel model;
        ceinms::ExponentialActivation act(getDefaultActivation());
        MultiInputMultiOutput<Excitation, Excitation> mimo(2,3);
        mimo.setName("mimo");
        model.addComponent<MultiInputMultiOutput<Excitation, Excitation>>(mimo);
        for (int i(0); i < 3; ++i) {
            string name = "mtu" + to_string(i);
            act.setName(name);
            model.addComponent(act);
            model.addInput<Excitation>(name);
        }
        for (int i(0); i < 3; ++i) {
            string name = "to_not_connect" + to_string(i);
            model.addInput<Excitation>(name);
        }
        cout << "# Connect Excitation input to MultiInputMultiOutput component" << endl;
        model.connect<Excitation, MultiInputMultiOutput<Excitation, Excitation>>(
            "mtu0", { "mimo", 0 });
        model.connect<Excitation, MultiInputMultiOutput<Excitation, Excitation>>(
            "mtu1", { "mimo", 1 });


        cout << "# Connect MultiInputMultiOutput component to ExponentialActivation" << endl;
        model.connect<MultiInputMultiOutput<Excitation, Excitation>, ceinms::ExponentialActivation>(
            { "mimo", 0 }, "mtu0" );
        model.connect<MultiInputMultiOutput<Excitation, Excitation>, ceinms::ExponentialActivation>(
            { "mimo", 1 }, "mtu1");
        model.connect<MultiInputMultiOutput<Excitation, Excitation>, ceinms::ExponentialActivation>(
           { "mimo", 2 }, "mtu2");

    } catch (...) { return 1; }

    return 0;
}

int testConcepts() { 
    static_assert(is_same<typename Excitation::concept_t, input_t>::value);
    static_assert(is_same<typename MusculotendonLength::concept_t, input_t>::value);
    static_assert(is_same<typename MomentArm::concept_t, input_t>::value);
    static_assert(is_same<typename ceinms::ExponentialActivation::concept_t, component_t>::value);
    static_assert(is_same<typename ceinms::Lloyd2019Muscle::concept_t, component_t>::value);
    static_assert(is_same<typename Stage<ceinms::Lloyd2019Muscle>::concept_t, stage_t>::value);
    static_assert(is_same<typename Stage<ceinms::ExponentialActivation>::concept_t, stage_t>::value);
    static_assert(is_same<typename Source<Excitation>::concept_t, source_t>::value);
    static_assert(is_same<typename Source<MusculotendonLength>::concept_t, source_t>::value);
    static_assert(
        is_same<typename MultiInputMultiOutput<Excitation, MusculotendonLength>::concept_t,
            component_mimo_t>::value);
    return 0;
}

int main() {
    testConcepts();
    testConnections();
    return 0;
}