#include <ceinms2/Lloyd2019Muscle.h>
#include <ceinms2/ExponentialActivation.h>
#include <ceinms2/NMSmodel.h>
#include <iostream>

using namespace std;
using std::cout;
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

    const vector<DoubleT> vel_x{ -1, -0.624, -0.312, -0.104, 0, 0.104, 0.312, 0.624, 0.832, 0.988 };
    const vector<DoubleT> vel_y{
        0, 0.0939, 0.3174, 0.6162, 1, 1.2278, 1.2972, 1.3507, 1.3823, 1.4
    };

    const vector<DoubleT> ten_x{
        0, 0.01, 0.02, 0.06, 0.1, 0.14, 0.18, 0.24, 0.32, 0.38, 0.42, 0.44, 0.46, 0.48, 0.5
    };
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

namespace ceinms {

template<int N>
class BypassComponent {
  public:
    struct Output {
        DoubleT output{ 0 };
        DoubleT getPrimary() const { return output; }
    };

  private:
    DoubleT excitation_{ 0 }, activation_{ 0 }, musculotendonLength_{ 0 };
    Output o_;
    string name_;

  public:
    static constexpr int value = N;
    using concept_t = component_t;
    static constexpr string_view class_name = "BypassComponent";
    BypassComponent() = default;
    [[nodiscard]] std::string getName() const { return name_; }
    void setName(std::string name) { name_ = name; }
    void setInput(Excitation input) { excitation_ = input.value; }
    void setInput(MusculotendonLength input) { musculotendonLength_ = input.value; }
    void setInput(Activation input) { activation_ = input.value; }
    void evaluate(DoubleT) {
        if constexpr (value == 0) {
            o_.output = excitation_ + 1;
            cout << "input excitation: " << excitation_ << " output: " << o_.output << endl;
        } else if constexpr (value == 1) {
            o_.output = activation_ + 10;
            cout << "input activation: " << activation_ << " output: " << o_.output << endl;
        } else if constexpr (value == 2) {
            o_.output = musculotendonLength_ + 100;
            cout << "input musculotendonLength: " << musculotendonLength_
                 << " output: " << o_.output << endl;
        } else {
            cout << "Error, bad connection" << endl;
            o_.output = 0;
        }
    }
    const Output &getOutput() const { return o_; }
};

void connectSocket(const Excitation &parent, BypassComponent<0> &child) {
    child.setInput(parent);
}

void connectSocket(const Activation &parent, BypassComponent<1> &child) {
    child.setInput(parent);
}

void connectSocket(const MusculotendonLength &parent, BypassComponent<2> &child) {
    child.setInput(parent);
}

void connectSocket(const BypassComponent<0> &parent, BypassComponent<1> &child) {
    child.setInput(Activation(parent.getOutput().getPrimary()));
}
}// namespace ceinms

int testConnectionFlow() {
    cout << "#### TEST CONNECTION FLOW ####\n";
    using ExcitationBypass = BypassComponent<0>;
    using ActivationBypass = BypassComponent<1>;
    using MTLBypass = BypassComponent<2>;

    try {
        cout << "------------TEST 1.1-----------\n";
        cout << "Connecting a single input to a single component and testing the correctness of "
                "the output\n";
        using MyNMSmodel = NMSmodel<ExcitationBypass>;
        MyNMSmodel model;
        ExcitationBypass exct;
        exct.setName("muscle1");
        model.addComponent<ExcitationBypass>(exct);
        model.addInput<Excitation>("muscle1");
        model.connect<Excitation, ExcitationBypass>();
        model.setInput(vector{ Excitation{ 1 } });
        model.evaluate(0.1);
        if (model.getOutput<ExcitationBypass>().at(0).getPrimary() != 2)
            throw std::logic_error("Wrong connection flow\n");
        model.setInput(vector{ Excitation{ 2 } });
        model.evaluate(0.1);
        if (model.getOutput<ExcitationBypass>().at(0).getPrimary() != 3)
            throw std::logic_error("Wrong connection flow\n");

    } catch (const std::exception &e) {
        std::cout << e.what();
        return 1;
    }
    try {
        cout << "------------TEST 1.2-----------\n";
        cout << "Connecting a single input to a single component and testing the correctness of "
                "the output\n";
        using MyNMSmodel = NMSmodel<ExcitationBypass, ActivationBypass>;
        MyNMSmodel model;
        ExcitationBypass exct;
        ActivationBypass act;
        exct.setName("muscle1");
        act.setName("muscle1");
        model.addComponent<ExcitationBypass>(exct);
        model.addComponent<ActivationBypass>(act);
        model.addInput<Excitation>("muscle1");
        model.connect<Excitation, ExcitationBypass>();
        model.connect<ExcitationBypass, ActivationBypass>();
        model.setInput(vector{ Excitation{ 1 } });
        model.evaluate(0.1);
        if (model.getOutput<ActivationBypass>().at(0).getPrimary() != 12)
            throw std::logic_error("Wrong connection flow\n");
        model.setInput(vector{ Excitation{ 3 } });
        model.evaluate(0.1);
        if (model.getOutput<ActivationBypass>().at(0).getPrimary() != 14)
            throw std::logic_error("Wrong connection flow\n");

    } catch (const std::exception &e) {
        std::cout << e.what();
        return 1;
    }
    try {
        cout << "------------TEST 1.3-----------\n";
        cout << "Connecting a single input to a single component and testing the correctness of "
                "the output\n";
        using MyNMSmodel = NMSmodel<ExcitationBypass, ActivationBypass, MTLBypass>;
        MyNMSmodel model;
        ExcitationBypass exct;
        ActivationBypass act;
        MTLBypass mtl;
        exct.setName("muscle1");
        act.setName("muscle1");
        mtl.setName("muscle1");
        model.addComponent<ExcitationBypass>(exct);
        model.addComponent<ActivationBypass>(act);
        model.addComponent<MTLBypass>(mtl);
        model.addInput<Excitation>("muscle1");
        model.addInput<MusculotendonLength>("muscle1");
        model.connect<Excitation, ExcitationBypass>();
        model.connect<ExcitationBypass, ActivationBypass>();
        model.connect<MusculotendonLength, MTLBypass>();
        model.setInput(vector{ Excitation{ 1 } });
        model.setInput(vector{ MusculotendonLength{ 1.1 } });
        model.evaluate(0.1);
        if (model.getOutput<ActivationBypass>().at(0).getPrimary() != 12)
            throw std::logic_error("Wrong connection flow\n");
        if (model.getOutput<MTLBypass>().at(0).getPrimary() != 101.1)
            throw std::logic_error("Wrong connection flow\n");
    } catch (const std::exception &e) {
        std::cout << e.what();
        return 1;
    }
    try {
        cout << "------------TEST 2------------\n";
        cout << "Connecting a single input to multiple components and testing the correctness of "
                "the output\n";
        using MyNMSmodel = NMSmodel<ExcitationBypass, ActivationBypass, MTLBypass>;
        MyNMSmodel model;
        ExcitationBypass exct;
        ActivationBypass act;
        MTLBypass mtl;
        exct.setName("muscle1");
        model.addComponent<ExcitationBypass>(exct);
        act.setName("muscle1a");
        model.addComponent<ActivationBypass>(act);
        act.setName("muscle1b");
        model.addComponent<ActivationBypass>(act);
        mtl.setName("muscle1a");
        model.addComponent<MTLBypass>(mtl);
        mtl.setName("muscle1b");
        model.addComponent<MTLBypass>(mtl);

        model.addInput<Excitation>("muscle1");
        model.addInput<MusculotendonLength>("muscle1a");
        model.addInput<MusculotendonLength>("muscle1b");

        model.connect<Excitation, ExcitationBypass>();
        model.connect<ExcitationBypass, ActivationBypass>("muscle1", "muscle1a");
        model.connect<ExcitationBypass, ActivationBypass>("muscle1", "muscle1b");
        model.connect<MusculotendonLength, MTLBypass>();
        model.setInput(vector{ Excitation{ 1 } });
        model.setInput(vector{ MusculotendonLength{ 1.1 }, MusculotendonLength{ 3.1 } });
        model.evaluate(0.1);
        auto out1 = model.getOutput<ActivationBypass>();
        if (out1.at(0).getPrimary() != 12 && out1.at(1).getPrimary())
            throw std::logic_error("Wrong connection flow\n");
        auto out2 = model.getOutput<MTLBypass>();
        if (out2.at(0).getPrimary() != 101.1 && out2.at(1).getPrimary() != 103.1)
            throw std::logic_error("Wrong connection flow\n");
    } catch (const std::exception &e) {
        std::cout << e.what();
        return 1;
    }
    return 0;
}

int testConnections() {
    using MyNMSmodel = NMSmodel<ExponentialActivation,
        Lloyd2019Muscle,
        MultiInputMultiOutput<Excitation, Excitation>>;
    try {
        /*Test if automatic connection between input and component and between component and
         * component works well. The matching is based on naming.
         */
        cout << "------------TEST 1------------\n";
        MyNMSmodel model;
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
    } catch (const std::exception &e) {
        std::cout << e.what();
        return 1;
    }
    try {
        /*Test if automatic connection between input and component and between component and
         * component works well. The matching is based on naming. Create a model with multiple input
         * and components
         */
        cout << "------------TEST 2------------\n";
        MyNMSmodel model;
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
    } catch (const std::exception &e) {
        std::cout << e.what();
        return 1;
    }
    try {
        /*Test if automatic connection between input and component and between component and
         * component works well. In this test we are going to create a model with a different number
         * of input and components.
         */
        cout << "------------TEST 3------------\n";
        MyNMSmodel model;
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
    } catch (const std::exception &e) {
        std::cout << e.what();
        return 1;
    }
    try {
        /*Test if conncetion with multi input multi output component
         */
        cout << "------------TEST 4------------\n";
        MyNMSmodel model;
        ceinms::ExponentialActivation act(getDefaultActivation());
        MultiInputMultiOutput<Excitation, Excitation> mimo(2, 3);
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
            { "mimo", 0 }, "mtu0");
        model.connect<MultiInputMultiOutput<Excitation, Excitation>, ceinms::ExponentialActivation>(
            { "mimo", 1 }, "mtu1");
        model.connect<MultiInputMultiOutput<Excitation, Excitation>, ceinms::ExponentialActivation>(
            { "mimo", 2 }, "mtu2");

    } catch (const std::exception &e) {
        std::cout << e.what();
        return 1;
    }

    return 0;
}

int testConcepts() {
    static_assert(is_same<typename Excitation::concept_t, input_t>::value);
    static_assert(is_same<typename MusculotendonLength::concept_t, input_t>::value);
    static_assert(is_same<typename MomentArm::concept_t, input_t>::value);
    static_assert(is_same<typename ceinms::ExponentialActivation::concept_t, component_t>::value);
    static_assert(is_same<typename ceinms::Lloyd2019Muscle::concept_t, component_t>::value);
    static_assert(is_same<typename Stage<ceinms::Lloyd2019Muscle>::concept_t, stage_t>::value);
    static_assert(
        is_same<typename Stage<ceinms::ExponentialActivation>::concept_t, stage_t>::value);
    static_assert(is_same<typename Source<Excitation>::concept_t, source_t>::value);
    static_assert(is_same<typename Source<MusculotendonLength>::concept_t, source_t>::value);
    static_assert(
        is_same<typename MultiInputMultiOutput<Excitation, MusculotendonLength>::concept_t,
            component_mimo_t>::value);
    return 0;
}

template<typename T>
bool runTest(T testFun, string name) {
    bool failed = testFun();
    cout << name << ": ";
    if (failed)
        cout << "FAILED\n";
    else
        cout << "PASSED\n";
    return failed;
}

int main() {
    runTest(&testConcepts, "      TestConcepts");
    runTest(&testConnections, "   TestConnections");
    runTest(&testConnectionFlow, "TestConnectionFlow");
    return 0;
}