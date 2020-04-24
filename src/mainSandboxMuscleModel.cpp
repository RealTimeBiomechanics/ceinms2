#include <ceinms2/Lloyd2019Muscle.h>
#include <ceinms2/ExponentialActivation.h>

#include <vector>
#include <map>
#include <memory>
#include <tuple>
#include <functional>
#include <utility>
#include <algorithm>
#include <stdexcept>
#include <algorithm>
#include <string>
using namespace std;
using DoubleT = double;

enum class DataType { Excitation, MusculotendonLength, MomentArm };

template<DataType>
struct Input {
    using type = DataType;
    static constexpr std::string_view class_name = "Input";
    Input()
        : value(0) {}
    explicit Input(DoubleT val)
        : value(val) {}
    DoubleT value;
};

using Excitation = Input<DataType::Excitation>;
using MusculotendonLength = Input<DataType::MusculotendonLength>;
using MomentArm = Input<DataType::MomentArm>;

class Socket {
  public:
    Socket(char name[])
        : name_(name)
        , hasSlot_(false)
        , slot_(0) {}
    Socket(const string &name)
        : name_(name)
        , hasSlot_(false)
        , slot_(0) {}

    Socket(const string &name, size_t slot)
        : name_(name)
        , hasSlot_(true)
        , slot_(slot) {}
    [[nodiscard]] string getName() const { return name_; }
    [[nodiscard]] size_t getSlot() const { return slot_; }
    [[nodiscard]] bool hasSlot() const { return hasSlot_; }

  private:
    string name_;
    size_t slot_;
    bool hasSlot_;
};


auto getDefaultMuscle() {
    vector<DoubleT> act_x{ 0.4241,
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
    vector<DoubleT> act_y{ 0,
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

    vector<DoubleT> pas_x{ 1, 1.116, 1.22, 1.324, 1.428, 1.61, 2 };
    vector<DoubleT> pas_y{ 0, 0.0208, 0.0604, 0.1536, 0.3117, 0.7448, 1.8571 };

    vector<DoubleT> vel_x{
        -1, -0.624, -0.312, -0.104, 0, 0.104, 0.312, 0.624, 0.832, 0.988
    };
    vector<DoubleT> vel_y{
        0, 0.0939, 0.3174, 0.6162, 1, 1.2278, 1.2972, 1.3507, 1.3823, 1.4
    };

    vector<DoubleT> ten_x{ 0,
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
    vector<DoubleT> ten_y{ 0,
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

    CurveOffline act(act_x, act_y);
    CurveOffline pas(pas_x, pas_y);
    CurveOffline vel(vel_x, vel_y);
    CurveOffline ten(ten_x, ten_y);

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

template<typename, typename>
class MultiInputMultiOutput;

auto getDefaultActivation() {
    ceinms::ExponentialActivation act;
    act.setName("activation");
    return act;
}


void connectSocket(const ceinms::ExponentialActivation &parent,
    ceinms::Lloyd2019Muscle &child) {
    child.setActivation(parent.getOutput().activation);
}

void connectSocket(const MusculotendonLength &parent,
    ceinms::Lloyd2019Muscle &child) {
    child.setMusculotendonLength(parent.value);
}

void connectSocket(const Excitation &parent,
    ceinms::ExponentialActivation &child) {
    child.setInput(parent.value);
}


template<typename MultiInMultiOutT, typename Component, 
enable_if_t<!is_same<DataType, typename Component::type>::value, int> = 0>
void connectSocket(const MultiInMultiOutT &parent,
    Component &child,
    size_t parentSlot) {
    connectSocket(parent.getOutput().at(parentSlot), child);
}

template<typename InputT,
    typename MultiInMultiOutT,
    enable_if_t<is_same<DataType, typename InputT::type>::value, int> = 0>
     void connectSocket(const InputT &parent, MultiInMultiOutT &child,
    size_t childSlot) {
    connectSocket(parent, child.getInput().at(childSlot));
}


template<typename T>
class Stage {
  private:
    vector<unique_ptr<T>> components_;
    map<string, vector<function<void(T &)>>> connections_;

  public:
    using type = T;
    /*Makes an internal copy of `component`*/
    void addComponent(const T &component) noexcept { components_.emplace_back(new T(component)); }

    vector<DoubleT> evaluate(DoubleT dt) noexcept {
        vector<DoubleT> ret;
        for (const auto &c : components_) {
            for (auto &f : connections_[c->getName()])
                f(*c);
            c->evaluate(dt);
            ret.emplace_back(c->getOutput().getPrimary());
        }
        return ret;
    }

    template<typename U>
    void connectToParent(Socket childSocket, const U &parent) noexcept {
        struct Command {
            Command(const U &parent, Socket childSocket)
                : p(parent)
                , s(childSocket) {}
            const U &p;
            const Socket s;
            void operator()(T &child) {
                if (s.hasSlot())
                    connectSocket(p, child);
                else
                    connectSocket<U, T>(p, child, s.getSlot());
            }
        };
        connections_[childSocket.getName()].emplace_back(
            function<void(T &)>(Command(parent, childSocket)));
    }


    template<typename U>
    void connectToParent(string childName, const U &parent, size_t subslot) noexcept {
        struct Command {
            Command(const U &parent)
                : p(parent) {}
            const U &p;
            void operator()(T &child) { connectSocket(p, child); }
        };
        connections_[childName].emplace_back(function<void(T &)>(Command(parent)));
    }


    [[nodiscard]] const T &get(string name) const {
        auto it = find_if(cbegin(components_), cend(components_), [&name](const auto &e) {
            return e->getName() == name;
        });
        if (it == cend(components_)) throw invalid_argument("Component " + name + " not found");
        return **it;
    }

    [[nodiscard]] T &get(string name) {
        auto it = find_if(begin(components_), end(components_), [&name](const auto &e) {
            return e->getName() == name;
        });
        if (it == end(components_)) throw invalid_argument("Component " + name + " not found");
        return **it;
    }

    [[nodiscard]] auto getOutput() const noexcept {
        vector<typename T::Output> res;
        for (const auto &e : components_)
            res.emplace_back(e->getOutput());
        return res;
    }

    auto getNames() const noexcept {
        vector<string> names;
        for (const auto &c : components_)
            names.emplace_back(c->getName());
        return names;
    }
};


template<typename T>
class Source {
  private:
    vector<T> input_;
    vector<string> names_;
    size_t nInput_;

  public:
    using type = T;
    Source()
        : nInput_(0) {}

    void addInput(string name) noexcept {
        names_.emplace_back(name);
        input_.emplace_back(0.);
        ++nInput_;
    }

    void set(vector<T> values) noexcept { input_ = values; }

    void setInputNames(const vector<string>& names) noexcept {
        names_ = names;
        nInput_ = names.size();
    }

    [[nodiscard]] T &get(string name) {
        auto it = find(begin(names_), end(names_), name);
        if (it == end(names_))
            throw invalid_argument("Input " + name + " not found");
        size_t i = distance(begin(names_), it);
        return input_.at(i);
    }

    [[nodiscard]] const T &get(string name) const {
        auto it = find(cbegin(names_), cend(names_), name);
        if (it == cend(names_))
            throw invalid_argument("Input " + name + " not found");
        size_t i = distance(cbegin(names_), it);
        return input_.at(i);
    }

    auto getNames() const noexcept { return names_; }

    size_t getSize() const noexcept { return nInput_; }
};

template<typename InT, typename OutT>
class MultiInputMultiOutput {
  private:
    using type = nullptr_t;
    using in_type = InT;
    using out_type = OutT;
    function<vector<OutT>(const vector<InT>&)> fun_;
    size_t nInput_, nOutput_;
    vector<InT> input_;
    vector<OutT> output_;
    string name_;

  public:
    MultiInputMultiOutput(size_t nInput, size_t nOutput)
        : nInput_(nInput)
        , nOutput_(nOutput)
        , fun_([](const vector<OutT> &) { return vector<InT>{};
        }) {
        input_.resize(nInput_);
        output_.resize(nOutput_);
    }
    void setFunction(function<vector<OutT>(const vector<InT> &)> fun) { fun_ = fun; }
    void setInput(const vector<InT> &input);
    const vector<OutT> &getOutput() const { return output_; }
    vector<OutT> &getInput() { return input_; }

    void setName(string name) { name_ = name; }
    void evaluate(DoubleT) { output_ = fun_(input_); }
};


class NMSmodel {
  public:

    template<typename T>
    void addComponent(const T &component) noexcept;

    template<typename T>
    void addInput(const string &name) noexcept;

    template<typename T>
    [[nodiscard]] auto getOutput() const noexcept;

    template<typename Parent, typename Child>
    void connect(Socket parent, Socket child);

    template<typename Parent, typename Child>
    void connect();

    template<typename T>
    void setInput(const vector<T> &inputVector) noexcept;

    template<typename T>
    [[nodiscard]] const auto &getComponent(string name) const;

     template<typename T>
    [[nodiscard]] auto &getComponent(string name);

    void evaluate(DoubleT dt) noexcept;

  private:
    tuple<Stage<ceinms::ExponentialActivation>,
        Stage<ceinms::Lloyd2019Muscle>,
        Stage<MultiInputMultiOutput<Excitation, Excitation>>>
        stages_;
    tuple<Source<Excitation>, Source<MusculotendonLength>, Source<MomentArm>>
        input_;

    //`T` and `U` are of type `Stage`
    template<typename T,
        typename U,
        enable_if_t<is_same<T, Stage<typename T::type>>::value
                        && is_same<U, Stage<typename U::type>>::value,
            int> = 0>
    void connect_() {
        auto parentComponentNames = std::get<T>(stages_).getNames();
        auto childComponentNames = std::get<U>(stages_).getNames();

        std::sort(parentComponentNames.begin(), parentComponentNames.end());
        std::sort(childComponentNames.begin(), childComponentNames.end());

        std::vector<string> commonNames;

        std::set_intersection(parentComponentNames.begin(),
            parentComponentNames.end(),
            childComponentNames.begin(),
            childComponentNames.end(),
            std::back_inserter(commonNames));
        for (const string &name : commonNames) {
            connect<typename T::type, typename U::type>(name, name);
        }
        
    }

        //`T` is an Input and and `U` is a `Stage`
    template<typename T,
        typename U,
        enable_if_t<is_same<DataType, typename T::type>::value
                        && is_same<U, Stage<typename U::type>>::value,
            int> = 0>
    void connect_() {
        auto parentComponentNames = std::get<Source<T>>(input_).getNames();
        auto childComponentNames = std::get<U>(stages_).getNames();

        std::sort(parentComponentNames.begin(), parentComponentNames.end());
        std::sort(childComponentNames.begin(), childComponentNames.end());

        std::vector<string> commonNames;

        std::set_intersection(parentComponentNames.begin(),
            parentComponentNames.end(),
            childComponentNames.begin(),
            childComponentNames.end(),
            std::back_inserter(commonNames));
        for (const string &name : commonNames) {
            connect<T, typename U::type>(name, name);
        }
    }

    // Need to implement component concepts
    //`T` is of type `Source` and `U` is a Component
    template<typename T,
        typename U,
        enable_if_t<is_same<DataType, typename T::type>::value, int> = 0>
    void connect_(Socket parent, Socket child) {
        cout << "Connecting " << T::class_name << "." << parent.getName() << " -> " << U::class_name
             << "." << child.getName() << endl;

        const auto &input = get<Source<T>>(input_).get(parent.getName());
        std::get<Stage<U>>(stages_).connectToParent(child.getName(), input);
    }


    // Need to implement the concept of component, once that is available
    //`T` and `U` are both a Component
    template<typename T,
        typename U,
        enable_if_t<!is_same<DataType, typename T::type>::value
                        && !is_same<DataType, typename U::type>::value
                        && !is_same<T, Stage<typename T::type>>::value
                        && !is_same<U, Stage<typename U::type>>::value,
            int> = 0>
    void connect_(Socket parent, Socket child, ...) {
        cout << "Connecting  " << T::class_name<< "." << parent.getName()
             << " -> " << U::class_name << "." << child.getName() << endl;
        auto parentComponent = get<Stage<T>>(stages_).get(parent.getName());
        std::get<Stage<U>>(stages_).connectToParent(child.getName(), parentComponent);
    }

    template<typename T,
        typename U,
        enable_if_t<is_same<DataType, typename U::type>::value, int> = 0>
    void connect_(Socket parent, Socket child, ...) {
        static_assert(false, "It is not possible to connect to an Input");
    }
};

template<typename T>
void NMSmodel::setInput(const std::vector<T> &values) noexcept {
    //TODO: check size first
    std::get<Source<T>>(input_).set(values);
}

template<typename T>
void NMSmodel::addInput(const string &name) noexcept {
    //TODO: check for other input with same name
    // throw error if same component already exists
    std::get<Source<T>>(input_).addInput(name);
}

template<typename T>
void NMSmodel::addComponent(const T &component) noexcept {
    // TODO: check for other components with same name
    //throw error if same component already exists
    std::get<Stage<T>>(stages_).addComponent(component);
}

template<typename Parent, typename Child>
void NMSmodel::connect(Socket parent, Socket child) {
    connect_<Parent, Child>(parent, child);
}

template<typename Parent, typename Child>
void NMSmodel::connect() {
    if constexpr (is_same<DataType, typename Parent::type>::value)
        connect_<Parent, Stage<Child>>();
    else
        connect_<Stage<Parent>, Stage<Child>>();
}

void NMSmodel::evaluate(DoubleT dt) noexcept {
    std::get<Stage<ceinms::ExponentialActivation>>(stages_).evaluate(dt);
    std::get<Stage<ceinms::Lloyd2019Muscle>>(stages_).evaluate(dt);
}

template<typename T>
auto NMSmodel::getOutput() const noexcept {
    return std::get<Stage<T>>(stages_).getOutput();
}

template<typename Component>
const auto& NMSmodel::getComponent(string name) const {
    return std::get<Stage<Component>>(stages_).get(name);
}

template<typename Component>
auto& NMSmodel::getComponent(string name){
    return std::get<Stage<Component>>(stages_).get(name);
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
        cout << "Connect Excitation input to ExponentialActivation component" << endl;
        model.connect<Excitation, ceinms::ExponentialActivation>();
        cout << "Connect MusculotendonLength input to Lloyd2019Muscle component" << endl;
        model.connect<MusculotendonLength, ceinms::Lloyd2019Muscle>();
        cout << "Connect ExponentialActivation component to Lloyd2019Muscle component" << endl;
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
        cout << "Connect Excitation input to ExponentialActivation component" << endl;
        model.connect<Excitation, ceinms::ExponentialActivation>();
        cout << "Connect MusculotendonLength input to Lloyd2019Muscle component" << endl;
        model.connect<MusculotendonLength, ceinms::Lloyd2019Muscle>();
        cout << "Connect ExponentialActivation component to Lloyd2019Muscle component" << endl;
        model.connect<ceinms::ExponentialActivation, ceinms::Lloyd2019Muscle>();
    } catch (...) { return 1; }
    return 0;
}

int main() {
    testConnections();
    return 0;
}