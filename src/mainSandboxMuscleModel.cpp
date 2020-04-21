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
using namespace std;
using DoubleT = double;

enum class DataType {
    Excitation,
    MusculotendonLength,
    MomentArm
};

template<DataType>
struct Input {
    Input()
        : value(0) {}
    explicit Input(DoubleT val)
        : value(val) {}
    DoubleT value;
};

using Excitation = Input<DataType::Excitation>;
using MusculotendonLength = Input<DataType::MusculotendonLength>;
using MomentArm = Input<DataType::MomentArm>;


auto getDefaultMuscle() {
    vector<DoubleT> act_x{ 0.4241, 0.4641, 0.5441, 0.6241, 0.7041, 0.7441, 0.8641, 0.9441, 0.9841, 1.0641, 1.1041, 1.2241, 1.7841, 1.8241 };
    vector<DoubleT> act_y{ 0, 0.0036, 0.2531, 0.5362, 0.7739, 0.8204, 0.926, 0.9919, 1, 0.9907, 0.9524, 0.7902, 0.0117, 0 };

    vector<DoubleT> pas_x{ 1, 1.116, 1.22, 1.324, 1.428, 1.61, 2 };
    vector<DoubleT> pas_y{ 0, 0.0208, 0.0604, 0.1536, 0.3117, 0.7448, 1.8571 };

    vector<DoubleT> vel_x{ -1, -0.624, -0.312, -0.104, 0, 0.104, 0.312, 0.624, 0.832, 0.988 };
    vector<DoubleT> vel_y{ 0, 0.0939, 0.3174, 0.6162, 1, 1.2278, 1.2972, 1.3507, 1.3823, 1.4 };

    vector<DoubleT> ten_x{ 0, 0.01, 0.02, 0.06, 0.1, 0.14, 0.18, 0.24, 0.32, 0.38, 0.42, 0.44, 0.46, 0.48, 0.5 };
    vector<DoubleT> ten_y{ 0, 0.094, 0.2609, 1.3087, 2.4311, 3.5536, 4.676, 6.3597, 8.6046, 10.2883, 11.4107, 11.9719, 12.5332, 13.0944, 13.6556 };

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


auto getDefaultActivation() {
    ceinms::ExponentialActivation act;
    act.setName("activation");
    return act;
}


void connectSocket(const ceinms::ExponentialActivation& parent, ceinms::Lloyd2019Muscle& child) {
    child.setActivation(parent.getOutput().activation);
}


void connectSocket(const MusculotendonLength &parent, ceinms::Lloyd2019Muscle &child) {
    child.setMusculotendonLength(parent.value);
}

void connectSocket(const Excitation &parent, ceinms::ExponentialActivation &child) {
    child.setInput(parent.value);
}

template<typename T>
class Stage {
  private:
    vector<unique_ptr<T>> components_;
    map<string, vector<function<void(T &)>>> connections_;

  public:
    /*Makes an internal copy of `component`*/
    void addComponent(const T &component) {
        components_.emplace_back(new T(component));
    }

    vector<DoubleT> evaluate(DoubleT dt) {
        vector<DoubleT> ret;
        for (const auto &c : components_) {
            for (auto &f : connections_[c->getName()])
                f(*c);
            ret.emplace_back(c->evaluate(dt));
        }
        return ret;
    }

    template<typename U>
    void connectToParent(string childName, const U &parent) {
        struct Command {
            Command(const U &parent)
                : p(parent) {}
            const U &p;
            void operator()(T &child) {
                connectSocket(p, child);
            }
        };
        connections_[childName].emplace_back(function<void(T &)>(Command(parent)));
    }

    [[nodiscard]] const T &get(string name) const {
        auto it = find_if(cbegin(components_), cend(components_), [&name](const auto &e) {
            return e->getName() == name;
        });
        if (it == cend(components_))
            throw invalid_argument("Component " + name + " not found");
        return **it;
    }

    [[nodiscard]] T &get(string name) {
        auto it = find_if(begin(components_), end(components_), [&name](const auto &e) {
            return e->getName() == name;
        });
        if (it == end(components_))
            throw invalid_argument("Component " + name + " not found");
        return **it;
    }

    [[nodiscard]] auto getOutput() const {
        vector<typename T::Output> res;
        for (const auto &e : components_)
            res.emplace_back(e->getOutput());
        return res;
    }
};


template<typename T>
class Source {
  private:
    vector<T> input_;
    vector<string> names_;
    size_t nInput_;

  public:
    Source()
        : nInput_(0) {}

    void addInput(string name) {
        names_.emplace_back(name);
        input_.emplace_back(0.);
        ++nInput_;
    }

    void set(vector<T> values) {
        input_ = values;
    }

    void setInputNames(vector<string> names) {
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

    auto getNames() const {
        return names_;
    }

    size_t getNInput() const {
        return nInput_;
    }
};


struct NMSmodel {
    template<typename T>
    void addComponent(const T &component);

    template<typename T>
    void addInput(const string& name);

    template<typename T>
    [[nodiscard]] auto getOutput() const;

    template<typename Child, typename Parent>
    void connectChildToParent(string childName, string parentName);

    template<typename Component, typename Input>
    void connectToInput(string componentName, string inputName);

    template<typename T>
    void setInput(const vector<T> &inputVector);

    template<typename T>
    [[nodiscard]] const auto getComponent(string name) const ;

    DoubleT evaluate(DoubleT dt);
    tuple<Stage<ceinms::ExponentialActivation>, Stage<ceinms::Lloyd2019Muscle>> stages_;
    tuple<Source<Excitation>, Source<MusculotendonLength>, Source<MomentArm>> input_;
};

template<typename T>
void NMSmodel::setInput(const std::vector<T>& values) {
    //check size first
    std::get<Source<T>>(input_).set(values);
}

template<typename T>
void NMSmodel::addInput(const string &name) {
    std::get<Source<T>>(input_).addInput(name);
}

template<typename T>
void NMSmodel::addComponent(const T& component) {
    std::get<Stage<T>>(stages_).addComponent(component);
}

template<typename Child, typename Parent>
void NMSmodel::connectChildToParent(string childName, string parentName) {
    auto parent = get<Stage<Parent>>(stages_).get(parentName);
    std::get<Stage<Child>>(stages_).connectToParent(childName, parent);

}

template<typename Component, typename Input>
void NMSmodel::connectToInput(string componentName, string inputName) {

    const auto& input = get<Source<Input>>(input_).get(inputName);
    std::get<Stage<Component>>(stages_).connectToParent(componentName, input);
}

DoubleT NMSmodel::evaluate(DoubleT dt) {
 
    std::get<Stage<ceinms::ExponentialActivation>>(stages_).evaluate(dt);
    std::get<Stage<ceinms::Lloyd2019Muscle>>(stages_).evaluate(dt);
    return 0;
}

template<typename T>
auto NMSmodel::getOutput() const  {
    return std::get<Stage<T>>(stages_).getOutput();
}


template<typename Component>
const auto NMSmodel::getComponent(string name) const {
    return std::get<Stage<Component>>(stages_).get(name);
}


int main_() {

    ceinms::Lloyd2019Muscle muscle = getDefaultMuscle();
    ceinms::ExponentialActivation act = getDefaultActivation();
    NMSmodel myModel;
    myModel.addComponent(muscle);
    myModel.addComponent(act);
  //  myModel.connectChildToParent<ceinms::Lloyd2019Muscle, ceinms::ExponentialActivation>("muscle", "activation");

    act.setInput(1);
    auto a = act.evaluate(0.1);
    muscle.setActivation(a);
    muscle.setMusculotendonLength(1.);
  //  auto force = muscle.evaluate(0.1);
    return 0;
}

int main() {
    ceinms::Lloyd2019Muscle muscle(getDefaultMuscle());
    //ceinms::ExponentialActivation act = getDefaultActivation();
    ceinms::ExponentialActivation act(getDefaultActivation());
    NMSmodel model;
    model.addComponent(act);
    model.addComponent(muscle);
    model.addInput<MusculotendonLength>("mtl1");
    model.addInput<Excitation>("excitation1");
    model.connectChildToParent<ceinms::Lloyd2019Muscle, ceinms::ExponentialActivation>("muscle", "activation");
    model.connectToInput<ceinms::Lloyd2019Muscle, MusculotendonLength>("muscle", "mtl1");
    model.connectToInput<ceinms::ExponentialActivation, Excitation>("activation", "excitation1");

        
    Excitation neuralInput{ 0.6 };
    MusculotendonLength musculotendonLengthInput{ 2 };

    model.setInput(vector<Excitation>{ neuralInput });
    model.setInput(vector<MusculotendonLength>{ musculotendonLengthInput});
    model.evaluate(0.01);
   for (auto &e : model.getOutput<ceinms::ExponentialActivation>())
        cout << e.activation << endl;
    return 0;
}