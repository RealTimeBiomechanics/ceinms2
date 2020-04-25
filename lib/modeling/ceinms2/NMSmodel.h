#ifndef ceinms2_NMSmodel_h
#define ceinms2_NMSmodel_h
#include <string>
#include <string_view>
#include <vector>
#include <iostream>
#include <tuple>
#include <functional>
#include <type_traits>
#include <memory>
#include <map>
#include "ceinms2/Types.h"
#include "ceinms2/ExponentialActivation.h"
#include "ceinms2/Lloyd2019Muscle.h"

namespace ceinms {

class Socket;

template<typename T>
class Stage;

template<typename T>
class Source;

template<typename InT, typename OutT>
class MultiInputMultiOutput;

class NMSmodel;

template<typename MultiInMultiOutT,
    typename Component,
    std::enable_if_t<std::is_same<typename MultiInMultiOutT::concept_t, component_mimo_t>::value
                    && std::is_same<typename Component::concept_t, component_t>::value,
        int> = 0>
void connectSocket(const MultiInMultiOutT &parent, Component &child, size_t parentSlot) {
    connectSocket(parent.getOutput().at(parentSlot), child);
}

template<typename InputT,
    typename MultiInMultiOutT,
    std::enable_if_t<std::is_same<DataType, typename InputT::type>::value
                    && std::is_same<typename MultiInMultiOutT::concept_t, component_mimo_t>::value,
        int> = 0>
void connectSocket(const InputT &parent, MultiInMultiOutT &child, size_t childSlot) {
    connectSocket(parent, child.getInput().at(childSlot));
}

template<typename T, typename U>
void connectSocket(const T &, U &, ...) {
    cout << "No socket connection" << endl;
}

class Socket {
  public:
    Socket(char name[])
        : name_(name)
        , hasSlot_(false)
        , slot_(0) {}
    Socket(const std::string &name)
        : name_(name)
        , hasSlot_(false)
        , slot_(0) {}

    Socket(const std::string &name, size_t slot)
        : name_(name)
        , hasSlot_(true)
        , slot_(slot) {}
    [[nodiscard]] std::string getName() const { return name_; }
    [[nodiscard]] size_t getSlot() const { return slot_; }
    [[nodiscard]] bool hasSlot() const { return hasSlot_; }

  private:
    std::string name_;
    size_t slot_;
    bool hasSlot_;
};

std::ostream &operator<<(std::ostream &os, const Socket &sock) {
    os << sock.getName();
    if (sock.hasSlot()) { os << "." << sock.getSlot(); }
    return os;
}

template<typename T>
class Stage {
  private:
    std::vector<std::unique_ptr<T>> components_;
    std::map<std::string, std::vector<std::function<void(T &)>>> connections_;

  public:
    using type = T;
    using concept_t = stage_t;
    /*Makes an internal copy of `component`*/
    void addComponent(const T &component) noexcept { components_.emplace_back(new T(component)); }

    std::vector<DoubleT> evaluate(DoubleT dt) noexcept {
        std::vector<DoubleT> ret;
        for (const auto &c : components_) {
            for (auto &f : connections_[c->getName()]) {
                f(*c);
            }
            c->evaluate(dt);
            ret.emplace_back(c->getOutput().getPrimary());
        }
        return ret;
    }

    /*Need to check what slot is used to connect to a parent with multi input multi output*/
    template<typename U>
    void connectToParent(Socket childSocket, const U &parent) noexcept {
        struct Command {
            Command(const U &parent)
                : p(parent) {}
            const U &p;
            void operator()(T &child) { connectSocket(p, child); }
        };

        struct CommandWithSlot {
            CommandWithSlot(const U &parent, Socket childSocket)
                : p(parent)
                , s(childSocket) {}
            const U &p;
            const Socket s;
            void operator()(T &child) { connectSocket<U, T>(p, child, s.getSlot()); }
        };

        if (!childSocket.hasSlot()) {
            connections_[childSocket.getName()].emplace_back(function<void(T &)>(Command(parent)));
        } else {
            connections_[childSocket.getName()].emplace_back(
                function<void(T &)>(CommandWithSlot(parent, childSocket)));
        }
    }

    [[nodiscard]] const T &get(std::string name) const {
        auto it = std::find_if(std::cbegin(components_), std::cend(components_), [&name](const auto &e) {
            return e->getName() == name;
        });
        if (it == std::cend(components_)) throw std::invalid_argument("Component " + name + " not found");
        return **it;
    }

    [[nodiscard]] T &get(std::string name) {
        auto it = std::find_if(std::begin(components_), std::end(components_), [&name](const auto &e) {
            return e->getName() == name;
        });
        if (it == std::end(components_)) { throw std::invalid_argument("Component " + name + " not found"); }
        return **it;
    }

    [[nodiscard]] auto getOutput() const noexcept {
        std::vector<typename T::Output> res;
        for (const auto &e : components_) {
            res.emplace_back(e->getOutput());
        }
        return res;
    }

    [[nodiscard]] auto getNames() const noexcept {
        std::vector<std::string> names;
        for (const auto &c : components_) {
            names.emplace_back(c->getName());
        }
        return names;
    }
};


template<typename T>
class Source {
  private:
    std::vector<T> input_;
    std::vector<std::string> names_;
    size_t nInput_{ 0 };

  public:
    using type = T;
    using concept_t = source_t;
    Source() = default;

    void addInput(std::string name) noexcept {
        names_.emplace_back(name);
        input_.emplace_back(0.);
        ++nInput_;
    }

    void set(std::vector<T> values) noexcept { input_ = values; }

    void setInputNames(const std::vector<std::string> &names) noexcept {
        names_ = names;
        nInput_ = names.size();
    }

    [[nodiscard]] T &get(std::string name) {
        auto it = std::find(begin(names_), std::end(names_), name);
        if (it == std::end(names_)) { throw std::invalid_argument("Input " + name + " not found"); }
        size_t i = distance(begin(names_), it);
        return input_.at(i);
    }

    [[nodiscard]] const T &get(std::string name) const {
        auto it = std::find(std::cbegin(names_), std::cend(names_), name);
        if (it == std::cend(names_)) throw std::invalid_argument("Input " + name + " not found");
        size_t i = std::distance(std::cbegin(names_), it);
        return input_.at(i);
    }

    [[nodiscard]] auto getNames() const noexcept { return names_; }

    [[nodiscard]] size_t getSize() const noexcept { return nInput_; }
};

template<typename InT, typename OutT>
class MultiInputMultiOutput {
  private:
    std::function<std::vector<OutT>(const std::vector<InT> &)> fun_;
    size_t nInput_, nOutput_;
    std::vector<InT> input_;
    std::vector<OutT> output_;
    std::string name_;

  public:
    using concept_t = component_mimo_t;
    static constexpr std::string_view class_name = "MultiInputMultiOutput";
    MultiInputMultiOutput(size_t nInput, size_t nOutput)
        : nInput_(nInput)
        , nOutput_(nOutput)
        , fun_([](const std::vector<OutT> &) { return std::vector<InT>{}; }) {
        input_.resize(nInput_);
        output_.resize(nOutput_);
    }
    void setFunction(std::function<std::vector<OutT>(const std::vector<InT> &)> fun) { fun_ = fun; }
    void setInput(const std::vector<InT> &input);
    [[nodiscard]] const std::vector<OutT> &getOutput() const { return output_; }
    std::vector<OutT> &getInput() { return input_; }

    void setName(std::string name) { name_ = name; }
    [[nodiscard]] std::string getName() const { return name_; }
    void evaluate(DoubleT) { output_ = fun_(input_); }
};


class NMSmodel {
  public:
    template<typename T>
    void addComponent(const T &component) noexcept;

    template<typename T>
    void addInput(const std::string &name) noexcept;

    template<typename T>
    [[nodiscard]] auto getOutput() const noexcept;

    template<typename Parent, typename Child>
    void connect(Socket parent, Socket child);

    template<typename Parent, typename Child>
    void connect();

    template<typename T>
    void setInput(const std::vector<T> &inputVector) noexcept;

    template<typename T>
    [[nodiscard]] const auto &getComponent(std::string name) const;

    template<typename T>
    [[nodiscard]] auto &getComponent(std::string name);

    void evaluate(DoubleT dt) noexcept;

  private:
    std::tuple<Stage<ceinms::ExponentialActivation>,
        Stage<ceinms::Lloyd2019Muscle>,
        Stage<MultiInputMultiOutput<Excitation, Excitation>>>
        stages_;
    std::tuple<Source<Excitation>, Source<MusculotendonLength>, Source<MomentArm>> input_;

    //`T` and `U` are of type `Stage`
    template<typename T,
        typename U,
        std::enable_if_t<std::is_same<typename T::concept_t, stage_t>::value
                        && std::is_same<typename U::concept_t, stage_t>::value,

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
        std::enable_if_t<std::is_same<typename T::concept_t, input_t>::value
                        && std::is_same<typename U::concept_t, stage_t>::value,
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

    //`T` is an `input_t` and `U` is a `component_t` or a `component_mimo_t`
    template<typename T,
        typename U,
        std::enable_if_t<std::is_same<typename T::concept_t, input_t>::value
                        && (std::is_same<typename U::concept_t, component_t>::value
                            || std::is_same<typename U::concept_t, component_mimo_t>::value),
            int> = 0>
    void connect_(Socket parent, Socket child) {
        cout << "Connecting " << T::class_name << "." << parent << " -> " << U::class_name << "."
             << child << endl;

        const auto &input = get<Source<T>>(input_).get(parent.getName());
        std::get<Stage<U>>(stages_).connectToParent(child, input);
    }


    //`T` and `U` are either a `component_t` or a `component_mimmo_t`
    template<typename T,
        typename U,
        std::enable_if_t<(std::is_same<typename T::concept_t, component_t>::value
                        || std::is_same<typename T::concept_t, component_mimo_t>::value)
                        && (std::is_same<typename U::concept_t, component_t>::value
                            || std::is_same<typename U::concept_t, component_mimo_t>::value),
            int> = 0>
    void connect_(Socket parent, Socket child) {
        cout << "Connecting " << T::class_name << "." << parent << " -> " << U::class_name << "."
             << child << endl;
        auto parentComponent = get<Stage<T>>(stages_).get(parent.getName());
        std::get<Stage<U>>(stages_).connectToParent(child, parentComponent);
    }

    template<typename T,
        typename U,
        std::enable_if_t<std::is_same<DataType, typename U::type>::value, int> = 0>
    void connect_(Socket parent, Socket child, ...) {
        static_assert(false, "It is not possible to connect to an Input");
    }
};


template<typename T>
void NMSmodel::setInput(const std::vector<T> &values) noexcept {
    // TODO: check size first
    std::get<Source<T>>(input_).set(values);
}

template<typename T>
void NMSmodel::addInput(const std::string &name) noexcept {
    // TODO: check for other input with same name
    // throw error if same component already exists
    std::get<Source<T>>(input_).addInput(name);
}

template<typename T>
void NMSmodel::addComponent(const T &component) noexcept {
    // TODO: check for other components with same name
    // throw error if same component already exists
    std::get<Stage<T>>(stages_).addComponent(component);
}

template<typename Parent, typename Child>
void NMSmodel::connect(Socket parent, Socket child) {
    connect_<Parent, Child>(parent, child);
}

template<typename Parent, typename Child>
void NMSmodel::connect() {
    if constexpr (is_same<DataType, typename Parent::type>::value) {
        connect_<Parent, Stage<Child>>();
    } else {
        connect_<Stage<Parent>, Stage<Child>>();
    }
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
const auto &NMSmodel::getComponent(std::string name) const {
    return std::get<Stage<Component>>(stages_).get(name);
}

template<typename Component>
auto &NMSmodel::getComponent(std::string name) {
    return std::get<Stage<Component>>(stages_).get(name);
}


}// namespace ceinms

#endif