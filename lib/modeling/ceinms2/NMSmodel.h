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
#include <algorithm>
#include "ceinms2/Types.h"

namespace ceinms {

class Socket;

template<typename T>
class Stage;

template<typename T>
class Source;

template<typename InT, typename OutT>
class MultiInputMultiOutput;

template<typename...>
class NMSmodel;

template<typename MultiInMultiOutT,
    typename Component,
    std::enable_if_t<std::is_same<typename MultiInMultiOutT::concept_t, component_mimo_t>::value
                         && std::is_same<typename Component::concept_t, component_t>::value,
        int> = 0>
void connectSocket(const MultiInMultiOutT &parent, Component &child, size_t parentSlot) {
    child.setInput(parent.getOutput(parentSlot));
}

template<typename InputT,
    typename MultiInMultiOutT,
    std::enable_if_t<
        std::is_same<typename InputT::concept_t, data_t>::value
            && std::is_same<typename MultiInMultiOutT::concept_t, component_mimo_t>::value,
        int> = 0>
void connectSocket(const InputT &parent, MultiInMultiOutT &child, size_t childSlot) {
    child.setInput(childSlot, parent);
}

template<typename T,
    typename U,
    std::enable_if_t<std::is_same<typename T::concept_t, data_t>::value
                         && std::is_same<typename U::concept_t, data_t>::value,
        int> = 0>
void connectSocket(const T &parent, U &child) {
    connectSocket(parent, child.get());
}

template<typename T, typename U>
void connectSocket(const T &, U &, ...) {
    std::cout << "No socket connection" << std::endl;
}


class Socket {
  public:
    Socket(const std::string& name)
        : name_(name) {}
    Socket(const std::string& name, size_t slot)
        : name_(name)
        , slot_(slot)
        , hasSlot_(true) {}
    Socket(const char* name)
        : Socket(std::string(name)) {}
    Socket(const char *name, size_t slot)
        : Socket(std::string(name), slot) {}

    [[nodiscard]] std::string getName() const { return name_; }
    [[nodiscard]] size_t getSlot() const { return slot_; }
    [[nodiscard]] bool hasSlot() const { return hasSlot_; }

  private:
    std::string name_;
    size_t slot_{ 0 };
    bool hasSlot_{ false };
};

std::ostream &operator<<(std::ostream &os, const Socket &sock) {
    os << sock.getName();
    if (sock.hasSlot()) { os << "." << sock.getSlot(); }
    return os;
}

template<typename T>
class Stage {
  private:
    std::vector<std::shared_ptr<T>> components_;
    std::map<std::string, std::vector<std::function<void(T&)>>> connections_;

  public:
    using type = T;
    using concept_t = stage_t;
    /*Makes an internal copy of `component`*/
    void addComponent(const T &component) noexcept { components_.emplace_back(new T(component)); }

    void evaluate(DoubleT dt) noexcept {
        for (const auto &c : components_) {
            for (auto &f : connections_[c->getName()]) {
                f(*c);
            }
            c->evaluate(dt);
        }
    }

    /*Need to check what slot is used to connect to a parent with multi input multi output*/
    template<typename U>
    void connectToParent(Socket childSocket, const std::shared_ptr<U> parent) noexcept {

        struct Command {
            Command(const std::shared_ptr<U> parent)
                : p(parent) {}
            const std::shared_ptr<U> p;
            void operator()(T &child) { connectSocket(*p, child); }
        };

        struct CommandWithSlot {
            CommandWithSlot(const std::shared_ptr<U> parent, Socket childSocket)
                : p(parent)
                , s(childSocket) {}
            const std::shared_ptr<U> p;
            const Socket s;
            void operator()(T &child) { connectSocket(*p, child, s.getSlot()); }
        };

        if (!childSocket.hasSlot()) {
            connections_[childSocket.getName()].emplace_back(std::function<void(T&)>(Command(parent)));
        } else {
            connections_[childSocket.getName()].emplace_back(
                std::function<void(T&)>(CommandWithSlot(parent, childSocket)));
        }
    }

    [[nodiscard]] const T &get(std::string name) const {
        auto idx{ getIndex(name) };
        return *components_.at(idx);
    }

    [[nodiscard]] T &get(std::string name) {
        auto idx{ getIndex(name) };
        return *components_.at(idx);
    }

    [[nodiscard]] auto getPtr(std::string name) {
        auto idx{ getIndex(name) };
        return components_.at(idx);
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

    private:
    [[nodiscard]] size_t getIndex(std::string name) const {
          auto it = std::find_if(std::cbegin(components_),
              std::cend(components_),
              [&name](const auto &e) { return e->getName() == name; });
          if (it == std::cend(components_))
              throw std::invalid_argument("Component " + name + " of type " + std::string(T::class_name) + " not found");
          return static_cast<size_t>(std::distance(std::cbegin(components_), it));
    }
};


template<typename T>
class Source {
  private:
    std::vector<std::shared_ptr<T>> input_;
    std::vector<std::string> names_;
    size_t nInput_{ 0 };

  public:
    using type = T;
    using concept_t = source_t;
    Source() = default;

    void addInput(std::string name) noexcept {
        names_.emplace_back(name);
        input_.emplace_back(std::make_shared<T>());
        ++nInput_;
    }

    void set(std::vector<T> values) { 
        if (nInput_ != values.size())
            throw std::invalid_argument(
                "Provided " + std::to_string(values.size()) + " input, but " + std::to_string(nInput_) + " were expected.");
        for (size_t i{ 0 }; i < nInput_; ++i) {
            *input_.at(i) = values.at(i);
        }
    }

    void setInputNames(const std::vector<std::string> &names) noexcept {
        names_ = names;
        nInput_ = names.size();
    }

    [[nodiscard]] T &get(std::string name) {
        size_t idx{ getIndex(name) };
        return *input_.at(idx);
    }

    [[nodiscard]] const T &get(std::string name) const {
        size_t idx{ getIndex(name) };
        return *input_.at(idx);
    }

    [[nodiscard]] auto getPtr(std::string name) const {
        size_t idx{ getIndex(name) };
        return input_.at(idx);
    }

    [[nodiscard]] auto getNames() const noexcept { return names_; }

    [[nodiscard]] size_t getSize() const noexcept { return nInput_; }

    private:
    [[nodiscard]] size_t getIndex(std::string name) const {
        auto it = std::find_if(std::cbegin(names_), std::cend(names_),
            [&name](const auto &e) { return e == name; });
        if (it == std::cend(names_))
            throw std::invalid_argument("Input " + name + " of type " + std::string(T::class_name) + " not found");
        return static_cast<size_t>(std::distance(std::cbegin(names_), it));
    }
};

template<typename InT, typename OutT>
class MultiInputMultiOutput {
  private:
    std::function<std::vector<OutT>(const std::vector<InT> &)> fun_{([](const std::vector<InT> &) { return std::vector<OutT>{}; })};
    size_t nInput_{0}, nOutput_{0};
    //consider having these input_ and output_ as shared_ptr, same as 
    // per `Source`
    std::vector<InT> input_;
    std::vector<OutT> output_;
    std::string name_;

  public:
    using concept_t = component_mimo_t;
    static constexpr std::string_view class_name = "MultiInputMultiOutput";
    MultiInputMultiOutput(size_t nInput, size_t nOutput)
        : nInput_(nInput)
        , nOutput_(nOutput) {
        input_.resize(nInput_);
        output_.resize(nOutput_);
    }
    void setFunction(std::function<std::vector<OutT>(const std::vector<InT> &)> fun) { fun_ = fun; }
    void setInput(size_t index, const InT& value) { input_.at(index) = value; }

    void setInput(const std::vector<InT>& input) {
        if (nInput_ != input.size())
            throw std::invalid_argument("Provided " + std::to_string(input.size()) + " input, but "
                                        + std::to_string(nInput_) + " were expected.");
        input_ = input;
    }

    [[nodiscard]] const std::vector<OutT> &getOutput() const { return output_; }
    [[nodiscard]] OutT getOutput(size_t index) const { return output_.at(index); }

    void setName(std::string name) { name_ = name; }
    [[nodiscard]] std::string getName() const { return name_; }
    void evaluate(DoubleT) { output_ = fun_(input_); }
};

template<typename... Args>
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
    void setInput(const std::vector<T> &inputVector);

    template<typename T>
    [[nodiscard]] const auto &getComponent(std::string name) const;

    template<typename T>
    [[nodiscard]] auto &getComponent(std::string name);

    void evaluate(DoubleT dt) noexcept;

  private:
    std::tuple<Stage<Args>...> stages_;
    // This can be improved so that the type of the sources is automatically derived from the input
    // types of computational stages
    std::tuple<Source<Excitation>, Source<Activation>, Source<MusculotendonLength>, Source<MomentArmsOnDof>> sources_;

    //`T` and `U` are of type `Stage`
    template<typename T,
        typename U,
        std::enable_if_t<std::is_same<typename T::concept_t, stage_t>::value
                             && std::is_same<typename U::concept_t, stage_t>::value,

            int> = 0>
    void connectStages() {
        auto parentComponentNames = std::get<T>(stages_).getNames();
        auto childComponentNames = std::get<U>(stages_).getNames();

        std::sort(parentComponentNames.begin(), parentComponentNames.end());
        std::sort(childComponentNames.begin(), childComponentNames.end());

        std::vector<std::string> commonNames;

        std::set_intersection(parentComponentNames.begin(),
            parentComponentNames.end(),
            childComponentNames.begin(),
            childComponentNames.end(),
            std::back_inserter(commonNames));
        for (const std::string &name : commonNames) {
            connectComponentToComponent<typename T::type, typename U::type>(name, name);
        }
    }

    //`T` is an Input and and `U` is a `Stage`
    template<typename T,
        typename U,
        std::enable_if_t<std::is_same<typename T::concept_t, data_t>::value
                             && std::is_same<typename U::concept_t, stage_t>::value,
            int> = 0>
    void connectInputToStages() {
        auto parentComponentNames = std::get<Source<T>>(sources_).getNames();
        auto childComponentNames = std::get<U>(stages_).getNames();

        std::sort(parentComponentNames.begin(), parentComponentNames.end());
        std::sort(childComponentNames.begin(), childComponentNames.end());

        std::vector<std::string> commonNames;

        std::set_intersection(parentComponentNames.begin(),
            parentComponentNames.end(),
            childComponentNames.begin(),
            childComponentNames.end(),
            std::back_inserter(commonNames));
        for (const std::string &name : commonNames) {
            connectInputToComponent<T, typename U::type>(name, name);
        }
    }

    //`T` is an `data_t` and `U` is a `component_t` or a `component_mimo_t`
    template<typename T,
        typename U,
        std::enable_if_t<std::is_same<typename T::concept_t, data_t>::value
                             && (std::is_same<typename U::concept_t, component_t>::value
                                 || std::is_same<typename U::concept_t, component_mimo_t>::value),
            int> = 0>
    void connectInputToComponent(Socket parent, Socket child) {
        std::cout << "Connecting " << T::class_name << "." << parent << " -> " << U::class_name
                  << "." << child << std::endl;

        auto input = std::get<Source<T>>(sources_).getPtr(parent.getName());
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
    void connectComponentToComponent(Socket parent, Socket child) {
        std::cout << "Connecting " << T::class_name << "." << parent << " -> " << U::class_name
                  << "." << child << std::endl;
        auto parentComponent = std::get<Stage<T>>(stages_).getPtr(parent.getName());
        std::get<Stage<U>>(stages_).connectToParent(child, parentComponent);
    }

    template<typename T,
        typename U,
        std::enable_if_t<std::is_same<typename T::concept_t, data_t>::value, int> = 0>
    void connect_(Socket parent, Socket child, ...) {
        //   static_assert(false, "It is not possible to connect to an Input");
    }
};

template<typename... Args>
template<typename T>
void NMSmodel<Args...>::setInput(const std::vector<T> &values) {
    // TODO: check size first
    std::get<Source<T>>(sources_).set(values);
}

template<typename... Args>
template<typename T>
void NMSmodel<Args...>::addInput(const std::string &name) noexcept {
    // TODO: check for other input with same name
    // throw error if same component already exists
    std::get<Source<T>>(sources_).addInput(name);
}

template<typename... Args>
template<typename T>
void NMSmodel<Args...>::addComponent(const T &component) noexcept {
    // TODO: check for other components with same name
    // throw error if same component already exists
    std::get<Stage<T>>(stages_).addComponent(component);
}

template<typename... Args>
template<typename Parent, typename Child>
void NMSmodel<Args...>::connect(Socket parent, Socket child) {
    if constexpr (std::is_same<typename Parent::concept_t, data_t>::value) {
        connectInputToComponent<Parent, Child>(parent, child);
    } else {
        connectComponentToComponent<Parent, Child>(parent, child);
    }
}


template<typename... Args>
template<typename Parent, typename Child>
void NMSmodel<Args...>::connect() {
    if constexpr (std::is_same<typename Parent::concept_t, data_t>::value) {
        connectInputToStages<Parent, Stage<Child>>();
    } else {
        connectStages<Stage<Parent>, Stage<Child>>();
    }
}

//Evaluate all the `Stage`s in the order in which they were defined
template<typename... Args>
void NMSmodel<Args...>::evaluate(DoubleT dt) noexcept {
    std::apply([dt](auto &&... args) { (args.evaluate(dt), ...); }, stages_);

}

template<typename... Args>
template<typename T>
auto NMSmodel<Args...>::getOutput() const noexcept {
    return std::get<Stage<T>>(stages_).getOutput();
}

template<typename... Args>
template<typename Component>
const auto &NMSmodel<Args...>::getComponent(std::string name) const {
    return std::get<Stage<Component>>(stages_).get(name);
}

template<typename... Args>
template<typename Component>
auto &NMSmodel<Args...>::getComponent(std::string name) {
    return std::get<Stage<Component>>(stages_).get(name);
}


}// namespace ceinms

#endif
