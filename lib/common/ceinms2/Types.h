#ifndef ceinms2_Types_h
#define ceinms2_Types_h
#include <string_view>
#include <concepts>
#include <vector>

template<typename T, typename U>
concept can_add = requires(T t, U u) {
    t + u;
};


namespace ceinms {
using DoubleT = double;
enum class Concept { Data, Component, Component_MIMO, Stage, Source };
template<Concept>
struct Selector;
using data_t = Selector<Concept::Data>;
using component_t = Selector<Concept::Component>;
using component_mimo_t = Selector<Concept::Component_MIMO>;
using stage_t = Selector<Concept::Stage>;
using source_t = Selector<Concept::Source>;

enum class DataType { Excitation, MusculotendonLength, MomentArm, Activation, Force, Torque, Frequency, NormalizedFiberLength };

template<DataType T>
constexpr std::string_view getDataTypeName() {
    if constexpr (T == DataType::Excitation) return "Excitation";
    if constexpr (T == DataType::MusculotendonLength) return "MusculotendonLength";
    if constexpr (T == DataType::MomentArm) return "MomentArm";
    if constexpr (T == DataType::Activation) return "Activation";
    if constexpr (T == DataType::Force) return "Force";
    if constexpr (T == DataType::Torque) return "Torque";
    if constexpr (T == DataType::Frequency) return "Frequency";
    if constexpr (T == DataType::NormalizedFiberLength) return "NormalizedFiberLength";
    return "Undefined";
}

template<DataType U, typename T = DoubleT>
class Data {
  private:
    T value_{ 0. };

  public:
    using concept_t = data_t;
    using type = DataType;
    static constexpr std::string_view class_name{getDataTypeName<U>() };
    Data() = default;
    Data(T val)
        : value_(val) {}
    [[nodiscard]] T &get() { return value_; }
    [[nodiscard]] const T &get() const { return value_; }
    operator T() const { return get(); }
};
using Excitation = Data<DataType::Excitation>;
using MusculotendonLength = Data<DataType::MusculotendonLength>;
using MomentArmsOnDof = Data<DataType::MomentArm, std::vector<DoubleT>>;
using MomentArm = Data<DataType::MomentArm>;
using Activation = Data<DataType::Activation>;
using Force = Data<DataType::Force>;
using Torque = Data<DataType::Torque>;
using Frequency = Data<DataType::Frequency>;
using NormalizedFiberLength = Data<DataType::NormalizedFiberLength>;


template<typename T>
concept ForceGenerator = requires(T m) {
    m.getOutput().force;
};

template<typename T>
concept ActivationGenerator = requires(T m) {
    m.getOutput().activation;
};

template<typename T>
concept ExcitationGenerator = requires(T m) {
    m.getOutput().excitation;
};



}// namespace ceinms
#endif