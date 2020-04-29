#ifndef ceinms2_Types_h
#define ceinms2_Types_h
#include <string_view>

using DoubleT = double;
enum class Concept {Data, Component, Component_MIMO, Stage, Source};
template<Concept>
struct Selector;
using data_t = Selector<Concept::Data>;
using component_t = Selector<Concept::Component>;
using component_mimo_t = Selector<Concept::Component_MIMO>;
using stage_t =  Selector<Concept::Stage>;
using source_t = Selector<Concept::Source>;

enum class DataType { Excitation, MusculotendonLength, MomentArm, Activation, Force };
template<DataType>
class Data {
  private:
    DoubleT value_{0.};
  public:
    using concept_t = data_t;
    using type = DataType;
    static constexpr std::string_view class_name = "Input";
    Data() = default;
    explicit Data(DoubleT val)
        : value_(val) {}
    [[nodiscard]] DoubleT &get() { return value_; }
    [[nodiscard]] const DoubleT &get() const { return value_; }
};
using Excitation = Data<DataType::Excitation>;
using MusculotendonLength = Data<DataType::MusculotendonLength>;
using MomentArm = Data<DataType::MomentArm>;
using Activation = Data<DataType::Activation>;
using Force = Data<DataType::Force>;

#endif