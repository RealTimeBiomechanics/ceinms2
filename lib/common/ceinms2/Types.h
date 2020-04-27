#ifndef ceinms2_Types_h
#define ceinms2_Types_h
#include <string_view>

using DoubleT = double;
enum class Concept {Input, Component, Component_MIMO, Stage, Source};
template<Concept>
struct Selector;
using input_t = Selector<Concept::Input>;
using component_t = Selector<Concept::Component>;
using component_mimo_t = Selector<Concept::Component_MIMO>;
using stage_t =  Selector<Concept::Stage>;
using source_t = Selector<Concept::Source>;

enum class DataType { Excitation, MusculotendonLength, MomentArm, Activation, Force };
template<DataType>
struct Input {
    using concept_t = input_t;
    using type = DataType;
    static constexpr std::string_view class_name = "Input";
    Input() = default;
    explicit Input(DoubleT val)
        : value(val) {}
    DoubleT value{0.};
};
using Excitation = Input<DataType::Excitation>;
using MusculotendonLength = Input<DataType::MusculotendonLength>;
using MomentArm = Input<DataType::MomentArm>;
using Activation = Input<DataType::Activation>;
using Force = Input<DataType::Force>;

#endif