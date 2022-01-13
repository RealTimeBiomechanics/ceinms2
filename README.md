# ceinms2
![Build Status](https://travis-ci.org/RealTimeBiomechanics/ceinms2.svg?branch=master)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/0dc038f806e54ed6a04cd4838feed0dd)](https://app.codacy.com/gh/RealTimeBiomechanics/ceinms2?utm_source=github.com&utm_medium=referral&utm_content=RealTimeBiomechanics/ceinms2&utm_campaign=Badge_Grade)

CEINMS2 is a software library implementing a modular and computationally efficient neuromusculoskeletal model. 

## Windows installation

CEINMS2 uses Micorsoft's VCPKG as package manage and to install dependencies and CMake as cross-platform building environment.

[Install CMake](https://cmake.org/)

[Install VCPKG](https://github.com/microsoft/vcpkg)

After installing VCPKG, create an environment variable called `VCPKG_ROOT` pointing to the install directory of vcpkg, this will enable CMake to automatically find all the packages installed on your system using VCPKG. 

Install boost libraries
```sh
./vcpkg install boost-math:x64-windows
```

## Linux installation

to be completed

## CEINMS2 architecture

CEINMS2 architecture is inspired by [OpenSim](https://github.com/opensim-org/opensim-core), but reducing the highly complex hierarchical structure of OpenSim's models to the bare minimum. CEINMS2 is designed to be computationally efficient and numerous optimizations happen at compile-time. 

Architectural requirements:
* Must be scalable and enable any type and number of subcomponents
* Must not restrict the user to a predefined set of input
* Must be thread safe
* Should avoid run time polymorphism
* Should be computationally efficient

### The `NMSmodel` class

Example of how to create a new model:

```cpp

#include <ceinms2.h>
using namespace ceinms;

int main() {
    //Create a model that implements two computational stages,
    // an exponential activation dynamics and Lloyd2003 muscle model
    NMSmodel<ExponentialActivation,Lloyd2003Muscle> model;
    
    // Add a new muscle to the model
    Lloyd2003Muscle muscle(getDefaultMuscle());
    muscle.setName("mtu1");
    model.addComponent(muscle);

    // Add a new activation dynamic to the model
    ExponentialActivation act(getDefaultActivation());
    act.setName("mtu1");
    model.addComponent(act);
    
    // add inputs to the model
    model.addInput<MusculotendonLength>("mtu1");
    model.addInput<Excitation>("mtu1");

    // Connect Input to computational Stages
    model.connect<Excitation, ExponentialActivation>();
    model.connect<MusculotendonLength, Lloyd2003Muscle>();
    model.connect<ExponentialActivation, Lloyd2003Muscle>();
    return 0;
}
```

### Computational stages

### Components

#### Muscle spindle
CEINMS implements a computational model of the muscle spindle as described in Mileusnic, M. P., I. E. Brown, N. Lan and G. E. Loeb (2006). 
"Mathematical models of proprioceptors. I. Control and transduction in the muscle spindle." J Neurophysiol 96(4): 1772-1788. The following is a validation of the implemented model, which replicates the tests in Mileusnic et al., where dots indicate experimental data.

Code example
```cpp

ceinms::Mileusnic2006MuscleSpindle spindle;
std::vector<double> stretchVelocities = { 0.11, 0.66, 1.55 }; // relative to normalized length
std::vector<std::pair<double, double>> stimulations = { { 0., 0. }, { 0., 70. }, { 70., 0. } };
double dt = 0.001;
const double startLength = 0.95;
const double stopLength = 1.08;
const double rampStartTime = 1.;
for (const auto& stim : stimulations) {
    for (const auto& vel : stretchVelocities) {
        double t = 0.;
        double l = startLength;
        for (int i = 0; i < 4000; ++i) {
            if (t > rampStartTime) l = startLength + vel * (t - rampStartTime);
            if (l > stopLength) l = stopLength;
            spindle.setNormalizedMuscleFiberLength(l);
            spindle.setFusimotorFrequency(stim.first, stim.second);
            spindle.evaluate(dt);
            // results are stored in ceinms::Mileusnic2006MuscleSpindle::Output structure
            auto output = spindle.getOutput();
    }
}
```

![Muscle spindle validation during ramp tests](/doc/fig/Mileusnic2006MuscleSpindle_rampStretches.png)



### Data sources





