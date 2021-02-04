# ceinms2

[![Codacy Badge](https://api.codacy.com/project/badge/Grade/0dc038f806e54ed6a04cd4838feed0dd)](https://app.codacy.com/gh/RealTimeBiomechanics/ceinms2?utm_source=github.com&utm_medium=referral&utm_content=RealTimeBiomechanics/ceinms2&utm_campaign=Badge_Grade)

![Build Status](https://travis-ci.org/RealTimeBiomechanics/ceinms2.svg?branch=master)

CEINMS2 is a software library implementing a modular and computationally efficient neuromusculoskeletal model. 

## Windows installation

CEINMS2 uses Micorsoft's VCPKG as package manage and to install dependencies and CMake as cross-platform building environment.

[Install CMake](https://cmake.org/)

[Install VCPKG](https://github.com/microsoft/vcpkg)

After installing VCPKG, create an environment variable called `VCPKG_ROOT` pointing to the install directory of vcpkg, this will enable CMake to automatically find all the packages installed on your system using VCPKG. 

Install boost libraries
```
./vcpkg install boost-math:x64-windows
```

## Linux installation

to be completed

## CEINMS2 architecture

CEINMS2 architecture is inspired by [OpenSim](https://github.com/opensim-org/opensim-core), but reducing the highly complex hierarchical structure of OpenSim's models to the bare minimum. CEINMS2 is designed to be computationally efficient and numerous optimizations happen at compile-time. 

Architectural requirements:
* Must be scalable and enable any type and number of subcomponents
* Must not restrict the user to a predefined set of input
* Must be computationally efficient
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

### Data sources





