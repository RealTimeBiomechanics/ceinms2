#include "ceinms2/Mileusnic2006MuscleSpindle.h"
#include "ceinms2/WDBsolver.h"
#include <cmath>
#include <boost/numeric/odeint.hpp>

namespace ceinms {

void Mileusnic2006IntrafusalFiber::setActivation(DoubleT fStatic, DoubleT fDynamic) {
    i_.fStatic = fStatic;
    i_.fDynamic = fDynamic;
}

void Mileusnic2006IntrafusalFiber::setNormalizedMuscleFiberLength(
    DoubleT normalizedMuscleFiberLength) {
    i_.normalisedMuscleFiberLength = normalizedMuscleFiberLength;
}

void Mileusnic2006IntrafusalFiber::setInput(Activation fStatic, Activation fDynamic) {
    i_.fStatic = fStatic.get();
    i_.fDynamic = fDynamic.get();
}

void Mileusnic2006IntrafusalFiber::setInput(Input input) {
    i_ = input;
}

DoubleT Mileusnic2006IntrafusalFiber::integratePolarRegionLength(DoubleT dt) {
    DoubleT minPolarRegionLength = 0.;
    DoubleT maxPolarRegionLength = 1.;

    auto f = [&](double polarRegionLength) -> double {
        const auto polarRegionVelocity =
            Mileusnic2006IntrafusalFiber::calculatePolarRegionVelocityFromPolarRegionLength(
                this->s_.polarRegionLength, polarRegionLength, dt);
        const auto forceSensoryRegion = Mileusnic2006IntrafusalFiber::calculateSensoryRegionForce(
            this->i_.normalisedMuscleFiberLength,
            polarRegionLength,
            this->p_.sensoryRegionStiffness,
            this->p_.sensoryRegionRestLength);

        const auto beta = calculateBeta(
            this->i_.fDynamic, this->i_.fStatic, this->p_.beta0, this->p_.beta1, this->p_.beta2);

        DoubleT coefficientOfAsymmetry = p_.coefficientOfAsymmetryLengthening;
        if (polarRegionVelocity < 0.) 
            coefficientOfAsymmetry = p_.coefficientOfAsymmetryShortening;
        const auto fvpr = calculateVelocityContributionToPolarRegionForce(polarRegionLength,
            polarRegionVelocity,
            beta,
            coefficientOfAsymmetry,
            this->p_.muscleFascicleSlackLength,
            this->p_.velocityPowerTerm);
        const auto fspr = calculateStiffnessContributionToPolarRegionForce(
            polarRegionLength, this->p_.polarRegionStiffness, this->p_.polarRegionRestLength);
        const auto gamma =
            calculateGamma(this->i_.fDynamic, this->i_.fStatic, this->p_.gamma1, this->p_.gamma2);
        const auto forcePolarRegion = fvpr + fspr + gamma;
        const auto err = forceSensoryRegion - forcePolarRegion;
        return err;
    };

    const auto currentLengthPolarRegion =
        wdbSolve(f, minPolarRegionLength, maxPolarRegionLength, 1e-9);
    return currentLengthPolarRegion;
}

void Mileusnic2006IntrafusalFiber::integrate(DoubleT dt) {
    sNew_.polarRegionLength = integratePolarRegionLength(dt);
    sNew_.polarRegionVelocity = calculatePolarRegionVelocityFromPolarRegionLength(
        s_.polarRegionLength, sNew_.polarRegionLength, dt);
}

void Mileusnic2006IntrafusalFiber::equilibrate() {
    const double tol = 1e-9;
    do {
        integrate(0.001);
        validateState();
    } while (sNew_.polarRegionVelocity > tol);
}

void Mileusnic2006IntrafusalFiber::validateState() {
    s_ = sNew_;
}

void Mileusnic2006IntrafusalFiber::evaluate(DoubleT dt) {
    integrate(dt);
    calculateOutput();
    validateState();
}

DoubleT Mileusnic2006IntrafusalFiber::calculateVelocityContributionToPolarRegionForce(
    DoubleT polarRegionLength,
    DoubleT polarRegionVelocity,
    DoubleT beta,
    DoubleT coefficientOfAsymmetry,
    DoubleT muscleFascicleSlackLength,
    DoubleT velocityPowerTerm) {
    return beta * coefficientOfAsymmetry * (polarRegionLength - muscleFascicleSlackLength)
           * std::signbit(polarRegionVelocity)
           * std::pow(std::abs(polarRegionVelocity), velocityPowerTerm);
}

DoubleT Mileusnic2006IntrafusalFiber::calculateStiffnessContributionToPolarRegionForce(
    DoubleT polarRegionLength,
    DoubleT polarRegionStiffness,
    DoubleT polarRegionRestLength) {
    return polarRegionStiffness * (polarRegionLength - polarRegionRestLength);
}

DoubleT Mileusnic2006IntrafusalFiber::calculateSensoryRegionForce(
    DoubleT normalisedMuscleFiberLength,
    DoubleT polarRegionLength,
    DoubleT sensoryRegionStiffness,
    DoubleT sensoryRegionRestLength) {
    return sensoryRegionStiffness
           * ((normalisedMuscleFiberLength - polarRegionLength) - sensoryRegionRestLength);
}

DoubleT Mileusnic2006IntrafusalFiber::calculateBeta(DoubleT fDynamic,
    DoubleT fStatic,
    DoubleT beta0,
    DoubleT beta1,
    DoubleT beta2) {
    return beta0 + beta1 * fDynamic + beta2 * fStatic;
}

DoubleT Mileusnic2006IntrafusalFiber::calculateGamma(DoubleT fDynamic,
    DoubleT fStatic,
    DoubleT gamma1,
    DoubleT gamma2) {
    return gamma1 * fDynamic + gamma2 * fStatic;
}

DoubleT Mileusnic2006IntrafusalFiber::calculatePrimaryAfferent() const {
    auto const tension = calculateSensoryRegionForce(i_.normalisedMuscleFiberLength,
        sNew_.polarRegionLength,
        p_.sensoryRegionStiffness,
        p_.sensoryRegionRestLength);
    return p_.sensoryRegionStretchToPrimaryAfferentFiring
           * (tension / p_.sensoryRegionStiffness
               - (p_.sensoryRegionThresholdLength - p_.sensoryRegionRestLength));
}

DoubleT Mileusnic2006IntrafusalFiber::calculateSecondaryAfferent() const {
    auto const tension = calculateSensoryRegionForce(i_.normalisedMuscleFiberLength,
        sNew_.polarRegionLength,
        p_.sensoryRegionStiffness,
        p_.sensoryRegionRestLength);
    const auto primaryAfferent = calculatePrimaryAfferent();
    const auto term1 = p_.percentageSecondaryAfferentOnSensoryRegion
                       * p_.secondaryAfferentRestLength / p_.sensoryRegionRestLength
                       * primaryAfferent;
    const auto term2 = p_.sensoryRegionStretchToSecondaryAfferentFiring
                       * (1. - p_.percentageSecondaryAfferentOnSensoryRegion)
                       * p_.secondaryAfferentRestLength / p_.polarRegionRestLength
                       * (i_.normalisedMuscleFiberLength - tension / p_.sensoryRegionStiffness
                           - p_.sensoryRegionRestLength - p_.polarRegionThresholdLength);
    return term1 + term2;
}

void Mileusnic2006IntrafusalFiber::calculateOutput() {
    o_.beta = calculateBeta(i_.fDynamic, i_.fStatic, p_.beta0, p_.beta1, p_.beta2);
    o_.gamma = calculateGamma(i_.fDynamic, i_.fStatic, p_.gamma1, p_.gamma2);
    o_.fiberTension = calculateSensoryRegionForce(i_.normalisedMuscleFiberLength,
        sNew_.polarRegionLength,
        p_.sensoryRegionStiffness,
        p_.sensoryRegionRestLength);
    o_.primaryAfferent = calculatePrimaryAfferent();
    o_.secondaryAfferent = calculateSecondaryAfferent();
}

void Mileusnic2006IntrafusalFiberActivationDynamics::integrate(DoubleT dt) {
    using namespace boost::numeric::odeint;
    
    if (properties_.isChain) { 
        sNew_.activation = (i_.fusimotorFrequency * i_.fusimotorFrequency)
                           / (i_.fusimotorFrequency * i_.fusimotorFrequency + p_.freq * p_.freq);

    } else {
        runge_kutta4<DoubleT> stepper;
        stepper.do_step(
            [&](const DoubleT &x, DoubleT &dxdt, const DoubleT /* t */) {
                dxdt = (i_.fusimotorFrequency * i_.fusimotorFrequency
                               / (i_.fusimotorFrequency * i_.fusimotorFrequency + p_.freq * p_.freq)
                           - x)
                       / p_.tau;
            },
            s_.activation,
            0.,
            sNew_.activation,
            dt);
    }
}

//All parameters are from Mileusnic 2006 table 1.
Mileusnic2006MuscleSpindle::Mileusnic2006MuscleSpindle() {
    Mileusnic2006IntrafusalFiberActivationDynamics::Parameters dynamicsBag1Parameters;
    dynamicsBag1Parameters.freq = 60.;
    dynamicsBag1Parameters.tau = 0.149;
    dynamicsBag1_.getParameters() = dynamicsBag1Parameters;

    Mileusnic2006IntrafusalFiberActivationDynamics::Parameters dynamicsBag2Parameters;
    dynamicsBag1Parameters.freq = 60.;
    dynamicsBag1Parameters.tau = 0.205;
    dynamicsBag2_.getParameters() = dynamicsBag2Parameters;

    Mileusnic2006IntrafusalFiberActivationDynamics::Parameters dynamicsChainParameters;
    dynamicsBag1Parameters.freq = 90.;
    dynamicsBag1Parameters.tau = 1.;
    dynamicsChain_.getParameters() = dynamicsChainParameters;
    dynamicsChain_.getProperties().isChain = true;

    Mileusnic2006IntrafusalFiber::Parameters bag1Parameters;
    bag1Parameters.sensoryRegionStiffness = 10.4649;
    bag1Parameters.polarRegionStiffness = 0.1500;
    bag1Parameters.beta0 = 0.0605;
    bag1Parameters.beta1 = 0.2592;
    bag1Parameters.beta2 = 0.;
    bag1Parameters.gamma1 = 0.0289;
    bag1Parameters.gamma2 = 0.;
    bag1Parameters.coefficientOfAsymmetryLengthening = 1.;
    bag1Parameters.coefficientOfAsymmetryShortening = 0.4200;
    bag1Parameters.percentageSecondaryAfferentOnSensoryRegion = 0.;
    bag1Parameters.sensoryRegionThresholdLength = 0.0423;
    bag1Parameters.polarRegionThresholdLength = 0.;
    bag1Parameters.sensoryRegionStretchToPrimaryAfferentFiring = 20000;
    bag1Parameters.sensoryRegionStretchToSecondaryAfferentFiring = 0;
    bag1Parameters.velocityPowerTerm = 0.3;
    bag1Parameters.muscleFascicleSlackLength = 0.46;
    bag1Parameters.sensoryRegionRestLength = 0.04;
    bag1Parameters.polarRegionRestLength = 0.76;
    bag1Parameters.secondaryAfferentRestLength = 0.;
    bag1_.getParameters() = bag1Parameters;

    Mileusnic2006IntrafusalFiber::Parameters bag2Parameters;
    bag2Parameters.sensoryRegionStiffness = 10.4649;
    bag2Parameters.polarRegionStiffness = 0.1500;
    bag2Parameters.beta0 = 0.0822;
    bag2Parameters.beta1 = 0.;
    bag2Parameters.beta2 = -0.0460;
    bag2Parameters.gamma1 = 0.;
    bag2Parameters.gamma2 = 0.0636;
    bag2Parameters.coefficientOfAsymmetryLengthening = 1.;
    bag2Parameters.coefficientOfAsymmetryShortening = 0.4200;
    bag2Parameters.percentageSecondaryAfferentOnSensoryRegion = 0.7;
    bag2Parameters.sensoryRegionThresholdLength = 0.0423;
    bag2Parameters.polarRegionThresholdLength = 0.89;
    bag2Parameters.sensoryRegionStretchToPrimaryAfferentFiring = 10000;
    bag2Parameters.sensoryRegionStretchToSecondaryAfferentFiring = 7250;
    bag2Parameters.velocityPowerTerm = 0.3;
    bag2Parameters.muscleFascicleSlackLength = 0.46;
    bag2Parameters.sensoryRegionRestLength = 0.04;
    bag2Parameters.polarRegionRestLength = 0.76;
    bag2Parameters.secondaryAfferentRestLength = 0.04;
    bag2_.getParameters() = bag2Parameters;
        
    Mileusnic2006IntrafusalFiber::Parameters chainParameters;
    chainParameters.sensoryRegionStiffness = 10.4649;
    chainParameters.polarRegionStiffness = 0.1500;
    chainParameters.beta0 = 0.0822;
    chainParameters.beta1 = 0.;
    chainParameters.beta2 = -0.0690;
    chainParameters.gamma1 = 0.;
    chainParameters.gamma2 = 0.0954;
    chainParameters.coefficientOfAsymmetryLengthening = 1.;
    chainParameters.coefficientOfAsymmetryShortening = 0.4200;
    chainParameters.percentageSecondaryAfferentOnSensoryRegion = 0.7;
    chainParameters.sensoryRegionThresholdLength = 0.0423;
    chainParameters.polarRegionThresholdLength = 0.89;
    chainParameters.sensoryRegionStretchToPrimaryAfferentFiring = 10000;
    chainParameters.sensoryRegionStretchToSecondaryAfferentFiring = 7250;
    chainParameters.velocityPowerTerm = 0.3;
    chainParameters.muscleFascicleSlackLength = 0.46;
    chainParameters.sensoryRegionRestLength = 0.04;
    chainParameters.polarRegionRestLength = 0.76;
    chainParameters.secondaryAfferentRestLength = 0.04;
    chain_.getParameters() = chainParameters;
}

void Mileusnic2006MuscleSpindle::setNormalizedMuscleFiberLength(
    DoubleT normalizedMuscleFiberLength) {
    i_.normalizedMuscleFiberLength = normalizedMuscleFiberLength;
    bag1_.setNormalizedMuscleFiberLength(normalizedMuscleFiberLength);
    bag2_.setNormalizedMuscleFiberLength(normalizedMuscleFiberLength);
    chain_.setNormalizedMuscleFiberLength(normalizedMuscleFiberLength);
}


void Mileusnic2006MuscleSpindle::setFusimotorFrequency(DoubleT gStatic, DoubleT gDynamic) {
    i_.gStatic = gStatic;
    i_.gDynamic = gDynamic;
    dynamicsBag1_.setFusimotorFrequency(gDynamic);
    dynamicsBag2_.setFusimotorFrequency(gStatic);
    dynamicsChain_.setFusimotorFrequency(gStatic);

}

void Mileusnic2006MuscleSpindle::setInput(Frequency gStatic, Frequency gDynamic) {
    setFusimotorFrequency(gStatic.get(), gDynamic.get());
}

void Mileusnic2006MuscleSpindle::setInput(NormalizedFiberLength normalizedMuscleFiberLength) {
    setNormalizedMuscleFiberLength(normalizedMuscleFiberLength.get());
}

void Mileusnic2006MuscleSpindle::setInput(Input input) {
    i_ = input;
    setNormalizedMuscleFiberLength(i_.normalizedMuscleFiberLength);
    setFusimotorFrequency(i_.gStatic, i_.gDynamic);
 }

void Mileusnic2006MuscleSpindle::integrate(DoubleT dt) {
     dynamicsBag1_.integrate(dt);
     dynamicsBag2_.integrate(dt);
     dynamicsChain_.integrate(dt);

     dynamicsBag1_.calculateOutput();
     dynamicsBag2_.calculateOutput();
     dynamicsChain_.calculateOutput();

     bag1_.setActivation(0., dynamicsBag1_.getOutput().activation);
     bag2_.setActivation(dynamicsBag2_.getOutput().activation, 0.);
     chain_.setActivation(dynamicsChain_.getOutput().activation, 0.);

     bag1_.integrate(dt);
     bag2_.integrate(dt);
     chain_.integrate(dt);
 }

void Mileusnic2006MuscleSpindle::validateState() {
     dynamicsBag1_.validateState();
    dynamicsBag2_.validateState();
     dynamicsChain_.validateState();
    bag1_.validateState();
    bag2_.validateState();
    chain_.validateState();
}

void Mileusnic2006MuscleSpindle::equilibrate() {
    bag1_.equilibrate();
    bag2_.equilibrate();
    chain_.equilibrate();
}

void Mileusnic2006MuscleSpindle::evaluate(DoubleT dt) {
    integrate(dt);
    calculateOutput();
    validateState();
}

void Mileusnic2006MuscleSpindle::calculateOutput() {
    dynamicsBag1_.calculateOutput();
    dynamicsBag2_.calculateOutput();
    dynamicsChain_.calculateOutput();
    bag1_.calculateOutput();
    bag2_.calculateOutput();
    chain_.calculateOutput();
    const auto sum = bag2_.o_.primaryAfferent + chain_.o_.primaryAfferent;
    if (bag1_.o_.primaryAfferent > sum)
        o_.primaryAfferent = bag1_.o_.primaryAfferent + p_.primaryAfferentPartialOcclusion * sum;
    else
        o_.primaryAfferent = p_.primaryAfferentPartialOcclusion * bag1_.o_.primaryAfferent + sum;
    o_.secondaryAfferent = bag2_.o_.secondaryAfferent + chain_.o_.secondaryAfferent;
}

}// namespace ceinms