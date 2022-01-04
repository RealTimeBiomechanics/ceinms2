#include "ceinms2/Mileusnic2006MuscleSpindle.h"
#include "ceinms2/WDBsolver.h"
#include <cmath>
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

        const auto fvpr = calculateVelocityContributionToPolarRegionForce(polarRegionLength,
            polarRegionVelocity,
            beta,
            this->p_.coefficientOfAsymmetry,
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
    return p_.sensoryRegionStretchToAfferentFiring
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
    const auto term2 = p_.sensoryRegionStretchToAfferentFiring
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

void Mileusnic2006MuscleSpindle::setActivation(DoubleT fStatic, DoubleT fDynamic) {
    i_.fStatic = fStatic;
    i_.fDynamic = fDynamic;

    bag1_.setActivation(0., fDynamic);
    bag2_.setActivation(fStatic, 0.);
    chain_.setActivation(fStatic, 0.);
}

void Mileusnic2006MuscleSpindle::setNormalizedMuscleFiberLength(
    DoubleT normalizedMuscleFiberLength) {
    i_.normalizedMuscleFiberLength = normalizedMuscleFiberLength;
    bag1_.setNormalizedMuscleFiberLength(normalizedMuscleFiberLength);
    bag2_.setNormalizedMuscleFiberLength(normalizedMuscleFiberLength);
    chain_.setNormalizedMuscleFiberLength(normalizedMuscleFiberLength);
}

void Mileusnic2006MuscleSpindle::setInput(Activation fStatic, Activation fDynamic) {
    i_.fStatic = fStatic.get();
    i_.fDynamic = fDynamic.get();

    bag1_.setActivation(0., i_.fDynamic);
    bag2_.setActivation(i_.fStatic, 0.);
    chain_.setActivation(i_.fStatic, 0.);
}
void Mileusnic2006MuscleSpindle::setInput(Input input) {
    i_ = input;
    setNormalizedMuscleFiberLength(i_.normalizedMuscleFiberLength);
    setActivation(i_.fStatic, i_.fDynamic);
}

void Mileusnic2006MuscleSpindle::integrate(DoubleT dt) {
    bag1_.integrate(dt);
    bag2_.integrate(dt);
    chain_.integrate(dt);
}

void Mileusnic2006MuscleSpindle::validateState() {
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