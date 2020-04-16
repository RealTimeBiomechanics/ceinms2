#ifndef ceinms_Lloyd2019Muscle_h
#define ceinms_Lloyd2019Muscle_h

#include "ceinms/Curve.h"
#include <ensmallen.hpp>
#include <iostream>
#include "ceinms/WDBsolver.h"
#include <autodiff/forward.hpp>
using autodiff::forward::sin;
using autodiff::forward::asin;
using autodiff::forward::cos;
using autodiff::forward::pow;

//using DoubleT = autodiff::dual;

using DoubleT = double;
using CurveOffline = ceinms::Curve<ceinms::CurveMode::Offline, DoubleT, ceinms::CurveMode::Interpolation::Cubic> ;
namespace ceinms {
	class ForceErrorFunction;

	class Lloyd2019Muscle {

	public:
		struct Parameters {
			DoubleT optimalFibreLength;
			DoubleT pennationAngleAtOptimalFibreLength;
			DoubleT tendonSlackLength;
			DoubleT maxContractionVelocity;
			DoubleT maxIsometricForce;
			DoubleT strengthCoefficient;
			DoubleT damping;
			DoubleT percentageChange;
			CurveOffline  forceVelocityCurve;
			CurveOffline  activeForceLengthCurve;
			CurveOffline  passiveForceLengthCurve;
			CurveOffline  tendonForceStrainCurve;
		};

		struct Input {
			Input() :
				activation(0),
				musculotendonLength(0) {}

			DoubleT activation;
			DoubleT musculotendonLength;
		};

		struct State {	
			State() :
				fibreLength(0),
				fibreVelocity(0) {}
				
			DoubleT fibreLength;
			DoubleT fibreVelocity;

		};
		struct Output {
			Output() :
				optimalFibreLengthAtT(0),
				pennationAngle(0),
				fibreForce(0),
				activeForce(0),
				passiveForce(0),
				dampingForce(0),
				normalisedFibreLength(0),
				normalisedFibreLengthAtT(0),
				normalisedFibreVelocity(0),
				tendonLength(0),
				tendonStrain(0),
				fa(0),
				fv(0),
				fp(0) {}

			DoubleT optimalFibreLengthAtT;
			DoubleT pennationAngle;
			DoubleT fibreForce;
			DoubleT activeForce;
			DoubleT passiveForce;
			DoubleT dampingForce;
			DoubleT normalisedFibreLength;
			DoubleT normalisedFibreLengthAtT;
			DoubleT normalisedFibreVelocity;
			DoubleT tendonLength;
			DoubleT tendonStrain;
			DoubleT tendonForce;
			DoubleT fa;
			DoubleT fv;
			DoubleT fp;
		};
	
		Lloyd2019Muscle(Parameters parameters) :
			p_(parameters) {
			s_.fibreLength = p_.optimalFibreLength;
			s_.fibreVelocity = 0.;
		}
		/* calculateXYZ functions perform calculations and return the calculated value*/
		void setInput(DoubleT activation, DoubleT musculotendonLength);
		void equilibrate();
		void setState(State state);
		//From input and current state calculate the new state of the system
		void integrate(DoubleT dt);
		//The temporary state calculated via `integrate` becomes the new state
		void validateState();
		//from the internal state of the system and the iput, calculate all the output;
		void calculateOutput();

		Parameters& updParameters() { return p_; }
		State& updState() { return s_; }
		const State& getState() { return s_; }
		Output& updOutput() { return o_; }

	
		static DoubleT calculateFibreVelocityFromFibreLength(DoubleT previousFibreLength, DoubleT currentFibreLength, DoubleT dt);
		static DoubleT calculateFibreForce(DoubleT activation, DoubleT fibreLength, DoubleT fibreVelocity, const Parameters& p);
		static DoubleT calculateTendonForce(DoubleT musculotendonLength, DoubleT fibreLength, const Parameters& p);
		static DoubleT calculatePennationAngle(DoubleT fibreLength, const Parameters& p);
		static DoubleT calculateOptimalFibreLengthAtT(DoubleT activation, const Parameters& p);
		static DoubleT calculateNormalisedFibreLengthAtT(DoubleT activation, DoubleT fibreLength, const Parameters& p);
		static DoubleT calculateNormalisedFibreLength(DoubleT fibreLength, const Parameters& p);
		static DoubleT calculateNormalisedFibreVelocity(DoubleT fibreVelocity, const Parameters& p);
		static DoubleT calculateTendonLength(DoubleT musculotendonLength, DoubleT fibreLength, const Parameters& p);
		static DoubleT calculateTendonStrain(DoubleT musculotendonLength, DoubleT fibreLength, const Parameters& p);

	private:
		DoubleT integrateFibreLength(DoubleT dt);
		Parameters p_;
		State s_, sNew_;
		Input i_;
		Output o_;

	};
	/*
	class ForceErrorFunction
	{
	public:
		// Construct the object with the given data matrix and responses.
		ForceErrorFunction(Lloyd2019Muscle::Input input, Lloyd2019Muscle::State s, Lloyd2019Muscle::Parameters p, DoubleT dt) :
			activation_(input.activation), musculotendonLength_(input.musculotendonLength), dt_(dt), s_(s), p_(p) { }

		// Return the objective function for model parameters x.
		double Evaluate(const arma::mat& x)
		{
			DoubleT fibreLength = x.at(0);
			DoubleT fibreLengthT1 = s_.fibreLength;
			DoubleT fibreVelocity = Lloyd2019Muscle::calculateFibreVelocityFromFibreLength(fibreLengthT1, fibreLength, dt_);
			double v = f(activation_, musculotendonLength_, fibreLength, fibreVelocity, p_).val;
			std::cout << "f(" << x << ") = " << v <<  std::endl;
			return v;
		}

		// Compute the gradient for model parameters x.
		void Gradient(const arma::mat& x, arma::mat& g)
		{
			DoubleT fibreLength = x.at(0);
			DoubleT fibreLengthT1{ s_.fibreLength };
			DoubleT fibreVelocity = (fibreLength - fibreLengthT1) / (2 * dt_);
			double dfdl = derivative(f, wrt(fibreLength), at(activation_, musculotendonLength_, fibreLength, fibreVelocity, p_));

			g.at(0) = dfdl;
			std::cout << "g(" << x << ") = " << g << std::endl;
		}

		/*	size_t NumConstraints() { return 4; }


			double EvaluateConstraint(const size_t i, const arma::mat& x) {
				const arma::mat lambda{0.2,1.8,-1,1};
				double normFibreLenght = x.at(0) / muscle_.updParameters().optimalFibreLength.val;
				double normFibreVel = x.at(1) / muscle_.updParameters().maxContractionVelocity.val;
				double val = DBL_MAX;
				if (i == 0) {
					if (x.at(0) > lambda.at(i))
						val = 0;
				}
				if (i == 1) {
					if (x.at(0) < lambda.at(i))
						val = 0;
				}
				if (i == 2) {
					if (x.at(1) > lambda.at(i))
						val = 0;
				}
				if (i == 3) {
					if (x.at(1) < lambda.at(i))
						val = 0;
				}
		//		std::cout << "constraint " << val << std::endl;
				return val;
			}

			void GradientConstraint(const size_t i, const arma::mat& x, arma::mat& g) {
				g.copy_size(x);
				g.at(0) = 0.;
				g.at(1) = 0.;
			}
			double EvaluateWithGradient(const arma::mat& x, arma::mat& g) {

				DoubleT fibreLength = x.at(0);
				DoubleT fibreVelocity = x.at(1);

				double dfdl = derivative(f, wrt(fibreLength), at(muscle_, activation_, musculotendonLength_, fibreLength, fibreVelocity));
				double dfdv = derivative(f, wrt(fibreVelocity), at(muscle_, activation_, musculotendonLength_, fibreLength, fibreVelocity));

				g.at(0) = dfdl;
				g.at(1) = dfdv;
				return f(muscle_, activation_, musculotendonLength_, fibreLength, fibreVelocity);
			}
			
	private:
		static DoubleT f(DoubleT activation, DoubleT musculotendonlength, DoubleT fibreLength, DoubleT fibreVelocity, Lloyd2019Muscle::Parameters& p) {
			DoubleT tendonForce = Lloyd2019Muscle::calculateTendonForce(musculotendonlength, fibreLength, p);
			DoubleT fibreForce = Lloyd2019Muscle::calculateFibreForce(activation, fibreLength, fibreVelocity, p);
			return pow(tendonForce - fibreForce, 2.0);
		}

		Lloyd2019Muscle::Parameters p_;
		Lloyd2019Muscle::State s_;
		DoubleT musculotendonLength_;
		DoubleT activation_;
		DoubleT dt_;
	};
	*/
	void Lloyd2019Muscle::setInput(DoubleT activation, DoubleT musculotendonLength) {
		i_.activation = activation;
		i_.musculotendonLength = musculotendonLength;
	}

	void Lloyd2019Muscle::equilibrate() {

		double diff = 1;
		double fLength = s_.fibreLength;
		double tol = 1e-9;
		while (diff > tol) {
			integrate(0.001);
			diff = fabs(sNew_.fibreLength - fLength);
			fLength = sNew_.fibreLength;
		}
		validateState();
	}

	void Lloyd2019Muscle::integrate(DoubleT dt) {
		sNew_.fibreLength = integrateFibreLength(dt);
		sNew_.fibreVelocity = calculateFibreVelocityFromFibreLength(s_.fibreLength, sNew_.fibreLength, dt);
	}

	void Lloyd2019Muscle::validateState() {
		s_ = sNew_;
	}

	void Lloyd2019Muscle::calculateOutput() {

		o_.normalisedFibreLengthAtT =  calculateNormalisedFibreLengthAtT(i_.activation, s_.fibreLength, p_);
		o_.normalisedFibreLength = calculateNormalisedFibreLength(s_.fibreLength, p_);
		o_.normalisedFibreVelocity = calculateNormalisedFibreVelocity(s_.fibreVelocity, p_);
		o_.pennationAngle = calculatePennationAngle(s_.fibreLength, p_);
		o_.fa = p_.activeForceLengthCurve.getValue(o_.normalisedFibreLengthAtT);
		o_.fp = p_.passiveForceLengthCurve.getValue(o_.normalisedFibreLength);
		o_.fv = p_.forceVelocityCurve.getValue(o_.normalisedFibreVelocity);
		o_.activeForce = p_.maxIsometricForce*p_.strengthCoefficient*o_.fa*o_.fv;
		o_.passiveForce = p_.maxIsometricForce*p_.strengthCoefficient*o_.fp;
		o_.dampingForce = p_.maxIsometricForce*p_.strengthCoefficient * p_.damping * o_.normalisedFibreVelocity;
		o_.fibreForce = calculateFibreForce(i_.activation, s_.fibreLength, s_.fibreVelocity, p_);
		o_.tendonStrain = calculateTendonStrain(i_.musculotendonLength, s_.fibreLength, p_);
		o_.tendonLength = calculateTendonLength(i_.musculotendonLength, s_.fibreLength, p_);
		o_.tendonForce = calculateTendonForce(i_.musculotendonLength, s_.fibreLength, p_);
	}

	DoubleT Lloyd2019Muscle::calculateFibreVelocityFromFibreLength(DoubleT previousFibreLength, DoubleT currentFibreLength, DoubleT dt) {

		return (currentFibreLength - previousFibreLength) / dt;
	}

	DoubleT Lloyd2019Muscle::integrateFibreLength(DoubleT dt) {

		DoubleT optimalFibreLengthAtT = calculateOptimalFibreLengthAtT(i_.activation, p_);

		DoubleT minFibreLength = 0.2*p_.optimalFibreLength;
		DoubleT maxFibreLength = 2 * p_.optimalFibreLength;
		double activation = i_.activation;
		double musculotendonLength = i_.musculotendonLength;

		auto f = [&](double fibreLength) -> double {
			DoubleT fibreVelocity = Lloyd2019Muscle::calculateFibreVelocityFromFibreLength(this->s_.fibreLength, fibreLength, dt);
			DoubleT tendonForce = Lloyd2019Muscle::calculateTendonForce(this->i_.musculotendonLength, fibreLength, this->p_);
			DoubleT fibreForce = Lloyd2019Muscle::calculateFibreForce(this->i_.activation, fibreLength, fibreVelocity, this->p_);
			DoubleT v = tendonForce - fibreForce;
			return v;
			};

//		ceinms::ForceErrorFunction fun(i_, s_, p_, 0.001);
		double currentFibreLength = wdbSolve(f, minFibreLength, maxFibreLength, 1e-9);
	//	ens::LineSearch optimizer;
//		arma::mat x1{minFibreLength.val }, x2{ maxFibreLength.val };
//		optimizer.Optimize(fun, x1, x2);
//		return x2.at(0);
		return currentFibreLength;
	}

	DoubleT Lloyd2019Muscle::calculateOptimalFibreLengthAtT(DoubleT activation, const Parameters& p) {

		return p.optimalFibreLength * (p.percentageChange *(1.0 - activation) + 1);
	}

	DoubleT Lloyd2019Muscle::calculateNormalisedFibreLengthAtT(DoubleT activation, DoubleT fibreLength, const Parameters& p) {
		
		return fibreLength / calculateOptimalFibreLengthAtT(activation, p);
	}

	DoubleT Lloyd2019Muscle::calculateNormalisedFibreLength(DoubleT fibreLength, const Parameters& p) {

		return fibreLength / p.optimalFibreLength;
	}

	DoubleT Lloyd2019Muscle::calculateNormalisedFibreVelocity(DoubleT fibreVelocity, const Parameters& p) {

		DoubleT normalisedFibreVelocity = fibreVelocity / (p.optimalFibreLength * p.maxContractionVelocity);
		if (normalisedFibreVelocity > 1) normalisedFibreVelocity = 1;
		if (normalisedFibreVelocity < -1) normalisedFibreVelocity = -1;
		return normalisedFibreVelocity;
	}

	DoubleT Lloyd2019Muscle::calculatePennationAngle(DoubleT fibreLength, const Parameters& p) {

		DoubleT value{ p.optimalFibreLength*sin(p.pennationAngleAtOptimalFibreLength) / fibreLength };
		if (value < 0) value = 0;
		if (value > 0.99) value = 0.99;
		return asin(value);
	}

	DoubleT Lloyd2019Muscle::calculateTendonLength(DoubleT musculotendonLength, DoubleT fibreLength, const Parameters& p) {
	
		return musculotendonLength - fibreLength * cos(calculatePennationAngle(fibreLength, p));
	}
	
	DoubleT Lloyd2019Muscle::calculateTendonStrain(DoubleT musculotendonLength, DoubleT fibreLength, const Parameters& p) {
		DoubleT tendonStrain = (calculateTendonLength(musculotendonLength, fibreLength, p) - p.tendonSlackLength) / p.tendonSlackLength;
		if(tendonStrain < 0.) tendonStrain = 0;
		return tendonStrain;
		
	}

	DoubleT Lloyd2019Muscle::calculateFibreForce(DoubleT activation, DoubleT fibreLength, DoubleT fibreVelocity, const Parameters& p) {

		DoubleT normalisedFibreLengthAtT{ calculateNormalisedFibreLengthAtT(activation, fibreLength, p) };
		DoubleT normalisedFibreLength{ calculateNormalisedFibreLength(fibreLength, p) };
		DoubleT normalisedFibreVelocity{ calculateNormalisedFibreVelocity(fibreVelocity, p) };
		DoubleT pennationAngle{ calculatePennationAngle(fibreLength, p) };
		DoubleT fa{ p.activeForceLengthCurve.getValue(normalisedFibreLengthAtT) };
		DoubleT fp{ p.passiveForceLengthCurve.getValue(normalisedFibreLength) };
		DoubleT fv{ p.forceVelocityCurve.getValue(normalisedFibreVelocity) };

		DoubleT fibreForce{ p.maxIsometricForce * p.strengthCoefficient*
			(fa*fv*activation + fp + p.damping * normalisedFibreVelocity)*
			cos(pennationAngle) };
		//clamp muscle force
		if (fibreForce < 0.) fibreForce = 0;
		return fibreForce;
	}


	DoubleT Lloyd2019Muscle::calculateTendonForce(DoubleT musculotendonLength, DoubleT fibreLength, const Parameters& p) {
			
		DoubleT tendonStrain{ calculateTendonStrain(musculotendonLength, fibreLength, p) };
		DoubleT tendonForce{ p.strengthCoefficient * p.maxIsometricForce*
			p.tendonForceStrainCurve.getValue(tendonStrain) };
		return tendonForce;
	}	

	
	
	
}

#endif