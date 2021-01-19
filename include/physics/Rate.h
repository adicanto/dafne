#pragma once

#include <hydra/Function.h>
#include <hydra/FunctorArithmetic.h>
#include <hydra/Lambda.h>

#include <physics/Mixing.h>

#include "RooMathM.h"

namespace dafne {

template<typename Amplitude>
auto rate(Amplitude const& amplitude)
{
	auto norm = hydra::wrap_lambda(
			  [=] __hydra_dual__ (hydra::complex<double> amp){
						 return hydra::norm( amp );
			  }
	);
	
	return hydra::compose(norm, amplitude);
}

template<typename ...Amplitudes>
auto rate(Amplitudes const& ... amplitudes)
{
	auto tot = hydra::sum(amplitudes...);
	
	auto norm = hydra::wrap_lambda(
			  [=] __hydra_dual__ (hydra::complex<double> amp){
						 return hydra::norm( amp );
			  }
	);
	
	return hydra::compose(norm, tot);
}

template<Flavor Tag, typename Time, typename ADIR, typename ABAR>
auto time_dependent_rate(hydra::Parameter const& tau, hydra::Parameter const& x, hydra::Parameter const& y, hydra::Parameter const& qop, hydra::Parameter const& phi, ADIR const& Adir, ABAR const& Abar)
{
	// Mixing functors
	auto gp  = MixingGp<Time>(tau,x,y);
	auto gm  = MixingGm<Time>(tau,x,y);
	
	// CP-violation functor
	auto rcp = QoverP<Tag>(qop,phi);
	
	// Time-dependent rate
	auto amp = gp * Adir + rcp * gm * Abar;
	return rate(amp);
}

template<typename T>
T _phi(T z)
{
	//return exp(z*z/2.) * 1./2. * std::complex<double>(RooMath::erfc(z/sqrt(2.))); 

	return exp(std::norm(z)/2.) * 1./2. * RooMathM::erfc(-z/sqrt(2.));

	//return (z.real() >= 0) ? exp(std::norm(z)/2.) * 1./2. * std::complex<double>(1. + RooMath::erf(z/sqrt(2.))) :
	//                         exp(std::norm(z)/2.) * 1./2. * std::complex<double>(1. - RooMath::erf(-z/sqrt(2.))) ;  
}

template<typename T>
T _psi(double x, T kappa)
{
	// return exp(-x*x/2.) * _phi<T>(x-kappa);
	return exp((-2*x*kappa+kappa*kappa)/2.) * 1./2. * (1. + RooMathM::erf((x - kappa)/sqrt(2.)));
}


template<Flavor Tag, typename MSqPlus, typename MSqMinus, typename Time, typename TimeError, typename ADIR, typename ABAR, typename PDFSIGMAT>
auto time_dependent_rate_with_time_resolution(hydra::Parameter const& tau, hydra::Parameter const& x, hydra::Parameter const& y, hydra::Parameter const& qop, hydra::Parameter const& phi, hydra::Parameter const& b, hydra::Parameter const& s, ADIR const& Adir, ABAR const& Abar, PDFSIGMAT const& pdf_sigma_t)
{
	auto rcp = QoverP<Tag>(qop,phi);

	auto T2 = hydra::wrap_lambda(
			  [=] __hydra_dual__ (MSqPlus m2p, MSqMinus m2m, Time t, TimeError sigma_t){
			  		// m2p = 1.5;
			  		// m2m = 0.8;


			  		double _tau = double(tau); // turn hydra::Parameter to double
			  		double _x = double(x);
			  		double _y = double(y);
			  		double _b = double(b);
			  		double _s = double(s);

			  		// A sum = (A + A-bar)/2
			  		std::complex<double> As = (Adir(m2p, m2m) + rcp()*Abar(m2p, m2m))/2; 


			  		// A difference = (A* - A-bar*)/2
			  		std::complex<double> Ad = (hydra::conj(Adir(m2p, m2m)) - hydra::conj(rcp()*Abar(m2p, m2m)))/2; 

			  		double Gamma = 1./_tau;
			  		double chi = (t-_b) / (_s*sigma_t);
			  		double kappa_p = (1.-_x)*Gamma*_s*sigma_t;
			  		double kappa_m = (1.+_x)*Gamma*_s*sigma_t;
			  		std::complex<double> kappa_i(Gamma*_s*sigma_t, -_y*Gamma*_s*sigma_t);

			  		double result = norm(As)*_psi(chi, kappa_p) + norm(Ad)*_psi(chi, kappa_m) + 2.*(As*Ad*_psi(chi, std::complex<double>(kappa_i))).real();
			  		result = result * pdf_sigma_t(sigma_t);

			  		// return norm(As)*_psi(chi, kappa_p) + norm(Ad)*_psi(chi, kappa_m) + 2.*(As*Ad*_psi(chi, std::complex<double>(kappa_i))).real();
			  		// return pdf_sigma_t(sigma_t); // for debugging
			  		return result;
			  }
	);
	
	return T2;
}



} // namespace dafne
