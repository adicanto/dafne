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

template<typename Amplitude>
auto real_part(Amplitude const& amplitude)
{
	auto _real = hydra::wrap_lambda(
			  [=] __hydra_dual__ (hydra::complex<double> amp){
						 return amp.real();
			  }
	);
	
	return hydra::compose(_real, amplitude);
}

template<typename Amplitude>
auto conjugate(Amplitude const& amplitude)
{
	auto _conj = hydra::wrap_lambda(
			  [=] __hydra_dual__ (hydra::complex<double> amp){
						 return hydra::conj(amp);
			  }
	);
	
	return hydra::compose(_conj, amplitude);
}

template<typename Amplitude>
auto divideBy(Amplitude const& amplitude, const double denominator)
{
	auto _divide = hydra::wrap_lambda(
			  [=] __hydra_dual__ (hydra::complex<double> amp){
						 return amp / denominator;
			  }
	);
	
	return hydra::compose(_divide, amplitude);
}

template<typename Amplitude>
auto rateOfDivideBy(Amplitude const& amplitude, const double denominator)
{
	auto _norm_divide = hydra::wrap_lambda(
			  [=] __hydra_dual__ (hydra::complex<double> amp){
						 return hydra::norm( amp / denominator );
			  }
	);
	
	return hydra::compose(_norm_divide, amplitude);
}


template<Flavor Tag, typename MSqPlus, typename MSqMinus, typename Time, typename TimeError, typename ADIR, typename ABAR, typename PDFSIGMAT>
auto time_dependent_rate_with_time_resolution(hydra::Parameter const& tau, hydra::Parameter const& x, hydra::Parameter const& y, hydra::Parameter const& qop, hydra::Parameter const& phi, hydra::Parameter const& b, hydra::Parameter const& s, ADIR const& Adir, ABAR const& Abar, PDFSIGMAT const& pdf_sigma_t)
{
	auto rcp = QoverP<Tag>(qop,phi);

	auto psi_p = MixingPsip<Time, TimeError>(y,tau,s,b);
	auto psi_m = MixingPsim<Time, TimeError>(y,tau,s,b);
	auto psi_i = MixingPsii<Time, TimeError>(x,tau,s,b);

	auto As = Adir + rcp * Abar; 
	auto Ad = conjugate( Adir - rcp * Abar ) ; 

	auto _arrange = hydra::wrap_lambda(
			  [=] __hydra_dual__ (hydra::tuple< hydra::complex<double>, hydra::complex<double>, hydra::complex<double>, hydra::complex<double>, hydra::complex<double>, double > input_amplitudes){
			  		auto _As = hydra::get<0>(input_amplitudes);
			  		auto _Ad = hydra::get<1>(input_amplitudes);
			  		auto _psi_p = hydra::get<2>(input_amplitudes);
			  		auto _psi_m = hydra::get<3>(input_amplitudes);
			  		auto _psi_i = hydra::get<4>(input_amplitudes);
			  		auto _pdf_sigma_t = hydra::get<5>(input_amplitudes);

			  		return  ( norm(_As/2.)*_psi_p.real() + norm(_Ad/2.)*_psi_m.real() + (_As*_Ad/2.*_psi_i).real() ) * _pdf_sigma_t ;
			  }
	);
	
	return hydra::compose(_arrange, As, Ad, psi_p, psi_m, psi_i, pdf_sigma_t);

}


} // namespace dafne
