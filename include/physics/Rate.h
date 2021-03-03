#pragma once

#include <hydra/Function.h>
#include <hydra/FunctorArithmetic.h>
#include <hydra/Lambda.h>
#include <hydra/Pdf.h>
#include <hydra/Plain.h>

#include <physics/Mixing.h>

#include <hydra/LogLikelihoodFCN.h>
#include <hydra/Vector4R.h>
#include <hydra/multivector.h>
#include <hydra/Integrator.h>

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
auto imag_part(Amplitude const& amplitude)
{
	auto _imag = hydra::wrap_lambda(
			  [=] __hydra_dual__ (hydra::complex<double> amp){
						 return amp.imag();
			  }
	);
	
	return hydra::compose(_imag, amplitude);
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
auto multiplyBy(Amplitude const& amplitude, const double factor)
{
	auto _multiply = hydra::wrap_lambda(
			  [=] __hydra_dual__ (hydra::complex<double> amp){
						 return amp * factor;
			  }
	);
	
	return hydra::compose(_multiply, amplitude);
}


template<typename BACKEND>
struct ConstantIntegrator;

template<hydra::detail::Backend BACKEND>
class ConstantIntegrator<hydra::detail::BackendPolicy<BACKEND>>: 
public hydra::Integral<ConstantIntegrator<hydra::detail::BackendPolicy<BACKEND>>>
{
	typedef hydra::detail::BackendPolicy<BACKEND> system_t;
public:
	ConstantIntegrator()=delete;

	ConstantIntegrator(double norm):fNorm(norm){};

	// double is not allow to be template parameter, therefore we need to treat the norm as member variables
	ConstantIntegrator( ConstantIntegrator<hydra::detail::BackendPolicy<BACKEND>> const& other):fNorm(other.fNorm){}

	ConstantIntegrator<hydra::detail::BackendPolicy<BACKEND>>&
	operator=( ConstantIntegrator<hydra::detail::BackendPolicy<BACKEND>> const& other){fNorm = other.fNorm;}

	template<hydra::detail::Backend BACKEND2>
	ConstantIntegrator( ConstantIntegrator<hydra::detail::BackendPolicy<BACKEND2>> const& other):fNorm(other.fNorm){}

	template<hydra::detail::Backend BACKEND2>
	ConstantIntegrator<hydra::detail::BackendPolicy<BACKEND>>&
    operator=( ConstantIntegrator<hydra::detail::BackendPolicy<BACKEND2>> const& other){fNorm = other.fNorm;}

    template<typename FUNCTOR> 
    inline std::pair<double, double> Integrate(FUNCTOR const& fFunctor) 
    {
    	return std::make_pair(fNorm, 0.0); // constant has no error
    }
private:
	const double fNorm;
};


// This function would build a pdf functor (the funtor is a pdf iself!) for time dependent amplitude with fixed Adir and Abar, and normalized pdf(sigma_t).
// Currently, q/p is not considered.
template<Flavor Tag, typename MSqPlus, typename MSqMinus, typename Time, typename TimeError, typename EFFICIENCY, typename ADIR, typename ABAR, typename PDFSIGMAT>
auto time_dependent_rate_with_time_resolution_pdf(hydra::Parameter const& tau, hydra::Parameter const& x, hydra::Parameter const& y, hydra::Parameter const& qop, hydra::Parameter const& phi, hydra::Parameter const& b, hydra::Parameter const& s, EFFICIENCY const& efficiency, ADIR const& Adir, ABAR const& Abar, PDFSIGMAT const& pdf_sigma_t, std::array<double,2> dalitzRangePlus, std::array<double,2> dalitzRangeMinus, std::array<double,2> timeRange, std::array<double,2> timeErrorRange, const long n_calls=50000000)
{

	// check if Adir and Abar are fixed
	std::cout << "Checking if Adir and Abar are fixed ... ..." << std::endl;

	// make a dummy fcn to easily get the parameter list
	hydra::multivector<hydra::tuple<MSqPlus,MSqMinus>, hydra::device::sys_t> data_dummy;
	data_dummy.push_back(hydra::make_tuple(MSqPlus(1.0),MSqMinus(1.0)));

	auto AdirDummyFCN = hydra::make_loglikehood_fcn( hydra::make_pdf( rate(Adir), ConstantIntegrator<hydra::device::sys_t>(1.0)) ,
													 data_dummy.begin(),
													 data_dummy.end());
	auto AbarDummyFCN = hydra::make_loglikehood_fcn( hydra::make_pdf( rate(Abar), ConstantIntegrator<hydra::device::sys_t>(1.0)) ,
													 data_dummy.begin(),
													 data_dummy.end());

    auto AdirParameters = AdirDummyFCN.GetParameters().GetVariables();
    auto AbarParameters = AbarDummyFCN.GetParameters().GetVariables();

    auto check_free_parameters = [](std::vector<hydra::Parameter* >& parameters, const char* amplitudeName) {
    	int n_free = 0;
	    std::string free_parameter_str("");
	    for (size_t i = 0; i < parameters.size(); ++i) {
	        if (!(parameters[i]->IsFixed())) {
	        	n_free++;
	        	free_parameter_str += " " + std::string(parameters[i]->GetName());
	        }
	    }

	    if (n_free != 0) {
	    	std::cout << "time_dependent_rate_with_time_resolution_pdf_type1: ";
	    	std::cout << amplitudeName << " should be fixed! But now it has " << n_free << " floating parameters: " << std::endl;
	    	std::cout << free_parameter_str << std::endl;
	        exit(-1);
    	}
    };

    check_free_parameters(AdirParameters, "Adir");
    check_free_parameters(AbarParameters, "Abar");

	std::cout << "Pass." << std::endl;

	// check if pdf(sigma_t) is normalized
	std::cout << "Checking if pdf(sigma_t) is normalized ... ..." << std::endl;

	hydra::Plain<1,  hydra::device::sys_t > PlainSigmatMC({timeErrorRange[0]}, {timeErrorRange[1]}, 5000000);
	double int_pdf_sigma_t = PlainSigmatMC.Integrate(pdf_sigma_t).first;
	std::cout << "Integrate(pdf(sigma_t)) on (timeErrorRange[0], timeErrorRange[1]) = " << int_pdf_sigma_t << std::endl;
	// double int_pdf_sigma_t = pdf_sigma_t.GetNorm();
	// std::cout << "pdf_sigma_t.GetNorm() = " << int_pdf_sigma_t << std::endl;
	if (abs(int_pdf_sigma_t - 1) > 0.01) {
		std::cout << "time_dependent_rate_with_time_resolution_pdf: ";
		std::cout << "The input pdf(sigma_t) should be normalized, but its integration seems not close enough to 1. ";
		std::cout << "Exit!" << std::endl; 
	}

	std::cout << "Pass." << std::endl;

	// integrate on dalitz plane and decay-time axis
	hydra::Plain<2,  hydra::device::sys_t > PlainDalitzMC({dalitzRangePlus[0], dalitzRangeMinus[0]}, 
														  {dalitzRangePlus[1], dalitzRangeMinus[1]}, 
														  n_calls);

	// numerically integrate for time-independent part
	// suppose the efficiency plane is time-independent and has not correlation with decay-time or sigmat.

	auto rcp = QoverP<Tag>(qop,phi);

	auto AlphaSq = efficiency * rate(Adir);
	auto BetaSq = efficiency * rate(rcp * Abar);
	auto AlphaBetaStar_real = efficiency * real_part(Adir * conjugate(rcp * Abar)); 
	auto AlphaBetaStar_imag = efficiency * imag_part(Adir * conjugate(rcp * Abar)); 
	
	double AlphaSq_int = PlainDalitzMC.Integrate(AlphaSq).first;
	double BetaSq_int = PlainDalitzMC.Integrate(BetaSq).first;
	double AlphaBetaStar_real_int = PlainDalitzMC.Integrate(AlphaBetaStar_real).first;
	double AlphaBetaStar_imag_int = PlainDalitzMC.Integrate(AlphaBetaStar_imag).first;

	double AsumD2Sq_int = 1. / 4. * ( AlphaSq_int + BetaSq_int + 2*AlphaBetaStar_real_int );
	double AdiffD2Sq_int = 1. / 4. * ( AlphaSq_int + BetaSq_int - 2*AlphaBetaStar_real_int );
	hydra::complex<double> AsumD2AdiffD2Star_int = 1. / 4. * hydra::complex<double>( AlphaSq_int - BetaSq_int , -2*AlphaBetaStar_imag_int );


	// build the final functor

	auto psi_p = MixingPsip<Time, TimeError>(y,tau,s,b);
	auto psi_m = MixingPsim<Time, TimeError>(y,tau,s,b);
	auto psi_i = MixingPsii<Time, TimeError>(x,tau,s,b);

	auto As = Adir + rcp * Abar; 
	auto Ad = conjugate( Adir - rcp * Abar ) ; 

	auto int_on_dalitzplane_decaytime = MixingNormType1<TimeError>(AlphaSq_int, BetaSq_int, AsumD2AdiffD2Star_int, 
																	  	 x, y, tau, s, b, timeRange);

	// the hydra::wrap_lambda and hydra::compose() below would handle the parameter list of Adir and Abar automatically
	auto _arrange = hydra::wrap_lambda(
			  [] __hydra_dual__ (hydra::tuple< double, hydra::complex<double>, hydra::complex<double>, hydra::complex<double>, hydra::complex<double>, hydra::complex<double>, double, double > input_amplitudes){
			  		auto _efficiency = hydra::get<0>(input_amplitudes);
			  		auto _As = hydra::get<1>(input_amplitudes);
			  		auto _Ad = hydra::get<2>(input_amplitudes);
			  		auto _psi_p = hydra::get<3>(input_amplitudes);
			  		auto _psi_m = hydra::get<4>(input_amplitudes);
			  		auto _psi_i = hydra::get<5>(input_amplitudes);
			  		auto _pdf_sigma_t = hydra::get<6>(input_amplitudes);
			  		auto _int_on_dalitzplane_decaytime = hydra::get<7>(input_amplitudes);

			  		return  _efficiency * ( norm(_As/2.)*_psi_p.real() + norm(_Ad/2.)*_psi_m.real() + (_As*_Ad/2.*_psi_i).real() ) * _pdf_sigma_t / _int_on_dalitzplane_decaytime;

			  }
	);
	
	return hydra::compose(_arrange, efficiency, As, Ad, psi_p, psi_m, psi_i, pdf_sigma_t, int_on_dalitzplane_decaytime);


	// auto _arrange = hydra::wrap_lambda(
	// 	  [=] __hydra_dual__ (hydra::tuple< hydra::complex<double>, hydra::complex<double>, hydra::complex<double>, hydra::complex<double>, hydra::complex<double>, double > input_amplitudes){
	// 	  		auto _As = hydra::get<0>(input_amplitudes);
	// 	  		auto _Ad = hydra::get<1>(input_amplitudes);
	// 	  		auto _psi_p = hydra::get<2>(input_amplitudes);
	// 	  		auto _psi_m = hydra::get<3>(input_amplitudes);
	// 	  		auto _psi_i = hydra::get<4>(input_amplitudes);
	// 	  		auto _pdf_sigma_t = hydra::get<5>(input_amplitudes);

	// 	  		return  ( norm(_As/2.)*_psi_p.real() + norm(_Ad/2.)*_psi_m.real() + (_As*_Ad/2.*_psi_i).real() ) * _pdf_sigma_t ;

	// 	  }
	// );
	
	// return hydra::compose(_arrange, As, Ad, psi_p, psi_m, psi_i, pdf_sigma_t);

}


// This function would build a pdf functor for time dependent amplitude with unfixed Adir and Abar, and normalized pdf(sigma_t).
// template<Flavor Tag, typename MSqPlus, typename MSqMinus, typename Time, typename TimeError, typename ADIR, typename ABAR, typename PDFSIGMAT>
// dummy time_dependent_rate_with_time_resolution_pdf_type2; // to be developed. Probably we need a new functor with _integrated_Adir 
// and _integrated_Abar as member variables, and update this member variables when it finds the fAdir/fAbar parameter list is changed, 
// while using this _integrated_Adir and _integrated_Abar in the evaluate() when calculating normalization factor. This class/function 
// would be developed in the future if needed






} // namespace dafne
