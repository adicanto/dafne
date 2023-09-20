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

#include <tools/FunctorTools.h>

#include "RooMathM.h"


namespace dafne {


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


// This function would build a pdf functor (the funtor is a pdf iself!) for time dependent amplitude with fixed Adir and Abar, and normalized pdf(sigma_t).
// Currently, q/p is not considered.
template<Flavor Tag, typename MSqPlus, typename MSqMinus, typename Time, typename TimeError, typename EFFICIENCY, typename ADIR, typename ABAR, typename PDFSIGMAT>
auto time_dependent_rate_with_time_resolution_pdf(hydra::Parameter const& tau, hydra::Parameter const& x, hydra::Parameter const& y, hydra::Parameter const& qop, hydra::Parameter const& phi, hydra::Parameter const& b, hydra::Parameter const& s, EFFICIENCY const& efficiency, ADIR const& Adir, ABAR const& Abar, PDFSIGMAT const& pdf_sigma_t, std::array<double,2> dalitzRangePlus, std::array<double,2> dalitzRangeMinus, std::array<double,2> timeRange, std::array<double,2> timeErrorRange, const long n_calls=50000000)
{
	// check if pdf(sigma_t) is normalized
	std::cout << "Checking if pdf(sigma_t) is normalized ... ..." << std::endl;

	hydra::Plain<1,  hydra::device::sys_t > sigma_t_integrator({timeErrorRange[0]}, {timeErrorRange[1]}, n_calls);
	double int_pdf_sigma_t = sigma_t_integrator.Integrate(pdf_sigma_t).first;
	std::cout << "Integrate(pdf(sigma_t)) on (timeErrorRange[0], timeErrorRange[1]) = " << int_pdf_sigma_t << std::endl;
	// double int_pdf_sigma_t = pdf_sigma_t.GetNorm();
	// std::cout << "pdf_sigma_t.GetNorm() = " << int_pdf_sigma_t << std::endl;
	if (abs(int_pdf_sigma_t - 1) > 0.01) {
		std::cout << "time_dependent_rate_with_time_resolution_pdf: ";
		std::cout << "The input pdf(sigma_t) should be normalized, but its integration seems not close enough to 1. ";
		std::cout << "Exit!" << std::endl; 
	}

	std::cout << "Pass." << std::endl;

	hydra::Plain<2,  hydra::device::sys_t > dalitz_integrator({dalitzRangePlus[0], dalitzRangeMinus[0]}, 
													  {dalitzRangePlus[1], dalitzRangeMinus[1]}, 
													  n_calls);

	auto psi_p = MixingPsip<Time, TimeError>(y,tau,s,b);
	auto psi_m = MixingPsim<Time, TimeError>(y,tau,s,b);
	auto psi_i = MixingPsii<Time, TimeError>(x,tau,s,b);
	auto rcp = QoverP<Tag>(qop,phi);
	return make_mixing_pdf_functor<MSqPlus, MSqMinus, Time, TimeError>(
								efficiency, Adir, Abar, psi_p, psi_m, psi_i, rcp,
					 			pdf_sigma_t, dalitz_integrator, timeRange
					 		); 


}


template<Flavor Tag, typename MSqPlus, typename MSqMinus, typename Time, typename TimeError, typename EFFICIENCY, typename ADIR, typename ABAR, typename PDFSIGMAT, typename FRND, typename FMISTAG, typename FCMB, typename PDFCMB>
auto time_dependent_rate_with_time_resolution_with_background_pdf(hydra::Parameter const& tau, hydra::Parameter const& x, hydra::Parameter const& y, hydra::Parameter const& qop, hydra::Parameter const& phi, hydra::Parameter const& b, hydra::Parameter const& s, EFFICIENCY const& efficiency, ADIR const& Adir, ABAR const& Abar, PDFSIGMAT const& pdf_sigma_t, std::array<double,2> dalitzRangePlus, std::array<double,2> dalitzRangeMinus, std::array<double,2> timeRange, std::array<double,2> timeErrorRange, const long n_calls, FRND const& f_rnd, FMISTAG const& f_mistag, FCMB const& f_cmb, PDFCMB const& pdf_cmb)
{
	// check if pdf(sigma_t) is normalized
	std::cout << "Checking if pdf(sigma_t) is normalized ... ..." << std::endl;

	hydra::Plain<1,  hydra::device::sys_t > sigma_t_integrator({timeErrorRange[0]}, {timeErrorRange[1]}, n_calls);
	double int_pdf_sigma_t = sigma_t_integrator.Integrate(pdf_sigma_t).first;
	std::cout << "Integrate(pdf(sigma_t)) on (timeErrorRange[0], timeErrorRange[1]) = " << int_pdf_sigma_t << std::endl;
	// double int_pdf_sigma_t = pdf_sigma_t.GetNorm();
	// std::cout << "pdf_sigma_t.GetNorm() = " << int_pdf_sigma_t << std::endl;
	if (abs(int_pdf_sigma_t - 1) > 0.01) {
		std::cout << "time_dependent_rate_with_time_resolution_pdf: ";
		std::cout << "The input pdf(sigma_t) should be normalized, but its integration seems not close enough to 1. ";
		std::cout << "Exit!" << std::endl; 
	}

	std::cout << "Pass." << std::endl;


	// integrator on dalitz plane
	hydra::Plain<2,  hydra::device::sys_t > dalitz_integrator({dalitzRangePlus[0], dalitzRangeMinus[0]}, 
														  {dalitzRangePlus[1], dalitzRangeMinus[1]}, 
														  n_calls);


	auto psi_p = MixingPsip<Time, TimeError>(y,tau,s,b);
	auto psi_m = MixingPsim<Time, TimeError>(y,tau,s,b);
	auto psi_i = MixingPsii<Time, TimeError>(x,tau,s,b);
	auto rcp = QoverP<Tag>(qop,phi);
	return make_mixing_pdf_with_background_functor<MSqPlus, MSqMinus, Time, TimeError>(
								efficiency, Adir, Abar, psi_p, psi_m, psi_i, rcp,
					 			pdf_sigma_t, dalitz_integrator, timeRange, 
					 			f_rnd, f_mistag, f_cmb, pdf_cmb
					 		); 


}



// This function would build a pdf functor (the funtor is a pdf iself!) for time dependent amplitude with fixed Adir and Abar, and normalized pdf(sigma_t).
// Currently, q/p is not considered.
template<Flavor Tag, typename MSqPlus, typename MSqMinus, typename Time, typename TimeError, typename EFFICIENCY, typename ADIR, typename ABAR, typename PDFSIGMAT>
auto time_dependent_rate_with_time_resolution_pdf_bak(hydra::Parameter const& tau, hydra::Parameter const& x, hydra::Parameter const& y, hydra::Parameter const& qop, hydra::Parameter const& phi, hydra::Parameter const& b, hydra::Parameter const& s, EFFICIENCY const& efficiency, ADIR const& Adir, ABAR const& Abar, PDFSIGMAT const& pdf_sigma_t, std::array<double,2> dalitzRangePlus, std::array<double,2> dalitzRangeMinus, std::array<double,2> timeRange, std::array<double,2> timeErrorRange, const long n_calls=50000000)
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

	hydra::Plain<1,  hydra::device::sys_t > sigma_t_integrator({timeErrorRange[0]}, {timeErrorRange[1]}, n_calls);
	double int_pdf_sigma_t = sigma_t_integrator.Integrate(pdf_sigma_t).first;
	std::cout << "Integrate(pdf(sigma_t)) on (timeErrorRange[0], timeErrorRange[1]) = " << int_pdf_sigma_t << std::endl;
	// double int_pdf_sigma_t = pdf_sigma_t.GetNorm();
	// std::cout << "pdf_sigma_t.GetNorm() = " << int_pdf_sigma_t << std::endl;
	if (abs(int_pdf_sigma_t - 1) > 0.01) {
		std::cout << "time_dependent_rate_with_time_resolution_pdf: ";
		std::cout << "The input pdf(sigma_t) should be normalized, but its integration seems not close enough to 1. ";
		std::cout << "Exit!" << std::endl; 
	}

	std::cout << "Pass." << std::endl;

	// integrate on dalitz plane
	hydra::Plain<2,  hydra::device::sys_t > PlainDalitzMC({dalitzRangePlus[0], dalitzRangeMinus[0]}, 
														  {dalitzRangePlus[1], dalitzRangeMinus[1]}, 
														  n_calls);

	// numerically integrate for time-independent part
	// suppose the efficiency plane is time-independent and has not correlation with decay-time or sigmat.

	auto rcp = QoverP<Tag>(qop,phi);

	auto AdirSq = efficiency * rate(Adir);
	auto AbarSq = efficiency * rate(rcp * Abar);
	auto AdirAbarStar_real = efficiency * real_part(Adir * conjugate(rcp * Abar)); 
	auto AdirAbarStar_imag = efficiency * imag_part(Adir * conjugate(rcp * Abar)); 
	
	double AdirSq_int = PlainDalitzMC.Integrate(AdirSq).first;
	double AbarSq_int = PlainDalitzMC.Integrate(AbarSq).first;
	double AdirAbarStar_real_int = PlainDalitzMC.Integrate(AdirAbarStar_real).first;
	double AdirAbarStar_imag_int = PlainDalitzMC.Integrate(AdirAbarStar_imag).first;

	double AsumD2Sq_int = 1. / 4. * ( AdirSq_int + AbarSq_int + 2*AdirAbarStar_real_int );
	double AdiffD2Sq_int = 1. / 4. * ( AdirSq_int + AbarSq_int - 2*AdirAbarStar_real_int );
	hydra::complex<double> AsumD2AdiffD2Star_int = 1. / 4. * hydra::complex<double>( AdirSq_int - AbarSq_int , -2*AdirAbarStar_imag_int );


	// build the final functor

	auto psi_p = MixingPsip<Time, TimeError>(y,tau,s,b);
	auto psi_m = MixingPsim<Time, TimeError>(y,tau,s,b);
	auto psi_i = MixingPsii<Time, TimeError>(x,tau,s,b);

	auto As = Adir + rcp * Abar; 
	auto Ad = conjugate( Adir - rcp * Abar ) ; 

	auto int_on_dalitzplane_decaytime = MixingNormType1<TimeError>(AsumD2Sq_int, AdiffD2Sq_int, AsumD2AdiffD2Star_int, 
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
