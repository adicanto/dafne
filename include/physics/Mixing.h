#pragma once

#include <hydra/Function.h>
#include <hydra/Complex.h>
#include <hydra/functions/JohnsonSUShape.h>

#include "RooMathM.h"

namespace dafne {

typedef enum Flavor_t { Negative=-1, Positive=+1 } Flavor;

template<Flavor Tag, typename Signature=hydra::complex<double>(void)>
class QoverP: public hydra::BaseFunctor<QoverP<Tag>, Signature, 2>
{
	typedef hydra::BaseFunctor<QoverP<Tag>, Signature, 2> super_type;
	using super_type::_par;

private:
	__hydra_dual__ inline
	double qop() const { return _par[0]; }

	__hydra_dual__ inline
	double phi() const { return _par[1]; }

public:
	QoverP() = delete;
     
	__hydra_dual__
	QoverP(hydra::Parameter const& qop, hydra::Parameter const& phi) : super_type({qop,phi})
	{}
	
	__hydra_dual__
	QoverP(QoverP<Tag> const& other) : super_type(other)
	{}
	
	__hydra_dual__ inline
	QoverP<Tag>& operator=(QoverP<Tag> const& other)
	{
		if(this==&other) return *this;
		super_type::operator=(other);
		return *this;
	}
	   
	__hydra_dual__ inline
	hydra::complex<double> Evaluate(void ) const
	{
		hydra::complex<double> i(0.,1.);
		return (Tag > 0) ? qop()*exp(i*phi()) : 1./(qop()*exp(i*phi()));
	}
};

template<bool Mixed, typename Time, typename Signature=hydra::complex<double>(Time)>
class MixingG: public hydra::BaseFunctor<MixingG<Mixed,Time>, Signature, 3>
{
	typedef hydra::BaseFunctor<MixingG<Mixed,Time>, Signature, 3> super_type;
	using super_type::_par;

private:
	__hydra_dual__ inline
	double tau() const { return _par[0]; }

	__hydra_dual__ inline
	double x() const { return _par[1]; }

	__hydra_dual__ inline
	double y() const { return _par[2]; }

public:
	MixingG() = delete;
     
	__hydra_dual__
	MixingG(hydra::Parameter const& tau, hydra::Parameter const& x, hydra::Parameter const& y) : super_type({tau,x,y})
	{}
	
	__hydra_dual__
	MixingG(MixingG<Mixed,Time> const& other) : super_type(other)
	{}
	
	__hydra_dual__ inline
	MixingG<Mixed,Time>& operator=( MixingG<Mixed,Time> const& other)
	{
		if(this==&other) return *this;
		super_type::operator=(other);
		return *this;
	}
	   
	__hydra_dual__ inline
	hydra::complex<double> Evaluate(Time time) const
	{
		// Mixing parameters z = -(y+ix);
		hydra::complex<double> z(-y(), -x());

		// t is in units of the average lifetime tau = 1/Gamma = 2/(Gamma_1 + Gamma_2)
		// g_+(t) = exp(-t/2) cosh(zt/2)
		// g_-(t) = exp(-t/2) sinh(zt/2)
		double t = time/tau()/2.;
		return exp(-t)*((Mixed) ? sinh(z*t) : cosh(z*t));
	}
};

template<typename Time>
auto MixingGp(hydra::Parameter const& tau, hydra::Parameter const& x, hydra::Parameter const& y)
{
	return MixingG<false,Time>(tau,x,y);
}

template<typename Time>
auto MixingGm(hydra::Parameter const& tau, hydra::Parameter const& x, hydra::Parameter const& y)
{
	return MixingG<true,Time>(tau,x,y);
}



typedef enum MixingPsiTag_t { PSI_P=+1, PSI_M=-1, PSI_I=0 } MixingPsiTag;

template<typename T>
__hydra_dual__ inline
T _phi(T z)
{
	return exp(std::norm(z)/2.) * 1./2. * RooMathM::erfc(-z/sqrt(2.));
}

template<typename T>
__hydra_dual__ inline
T _psi(double x, T kappa)
{
	return exp((-2*x*kappa+kappa*kappa)/2.) * 1./2. * (1. + RooMathM::erf((x - kappa)/sqrt(2.)));
}

template<typename T>
__hydra_dual__ inline
T _int_psi_dt(double sigma, double x, T kappa)
{
	return sigma/kappa * 0.5 * (RooMathM::erf(x/sqrt(2.)) - 2.*_psi(x,kappa)) ;
}

template<MixingPsiTag Tag, typename Time, typename TimeError, typename Signature=hydra::complex<double>(Time,TimeError)>
class MixingPsi: public hydra::BaseFunctor<MixingPsi<Tag,Time,TimeError>, Signature, 4>
{
	typedef hydra::BaseFunctor<MixingPsi<Tag,Time,TimeError>, Signature, 4> super_type;
	using super_type::_par;

private:
	__hydra_dual__ inline
	double suffix() const { return _par[0]; }

	__hydra_dual__ inline
	double tau() const { return _par[1]; }

	__hydra_dual__ inline
	double Gamma() const { return 1./_par[1]; }

	__hydra_dual__ inline
	double s() const { return _par[2]; }

	__hydra_dual__ inline
	double b() const { return _par[3]; }

public:
	MixingPsi() = delete;
     
	__hydra_dual__
	MixingPsi(hydra::Parameter const& suffix, hydra::Parameter const& tau, hydra::Parameter const& s, hydra::Parameter const& b) : super_type({suffix,tau,s,b})
	{}
	
	__hydra_dual__
	MixingPsi(MixingPsi<Tag,Time,TimeError> const& other) : super_type(other)
	{}
	
	__hydra_dual__ inline
	MixingPsi<Tag,Time,TimeError>& operator=( MixingPsi<Tag,Time,TimeError> const& other)
	{
		if(this==&other) return *this;
		super_type::operator=(other);
		return *this;
	}

	__hydra_dual__ inline
	std::complex<double> Evaluate(Time t, TimeError sigma_t) const
	{
		double chi = (t-b()) / (s()*sigma_t);
		if (Tag == PSI_P) {
			double y = suffix();
			double kappa_p = (1+y)*Gamma()*s()*sigma_t;
			double result = _psi(chi, kappa_p); // calling _psi<double>()
			return result; // convert double to std::complex<double> automatically, when returning.

		} else if (Tag == PSI_M) {
			double y = suffix();
			double kappa_m = (1-y)*Gamma()*s()*sigma_t;
			double result = _psi(chi, kappa_m);
			return result;

		} else if (Tag == PSI_I) {
			double x = suffix();
			std::complex<double> kappa_i(Gamma()*s()*sigma_t, x*Gamma()*s()*sigma_t); 
			std::complex<double> result = _psi(chi, kappa_i);
			return result;
		}

		std::cout << "MixingPsi: " << "The Tag(" << Tag << ") is not PSI_P(" << PSI_P << "), PSI_M(" << PSI_M << ") or PSI_I(" << PSI_I << ")" << std::endl;  
		std::cout << "Something must be wrong! Exit!" << std::endl;
		exit(-1);
		return -9999;
	}

	hydra::complex<double> AnalyticalIntegral(Time t0, Time t1, TimeError sigma_t) const 
	{
		double sigma = s()*sigma_t;
		double chi0 = (t0-b()) / (sigma);
		double chi1 = (t1-b()) / (sigma);


		if (Tag == PSI_P) {
			double y = suffix();
			double kappa_p = (1+y)*Gamma()*sigma;
			double result = _int_psi_dt(sigma, chi1, kappa_p) - _int_psi_dt(sigma, chi0, kappa_p);
			return result;
		} else if (Tag == PSI_M) {
			double y = suffix();
			double kappa_m = (1-y)*Gamma()*sigma;
			double result = (_int_psi_dt(sigma, chi1, kappa_m) - _int_psi_dt(sigma, chi0, kappa_m));
			return result;
		} else if (Tag == PSI_I) {
			double x = suffix();
			std::complex<double> kappa_i(Gamma()*sigma, x*Gamma()*sigma); 
			hydra::complex<double> result = _int_psi_dt(sigma, chi1, kappa_i) - _int_psi_dt(sigma, chi0, kappa_i);
			return result;
		}

		std::cout << "MixingPsi: " << "The Tag(" << Tag << ") is not PSI_P(" << PSI_P << "), PSI_M(" << PSI_M << ") or PSI_I(" << PSI_I << ")" << std::endl;  
		std::cout << "Something must be wrong! Exit!" << std::endl;
		exit(-1);
		return -9999;
	}
};

template<typename Time, typename TimeError>
auto MixingPsip(hydra::Parameter const& suffix, hydra::Parameter const& tau, hydra::Parameter const& s, hydra::Parameter const& b)
{
	return MixingPsi<PSI_P,Time,TimeError>(suffix,tau,s,b);
}

template<typename Time, typename TimeError>
auto MixingPsim(hydra::Parameter const& suffix, hydra::Parameter const& tau, hydra::Parameter const& s, hydra::Parameter const& b)
{
	return MixingPsi<PSI_M,Time,TimeError>(suffix,tau,s,b);
}

template<typename Time, typename TimeError>
auto MixingPsii(hydra::Parameter const& suffix, hydra::Parameter const& tau, hydra::Parameter const& s, hydra::Parameter const& b)
{
	return MixingPsi<PSI_I,Time,TimeError>(suffix,tau,s,b);
}


// the analytical normalization functor for time_dependent_rate_with_time_resolution_pdf_type1
template<typename TimeError, typename Signature=double(TimeError)>
class MixingNormType1: public hydra::BaseFunctor<MixingNormType1<TimeError>, Signature, 5>
{
	typedef hydra::BaseFunctor<MixingNormType1<TimeError>, Signature, 5> super_type;
	using super_type::_par;

private:
	__hydra_dual__ inline
	double x() const { return _par[0]; }

	__hydra_dual__ inline
	double y() const { return _par[1]; }

	__hydra_dual__ inline
	double tau() const { return _par[2]; }

	__hydra_dual__ inline
	double Gamma() const { return 1./_par[2]; }

	__hydra_dual__ inline
	double s() const { return _par[3]; }

	__hydra_dual__ inline
	double b() const { return _par[4]; }

public:
	MixingNormType1() = delete;
     
	__hydra_dual__
	MixingNormType1(const double AsumD2Sq_int, const double AdiffD2Sq_int, const std::complex<double> AsumD2AdiffD2Star_int, hydra::Parameter const& x, hydra::Parameter const& y, hydra::Parameter const& tau, hydra::Parameter const& s, hydra::Parameter const& b, std::array<double,2> const& timeRange) : super_type({x,y,tau,s,b}), fAsumD2Sq_int(AsumD2Sq_int), fAdiffD2Sq_int(AdiffD2Sq_int), fAsumD2AdiffD2Star_int(AsumD2AdiffD2Star_int), _timeRange(timeRange)
	{}
	
	__hydra_dual__
	MixingNormType1(MixingNormType1<TimeError> const& other) : super_type(other), fAsumD2Sq_int(other.fAsumD2Sq_int), fAdiffD2Sq_int(other.fAdiffD2Sq_int), fAsumD2AdiffD2Star_int(other.fAsumD2AdiffD2Star_int), _timeRange(other._timeRange)
	{}
	
	__hydra_dual__ inline
	MixingNormType1<TimeError>& operator=( MixingNormType1<TimeError> const& other)
	{
		if(this==&other) return *this;
		super_type::operator = (other);
		fAsumD2Sq_int = other.fAsumD2Sq_int; 
		fAdiffD2Sq_int = other.fAdiffD2Sq_int; 
		fAsumD2AdiffD2Star_int = other.fAsumD2AdiffD2Star_int;
		_timeRange = other._timeRange;
		return *this;
	}

	__hydra_dual__ inline
	double Evaluate(TimeError sigma_t) const 
	// The dalitz variables and decay-time are integrated, leaving only the sigma_t
	{

		double sigma = s()*sigma_t;
		double chi0 = (_timeRange[0]-b()) / (sigma);
		double chi1 = (_timeRange[1]-b()) / (sigma);

		double kappa_p = (1+y())*Gamma()*sigma;
		double kappa_m = (1-y())*Gamma()*sigma;
		std::complex<double> kappa_i(Gamma()*sigma, x()*Gamma()*sigma); 
		
		double first_term = fAsumD2Sq_int * (_int_psi_dt(sigma, chi1, kappa_p) - _int_psi_dt(sigma, chi0, kappa_p));
		double second_term = fAdiffD2Sq_int * (_int_psi_dt(sigma, chi1, kappa_m) - _int_psi_dt(sigma, chi0, kappa_m));
		double third_term = 2 * (fAsumD2AdiffD2Star_int * (_int_psi_dt(sigma, chi1, kappa_i) - _int_psi_dt(sigma, chi0, kappa_i))).real();

		return first_term + second_term + third_term;
	}

private:
	// member variables to store the integration for time dependent parts for the first, second and third term
	double fAsumD2Sq_int;
	double fAdiffD2Sq_int;
	std::complex<double> fAsumD2AdiffD2Star_int;

	std::array<double,2> _timeRange;
};


// We can not use:
// auto johnson_su = hydra::make_pdf(
// 					hydra::JohnsonSU<DecayTimeError>(johnson_gamma, johnson_delta, johnson_xi, johnson_lambda), 
// 					hydra::AnalyticalIntegral<hydra::JohnsonSU<DecayTimeError>>(phsp.TimeErrorMin(), phsp.TimeErrorMax()) 
// 					);
// when building the T(mSqP, mSqM, Time, TimeError), because the compose() could not take Pdf class as input
template<typename ArgType, typename Signature=double(ArgType)>
class NormalizedJohnsonSU: public hydra::BaseFunctor<NormalizedJohnsonSU<ArgType>, Signature, 4>
{
	using hydra::BaseFunctor<NormalizedJohnsonSU<ArgType>, Signature, 4>::_par;

public:

	NormalizedJohnsonSU()=delete;

	NormalizedJohnsonSU(hydra::Parameter const& gamma, hydra::Parameter const& delta, 
			hydra::Parameter const& xi, hydra::Parameter const& lambda, double const LowerLimit, double UpperLimit):
			hydra::BaseFunctor<NormalizedJohnsonSU<ArgType>, Signature, 4>({gamma, delta, xi, lambda}),
			fLowerLimit(LowerLimit), fUpperLimit(UpperLimit)
			{}

	__hydra_host__ __hydra_device__
	NormalizedJohnsonSU(NormalizedJohnsonSU<ArgType> const& other ):
			hydra::BaseFunctor<NormalizedJohnsonSU<ArgType>, Signature, 4>(other), 
			fLowerLimit(other.fLowerLimit), fUpperLimit(other.fUpperLimit)
			{}

	__hydra_host__ __hydra_device__
	NormalizedJohnsonSU<ArgType>&
	operator=(NormalizedJohnsonSU<ArgType> const& other ){
		if(this==&other) return  *this;
		hydra::BaseFunctor<NormalizedJohnsonSU<ArgType>, Signature, 4>::operator=(other);
		fLowerLimit = other.fLowerLimit;
		fUpperLimit = other.fUpperLimit;
		return  *this;
	}


	inline double cumulative( const double gamma,  const double delta,  const double xi,  const double lambda, const double x) const
	{
		//actually only 1/lambda is used
		double inverse_lambda = 1.0/lambda;

		// z =  (x-xi)/lambda
		double z    = (x-xi)*inverse_lambda;

		// C = {(\gamma + \delta * \asinh(z) )}
		double C = gamma + delta*::asinh(z);

		return 0.5*(1.0 + ::erf(C*hydra::math_constants::inverse_sqrt2));
	}


	__hydra_host__ __hydra_device__
	inline double Evaluate(ArgType x)  const
	{
		//gathering parameters
		double gamma  = _par[0];
		double delta  = _par[1];
		double xi     = _par[2];
		//actually only 1/lambda is used
		double inverse_lambda = 1.0/_par[3];

		// z =  (x-xi)/lambda
		double z    = (x-xi)*inverse_lambda;

		// A = \frac{\delta}{ \lambda * \sqrt{2\pi} }
		double A = inverse_lambda*delta*hydra::math_constants::inverse_sqrt2Pi;

		//B = \frac{1}{\sqrt{1 + z^2}}
		double B = 1.0/::sqrt( 1 + z*z);

		// C = {(\gamma + \delta * \asinh(z) )}^{2}
		double C = gamma + delta*::asinh(z); C *=C;

		double r = cumulative(_par[0], _par[1], _par[2], _par[3], fUpperLimit)
 				 - cumulative(_par[0], _par[1], _par[2], _par[3], fLowerLimit);

		double result = A*B*::exp(-0.5*C) / r;

		return result;

		// return CHECK_VALUE(result, "par[0]=%f, par[1]=%f, par[2]=%f, par[3]=%f", _par[0], _par[1], _par[2], _par[3]  );


	}

private:
	double fLowerLimit;
	double fUpperLimit;



};

// Psi0 and integration of Psi0 for the time dependent background
template<typename Time, typename TimeError, typename Signature=double(Time,TimeError)>
class Psi0: public hydra::BaseFunctor<Psi0<Time,TimeError>, Signature, 3>
{
	typedef hydra::BaseFunctor<Psi0<Time,TimeError>, Signature, 3> super_type;
	using super_type::_par;

private:
	__hydra_dual__ inline
	double tau() const { return _par[0]; }

	__hydra_dual__ inline
	double Gamma() const { return 1./_par[0]; }

	__hydra_dual__ inline
	double s() const { return _par[1]; }

	__hydra_dual__ inline
	double b() const { return _par[2]; }

public:
	Psi0() = delete;
     
	__hydra_dual__
	Psi0(hydra::Parameter const& tau, hydra::Parameter const& s, hydra::Parameter const& b) : super_type({tau,s,b})
	{}
	
	__hydra_dual__
	Psi0(Psi0<Time,TimeError> const& other) : super_type(other)
	{}
	
	__hydra_dual__ inline
	Psi0<Time,TimeError>& operator=( Psi0<Time,TimeError> const& other)
	{
		if(this==&other) return *this;
		super_type::operator=(other);
		return *this;
	}

	__hydra_dual__ inline
	double Evaluate(Time t, TimeError sigma_t) const
	{
		double chi = (t-b()) / (s()*sigma_t);
		double kappa_0 = (1+0)*Gamma()*s()*sigma_t;
		double result = _psi(chi, kappa_0); 
		return result; 
	}
};

// analytical integration of Psi0
template<typename TimeError, typename Signature=double(TimeError)>
class Psi0Integration: public hydra::BaseFunctor<Psi0Integration<TimeError>, Signature, 3>
{
	typedef hydra::BaseFunctor<Psi0Integration<TimeError>, Signature, 3> super_type;
	using super_type::_par;

private:
	__hydra_dual__ inline
	double tau() const { return _par[0]; }

	__hydra_dual__ inline
	double Gamma() const { return 1./_par[0]; }

	__hydra_dual__ inline
	double s() const { return _par[1]; }

	__hydra_dual__ inline
	double b() const { return _par[2]; }

public:
	Psi0Integration() = delete;
     
	__hydra_dual__
	Psi0Integration(hydra::Parameter const& tau, hydra::Parameter const& s, hydra::Parameter const& b, std::array<double,2> const& timeRange) : super_type({tau,s,b}), _timeRange(timeRange)
	{}
	
	__hydra_dual__
	Psi0Integration(Psi0Integration<TimeError> const& other) : super_type(other), _timeRange(other._timeRange)
	{}
	
	__hydra_dual__ inline
	Psi0Integration<TimeError>& operator=( Psi0Integration<TimeError> const& other)
	{
		if(this==&other) return *this;
		super_type::operator=(other);
		_timeRange = other._timeRange;
		return *this;
	}

	__hydra_dual__ inline
	double Evaluate(TimeError sigma_t) const
	{
		double sigma = s()*sigma_t;

		double chi0 = (_timeRange[0]-b()) / sigma;
		double chi1 = (_timeRange[1]-b()) / sigma;

		double kappa_0 = (1+0)*Gamma()*sigma;
		double result = _int_psi_dt(sigma, chi1, kappa_0) - _int_psi_dt(sigma, chi0, kappa_0); 
		return result; 
	}

private:
	std::array<double,2> _timeRange;
};





template<typename MSqPlus, typename MSqMinus, typename Time, typename TimeError, typename EFFICIENCY, typename ADIR, typename ABAR, typename PSI_P, typename PSI_M, typename PSI_I, typename RCP, typename PDFSIGMAT, typename INTEGRATOR, typename Signature=double(MSqPlus, MSqMinus, Time, TimeError)>
class MixingPdfFunctor: public hydra::BaseCompositeFunctor<
							MixingPdfFunctor<MSqPlus, MSqMinus, Time, TimeError, EFFICIENCY, ADIR, ABAR, PSI_P, PSI_M, PSI_I, RCP,  PDFSIGMAT, INTEGRATOR>, // typename Composite
							hydra_thrust::tuple<EFFICIENCY, ADIR, ABAR, PSI_P, PSI_M, PSI_I, RCP, PDFSIGMAT>, // typename FunctorList
								Signature>
{

	typedef hydra::BaseCompositeFunctor<
				MixingPdfFunctor<MSqPlus, MSqMinus, Time, TimeError, EFFICIENCY, ADIR, ABAR, PSI_P, PSI_M, PSI_I, RCP, PDFSIGMAT, INTEGRATOR>, // typename Composite
				hydra_thrust::tuple<EFFICIENCY, ADIR, ABAR, PSI_P, PSI_M, PSI_I, RCP, PDFSIGMAT>, // typename FunctorList
				Signature> super_type;


public:

	MixingPdfFunctor()=delete;

	MixingPdfFunctor(EFFICIENCY const& efficiency, ADIR const& Adir, ABAR const& Abar, 
					 PSI_P const& psi_p, PSI_M const& psi_m, PSI_I const& psi_i, RCP const& rcp,
					 PDFSIGMAT const& pdf_sigma_t, INTEGRATOR const& integrator, std::array<double,2> timeRange): 
	super_type(efficiency, Adir, Abar, psi_p, psi_m, psi_i, rcp, pdf_sigma_t),
	fIntegrator(integrator),
	fAsumD2SqNormCache(std::unordered_map<size_t, double>() ),
	fAdiffD2SqNormCache(std::unordered_map<size_t, double>() ),
	fAsumD2AdiffD2StarNormCache(std::unordered_map<size_t, hydra::complex<double>>() ),
	fTimeRange(timeRange)
	{NormalizeDalitzPlane();}

	__hydra_host__ __hydra_device__
	inline MixingPdfFunctor(MixingPdfFunctor<MSqPlus, MSqMinus, Time, TimeError, EFFICIENCY, ADIR, ABAR, PSI_P, PSI_M, PSI_I, RCP, PDFSIGMAT, INTEGRATOR> const& other):
    super_type(other),
    fIntegrator(other.fIntegrator),
    fAsumD2SqNormCache(other.fAsumD2SqNormCache),
    fAdiffD2SqNormCache(other.fAdiffD2SqNormCache),
    fAsumD2AdiffD2StarNormCache(other.fAsumD2AdiffD2StarNormCache),
	fTimeRange(other.fTimeRange)
	{NormalizeDalitzPlane();}

	__hydra_host__ __hydra_device__
	inline MixingPdfFunctor<MSqPlus, MSqMinus, Time, TimeError, EFFICIENCY, ADIR, ABAR, PSI_P, PSI_M, PSI_I, RCP, PDFSIGMAT, INTEGRATOR>& operator=(MixingPdfFunctor<MSqPlus, MSqMinus, Time, TimeError, EFFICIENCY, ADIR, ABAR, PSI_P, PSI_M, PSI_I, RCP, PDFSIGMAT, INTEGRATOR> const& other)
	{
		if(this==&other) return *this;
		super_type::operator=(other);
		fIntegrator=other.fIntegrator;
	    fAsumD2SqNormCache=other.fAsumD2SqNormCache;
	    fAdiffD2SqNormCache=other.fAdiffD2SqNormCache;
	    fAsumD2AdiffD2StarNormCache=other.fAsumD2AdiffD2StarNormCache;
		fTimeRange=other.fTimeRange;
		return *this;
	}

	inline size_t  GetDalitzParametersKeyC() const { // a const version of GetDalitzParametersKeyC to be called in NormalizeDalitzPlane() const

		std::vector<hydra::Parameter*> _parameters;
		auto _Adir = hydra::get<1>(this->GetFunctors());
		_Adir.AddUserParameters(_parameters);

		std::vector<double> _temp(_parameters.size());

		for(size_t i=0; i< _parameters.size(); i++)
			_temp[i]= *(_parameters[i]);

		size_t key = hydra::detail::hash_range(_temp.begin(), _temp.end() );

		return key;
	}

	// I guess that the result of normalization could not be shared between threads, therefore 
	// maybe the functor needs to normalize inside each thread. If the Evaluate() of the functor
	// could be called from the host right before getting into device, in each loop, then, this
	// problem could be solved (the host would normalize it, then store some information in
	// fNormCache).
	// One possible method is to follow the manner of hydra::SpilineFunctor: the functor only has
	// iterator members, and the device side iterables are in the main program.
	inline	void NormalizeDalitzPlane( ) const
	{
		auto efficiency =  hydra::get<0>(this->GetFunctors());
		auto Adir =  hydra::get<1>(this->GetFunctors());
		auto Abar =  hydra::get<2>(this->GetFunctors());
		auto rcp =  hydra::get<6>(this->GetFunctors());

		size_t key = GetDalitzParametersKeyC();

		double AsumD2SqNorm;
	    double AdiffD2SqNorm;
	    hydra::complex<double> AsumD2AdiffD2StarNorm;

		auto search1 = fAsumD2SqNormCache.find(key);
		if (search1 != fAsumD2SqNormCache.end() && fAsumD2SqNormCache.size()>0) {
			//std::cout << "found in cache "<< key << std::endl;
			AsumD2SqNorm = search1->second;

			auto search2 = fAdiffD2SqNormCache.find(key);
			if (search2 != fAdiffD2SqNormCache.end() && fAdiffD2SqNormCache.size()>0)
				AdiffD2SqNorm = search2->second;
			else {
				std::cout << "Cached AsumD2Sq and not Cached AdiffD2Sq, something must be wrong!" << std::endl;
				exit(-1);
			}

			auto search3 = fAsumD2AdiffD2StarNormCache.find(key);
			if (search3 != fAsumD2AdiffD2StarNormCache.end() && fAsumD2AdiffD2StarNormCache.size()>0)
				AsumD2AdiffD2StarNorm = search3->second;
			else {
				std::cout << "Cached AsumD2Sq and not Cached AsumD2AdiffD2Star, something must be wrong!" << std::endl;
				exit(-1);
			}

		} else {
			auto AdirSq = efficiency * rate(Adir);
			auto AbarSq = efficiency * rate(rcp * Abar);
			auto AdirAbarStar_real = efficiency * real_part(Adir * conjugate(rcp * Abar)); 
			auto AdirAbarStar_imag = efficiency * imag_part(Adir * conjugate(rcp * Abar)); 
	
			double AdirSq_int = fIntegrator(AdirSq).first;
			double AbarSq_int = fIntegrator(AbarSq).first;
			double AdirAbarStar_real_int = fIntegrator(AdirAbarStar_real).first;
			double AdirAbarStar_imag_int = fIntegrator(AdirAbarStar_imag).first;

			double AsumD2Sq_int = 1. / 4. * ( AdirSq_int + AbarSq_int + 2*AdirAbarStar_real_int );
			double AdiffD2Sq_int = 1. / 4. * ( AdirSq_int + AbarSq_int - 2*AdirAbarStar_real_int );
			hydra::complex<double> AsumD2AdiffD2Star_int = 1. / 4. * hydra::complex<double>( AdirSq_int - AbarSq_int , -2*AdirAbarStar_imag_int );


			AsumD2SqNorm = AsumD2Sq_int;
		    AdiffD2SqNorm = AdiffD2Sq_int;
		    AsumD2AdiffD2StarNorm = AsumD2AdiffD2Star_int;

			fAsumD2SqNormCache[key] = AsumD2SqNorm;
			fAdiffD2SqNormCache[key] = AdiffD2SqNorm;
			fAsumD2AdiffD2StarNormCache[key] = AsumD2AdiffD2StarNorm;
		}


		fAsumD2SqNormThis = AsumD2SqNorm;
		fAdiffD2SqNormThis = AdiffD2SqNorm;
		fAsumD2AdiffD2StarNormThis = AsumD2AdiffD2StarNorm;

	}


	 template<typename ...T>
  	__hydra_host__ __hydra_device__
  	inline typename  super_type::return_type Evaluate(MSqPlus m2p, MSqMinus m2m, Time t, TimeError sigma_t) const
  	{
  		NormalizeDalitzPlane();

  		auto efficiency =  hydra::get<0>(this->GetFunctors());
  		auto Adir =  hydra::get<1>(this->GetFunctors());
		auto Abar =  hydra::get<2>(this->GetFunctors());
		auto psi_p = hydra::get<3>(this->GetFunctors());
		auto psi_m = hydra::get<4>(this->GetFunctors());
		auto psi_i = hydra::get<5>(this->GetFunctors());
		auto rcp =  hydra::get<6>(this->GetFunctors());
		auto pdf_sigma_t =  hydra::get<7>(this->GetFunctors());

		double _efficiency = efficiency(m2p, m2m);
		hydra::complex<double> _Adir = Adir(m2p, m2m);
		hydra::complex<double> _Abar = Abar(m2p, m2m);
		double _psi_p = (psi_p(t, sigma_t)).real();
		double _psi_m = (psi_m(t, sigma_t)).real();
		hydra::complex<double> _psi_i = psi_i(t, sigma_t);
		hydra::complex<double> _rcp = rcp();
		double _pdf_sigma_t = pdf_sigma_t(sigma_t);

		hydra::complex<double> _As = _Adir + _rcp * _Abar; 
		hydra::complex<double> _Ad = hydra::conj( _Adir - _rcp * _Abar ) ; 

		double t0 = fTimeRange[0];
		double t1 = fTimeRange[1];
		double first_term = fAsumD2SqNormThis * (psi_p.AnalyticalIntegral(t0, t1, sigma_t).real());
		double second_term = fAdiffD2SqNormThis * (psi_m.AnalyticalIntegral(t0, t1, sigma_t).real());
		double third_term = 2 * (fAsumD2AdiffD2StarNormThis * psi_i.AnalyticalIntegral(t0, t1, sigma_t)).real();
		double _int_on_dalitzplane_decaytime = first_term + second_term + third_term;

		double result = _efficiency * ( hydra::norm(_As/2.)*_psi_p + hydra::norm(_Ad/2.)*_psi_m + (_As*_Ad/2.*_psi_i).real() ) * _pdf_sigma_t / _int_on_dalitzplane_decaytime;

  		return result;
  	}


private:
	mutable INTEGRATOR fIntegrator;
	// double fNorm; // fNorm is a father class member
	mutable double fAsumD2SqNormThis;  // define as mutable to be changed in the Evaluate
	mutable std::unordered_map<size_t, double> fAsumD2SqNormCache;
	mutable double fAdiffD2SqNormThis;  // define as mutable to be changed in the Evaluate
	mutable std::unordered_map<size_t, double> fAdiffD2SqNormCache;
	mutable hydra::complex<double> fAsumD2AdiffD2StarNormThis;  // define as mutable to be changed in the Evaluate
	mutable std::unordered_map<size_t, hydra::complex<double>> fAsumD2AdiffD2StarNormCache;
	std::array<double,2> fTimeRange;
};


template<typename MSqPlus, typename MSqMinus, typename Time, typename TimeError, typename EFFICIENCY, typename ADIR, typename ABAR, typename PSI_P, typename PSI_M, typename PSI_I, typename RCP, typename PDFSIGMAT, typename INTEGRATOR>
auto make_mixing_pdf_functor(EFFICIENCY const& efficiency, ADIR const& Adir, ABAR const& Abar, 
					 PSI_P const& psi_p, PSI_M const& psi_m, PSI_I const& psi_i, RCP const& rcp,
					 PDFSIGMAT const& pdf_sigma_t, INTEGRATOR const& integrator, std::array<double,2> timeRange)
{
	return MixingPdfFunctor<MSqPlus, MSqMinus, Time, TimeError, 
							EFFICIENCY, ADIR, ABAR, PSI_P, PSI_M, PSI_I, 
							RCP, PDFSIGMAT, INTEGRATOR>(
								efficiency, Adir, Abar, psi_p, psi_m, psi_i, 
								rcp, pdf_sigma_t, integrator, timeRange);
}



} // namespace dafne



// The RngFormula of NormalizedJohnsonSU and JohnsonSU should be the same
template<typename ArgType>
struct hydra::RngFormula< dafne::NormalizedJohnsonSU<ArgType> >
{

	typedef ArgType value_type;
	__hydra_host__ __hydra_device__
	inline unsigned NCalls( dafne::NormalizedJohnsonSU<ArgType>const&) const
	{
		return 1;
	}

	template< typename T>
	__hydra_host__ __hydra_device__
	inline unsigned NCalls( std::initializer_list<T>) const
	{
		return 1;
	}


	template<typename Engine>
	__hydra_host__ __hydra_device__
	inline value_type Generate(Engine& rng, dafne::NormalizedJohnsonSU<ArgType>const& functor) const
	{
		double gamma  = functor[0];
		double delta  = functor[1];
		double xi     = functor[2];
		double lambda = functor[3];

		return static_cast<value_type>(::sinh((nci(RngBase::uniform(rng)) -gamma)/delta )*lambda + xi);
	}

	template<typename Engine, typename T>
	__hydra_host__ __hydra_device__
	inline 	value_type Generate(Engine& rng, std::initializer_list<T> pars) const
	{
		double gamma  = pars[0];
		double delta  = pars[1];
		double xi     = pars[2];
		double lambda = pars[3];

		return static_cast<value_type>(::sinh((nci(RngBase::uniform(rng)) -gamma)/delta )*lambda + xi);

	}
private:
	__hydra_host__ __hydra_device__
	inline double nci(double x) const
		{
			static const double sqrt_two         = 1.4142135623730950488017;

			return sqrt_two *(hydra::erfinv(2*x-1));
		}

};




