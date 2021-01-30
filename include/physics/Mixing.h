#pragma once

#include <hydra/Function.h>
#include <hydra/Complex.h>

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
			std::complex<double> kappa_i(Gamma()*s()*sigma_t, -x*Gamma()*s()*sigma_t); 
			std::complex<double> result = _psi(chi, kappa_i);
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




} // namespace dafne
