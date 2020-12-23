#pragma once

#include <hydra/Function.h>
#include <hydra/Complex.h>

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

} // namespace dafne
