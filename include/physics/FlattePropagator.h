#pragma once

#include <hydra/Function.h>
#include <hydra/Complex.h>

namespace dafne {

template<typename R, size_t N=8, typename Signature=hydra::complex<double>(double)>
class FlattePropagator: public hydra::BaseFunctor<FlattePropagator<R>, Signature, N>
{
	typedef hydra::BaseFunctor<FlattePropagator<R>, Signature, N> super_type;
	using super_type::_par;
	
private:
	mutable R _resonance;

	__hydra_dual__ inline
	double gamma1() const { return _par[4]; }
	
	__hydra_dual__ inline
	double gamma2() const { return _par[5]; }
	
	__hydra_dual__ inline
	double m02a()  const { return _par[6]; }
	
	__hydra_dual__ inline
	double m02b()  const { return _par[7]; }
	
public:
	FlattePropagator() = delete;

	__hydra_dual__
	FlattePropagator(hydra::Parameter const &gamma1, hydra::Parameter const &gamma2, hydra::Parameter const &m02a, hydra::Parameter const &m02b, R const &resonance) : super_type({resonance.GetParameter(0), resonance.GetParameter(1), resonance.GetParameter(2), resonance.GetParameter(3), gamma1, gamma2, m02a, m02b}), _resonance(resonance)
	{
		HYDRA_STATIC_ASSERT( N==8, "Number of parameters must be 8.");
	}
	
	__hydra_dual__ inline
	FlattePropagator<R>& operator=(FlattePropagator<R> const &other)
	{
		if(this==&other) return *this;
		super_type::operator=(other);
		_resonance = other.Resonance();
		return *this;
	}
	
	__hydra_dual__ inline
	R Resonance() const { return _resonance; }
	
	inline
	void Update() 
	{
		_resonance.SetParameter(0, _par[0]);
		_resonance.SetParameter(1, _par[1]);
		_resonance.SetParameter(2, _par[2]);
		_resonance.SetParameter(3, _par[3]);
	}
	
	__hydra_dual__ inline
	hydra::complex<double> Evaluate(double mSqAB) const
	{
		const hydra::complex<double> i(0.,1.);
		
		double m = _resonance.Mass();
		double w = _resonance.Width();

		double g1 = gamma1() * gamma1() * _resonance.Rho(mSqAB) / _resonance.Rho0();

		double rho20 = sqrt( kallen(   m*m, m02a()*m02a(), m02b()*m02b() ) ) / mSqAB;
		double rho2  = sqrt( kallen( mSqAB, m02a()*m02a(), m02b()*m02b() ) ) / mSqAB;
		double g2    = gamma2() * gamma2() * rho2 / rho20;

		double fR = _resonance.BarrierFactor(mSqAB);

		return m*w*gamma1()*gamma1() / ( m*m - mSqAB - i*m*w*(g1+g2)*fR*fR );
	}
};

} // namespace dafne

