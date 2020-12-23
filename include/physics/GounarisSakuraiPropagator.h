#pragma once

#include <hydra/Function.h>
#include <hydra/Complex.h>

namespace dafne {

//========================================================================================
// Gounaris-Sakurai propagator as defined in PRL 21 (1968) 244
// (https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.21.244)
//
template<typename R, size_t N=4, typename Signature=hydra::complex<double>(double)>
class GounarisSakuraiPropagator: public hydra::BaseFunctor<GounarisSakuraiPropagator<R>, Signature, N>
{
	typedef hydra::BaseFunctor<GounarisSakuraiPropagator<R>, Signature, N> super_type;
	using super_type::_par;
	
private:
	mutable R _resonance;
	
	__hydra_dual__ inline
	double gsf( const double& mSqAB ) const
	{
		double w  = _resonance.Width();
		double m  = _resonance.Mass();
		double q0 = _resonance.Q0();
		double q  = _resonance.Q(mSqAB);
		
		double factor = w * m*m / q0;
		double first  = ( (q*q) / (q0*q0) ) * ( gsh(mSqAB) - gsh(m*m) );
		double second = ( m*m - mSqAB ) * gshprime(m*m);
		return factor * ( first + second );
	}
	
	__hydra_dual__ inline
	double gsh( const double& mSqAB ) const
	{
		double q   = _resonance.Q(mSqAB);
		double mAB = sqrt(mSqAB);
		return (2./PI) * (q/mAB) * log( ( mAB+2.*q ) / ( sqrt(_resonance.MSqA())+sqrt(_resonance.MSqB()) ) );
	}

	__hydra_dual__ inline
	double gshprime( const double& mSqAB ) const
	{
		double q = _resonance.Q(mSqAB);
	   double first  = 1. / ( 8. * q*q );
		double second = 1. / ( 2. * mSqAB);
		return ( first - second ) * gsh(mSqAB) + second / PI;
	}

public:
	GounarisSakuraiPropagator() = delete;

	__hydra_dual__
	GounarisSakuraiPropagator(R const &resonance) : super_type( {resonance.GetParameter(0), resonance.GetParameter(1),resonance.GetParameter(2),resonance.GetParameter(3)} ), _resonance(resonance)
	{
		HYDRA_STATIC_ASSERT(N==4, "Number of parameters must be 4.");
	}

	__hydra_dual__ inline
	GounarisSakuraiPropagator<R>& operator=(GounarisSakuraiPropagator<R> const &other)
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
		double m    = _resonance.Mass();
		double w    = _resonance.Width();
		double mSqA = _resonance.MSqA();
		double mSqB = _resonance.MSqB();
		double q0   = _resonance.Q0();
		
		double gsd = ( 3./PI ) * ( (mSqA+mSqB)/2.0/q0/q0 );
		gsd *= log( ( m+2.*q0 ) / ( sqrt(mSqA) + sqrt(mSqB) ) );
		gsd += m / ( 2.*PI*q0 ) - ( (mSqA+mSqB)*m ) / (2.0*PI*q0*q0*q0 );

		const hydra::complex<double> i(0.,1.);
		hydra::complex<double> prop = ( 1. + gsd * w / m ) / ( m*m - mSqAB + gsf(mSqAB) - i*m*_resonance.RunningWidth(mSqAB) );

		return prop;
	}
};

} // namespace dafne
