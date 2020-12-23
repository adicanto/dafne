#pragma once

#include <hydra/Function.h>
#include <hydra/Complex.h>

namespace dafne {

//========================================================================================
// LASS amplitude as defined in https://cds.cern.ch/record/2230910
//
template<typename R, size_t N=9, typename Signature=hydra::complex<double>(double)>
class LassPropagator: public hydra::BaseFunctor<LassPropagator<R>, Signature, N>
{
	typedef hydra::BaseFunctor<LassPropagator<R>, Signature, N> super_type;
	using super_type::_par;
	
private:
	mutable R _resonance;
	
	__hydra_dual__ inline
	double b1() const { return _par[4]; }

	__hydra_dual__ inline
	double b2() const { return _par[5]; }

	__hydra_dual__ inline
	double b3() const { return _par[6]; }
	
	__hydra_dual__ inline
	double r() const { return _par[7]; }
	
	__hydra_dual__ inline
	double a() const { return _par[8]; }

public:
	LassPropagator() = delete;

	__hydra_dual__
	LassPropagator(hydra::Parameter const &b1, hydra::Parameter const &b2, hydra::Parameter const &b3, hydra::Parameter const &r, hydra::Parameter const &a, R const &resonance) : super_type({resonance.GetParameter(0),resonance.GetParameter(1),resonance.GetParameter(2),resonance.GetParameter(3),b1,b2,b3,r,a}), _resonance(resonance)
	{
		HYDRA_STATIC_ASSERT( N==9, "Number of parameters must be 9.");
		if (_resonance.Spin()!=hydra::SWave)
			HYDRA_EXCEPTION("Lass propagator must be used only for spin-0 resonances");
	}
	
	__hydra_dual__ inline
	LassPropagator<R>& operator=(LassPropagator<R> const &other)
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
		double w = _resonance.RunningWidth(mSqAB);
		double q = _resonance.Q(mSqAB);

		double x = sqrt(mSqAB)/m;
		double f = exp( b1()*x + b2()*x*x + b3()*x*x*x );
		
		double delta = atan2(2.*a()*q , 2.+a()*r()*q*q);
		delta += atan2(m*w , m*m - mSqAB);
		
		return f * sin(delta) * exp(i*delta);
	}
};

} // namespace dafne

