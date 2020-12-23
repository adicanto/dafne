#pragma once

#include <hydra/Function.h>
#include <hydra/Complex.h>

namespace dafne {

template<typename R, size_t N=4, typename Signature=hydra::complex<double>(double)>
class BreitWignerPropagator: public hydra::BaseFunctor<BreitWignerPropagator<R>, Signature, N>
{
	typedef hydra::BaseFunctor<BreitWignerPropagator<R>, Signature, N> super_type;
	using super_type::_par;
	
private:
	mutable R _resonance;

public:
	BreitWignerPropagator() = delete;

	__hydra_dual__
	BreitWignerPropagator(R const &resonance) : super_type( {resonance.GetParameter(0), resonance.GetParameter(1),resonance.GetParameter(2),resonance.GetParameter(3)} ), _resonance(resonance)
	{
		HYDRA_STATIC_ASSERT(N==4, "Number of parameters must be 4.");
	}

	__hydra_dual__ inline
	BreitWignerPropagator<R>& operator=(BreitWignerPropagator<R> const &other)
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
		double m = _resonance.Mass();
		double w = _resonance.RunningWidth(mSqAB);
		return 1./hydra::complex<double>(m*m - mSqAB, -m*w);
	}
};

} // namespace dafne
