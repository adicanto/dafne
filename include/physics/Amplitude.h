#pragma once

#include <hydra/Function.h>
#include <hydra/Complex.h>

#include <physics/BaseAmplitude.h>

namespace dafne {

template<typename Arg1, typename Arg2, size_t N, typename P, typename Signature=hydra::complex<double>(Arg1,Arg2)>
class Amplitude: public BaseAmplitude, public hydra::BaseFunctor<Amplitude<Arg1,Arg2,N,P>, Signature, 2+N>
{
	typedef hydra::BaseFunctor<Amplitude<Arg1,Arg2,N,P>, Signature, 2+N> super_type;
	using super_type::_par;

private:
	mutable P _propagator;
	
	__hydra_dual__ inline
	hydra::complex<double> coeff() const { return hydra::complex<double>(_par[0],_par[1]); }

public:
	Amplitude() = delete;

	__hydra_dual__
	Amplitude(hydra::Parameter const &c_re, hydra::Parameter const &c_im, P const &propagator, const char *name="", const char *label="") : BaseAmplitude(name,label), super_type(), _propagator(propagator)
	{
		if (N != _propagator.GetNumberOfParameters())
			HYDRA_EXCEPTION("N must match the number of propagator parameters");

		std::array<hydra::Parameter,2+N> par;
		par[0] = c_re;
		par[1] = c_im;
		for (size_t i=0; i<N; ++i) par[2+i] = _propagator.GetParameter(i);
		_par = super_type(par);
	}

	__hydra_dual__
	Amplitude(Amplitude<Arg1,Arg2,N,P> const& other) : BaseAmplitude(other), super_type(other), _propagator(other.Propagator())
	{}

	__hydra_dual__ inline
	Amplitude<Arg1,Arg2,N,P>& operator=(Amplitude<Arg1,Arg2,N,P> const& other)
	{
		if(this==&other) return *this;
		BaseAmplitude::operator=(other);
		super_type::operator=(other);
		_propagator = other.Propagator();
		return *this;
	}

	__hydra_dual__ inline
	P Propagator() const { return _propagator; }

	__hydra_dual__ inline
	void FixParameters()
	{
		for (size_t i=0; i<2+N; ++i) this->Parameter(i).SetFixed(true);
	}

	__hydra_dual__ inline
	void ReleaseParameters()
	{
		for (size_t i=0; i<2+N; ++i) this->Parameter(i).SetFixed(false);
	}

	__hydra_dual__ inline
	ThreeBodyPhaseSpace PhaseSpace() const { return _propagator.Resonance().PhaseSpace(); }

	inline
	void Update()
	{
		for (size_t i=0; i<N; ++i) _propagator.SetParameter(i, _par[2+i]);
	}

	__hydra_dual__ inline
	hydra::complex<double> Evaluate(Arg1 mSq12, Arg2 mSq13) const
	{
		if (IsRemoved()) return hydra::complex<double>(0.,0.);

		if (coeff() == hydra::complex<double>(0.,0.)) return hydra::complex<double>(0.,0.);

		auto r = _propagator.Resonance();

		if ( !r.InPhaseSpace(mSq12,mSq13) ) return hydra::complex<double>(0.,0.);

		auto mSqAB = r.MSqAB(mSq12,mSq13);
		auto mSqAC = r.MSqAC(mSq12,mSq13);
		auto mSqBC = r.MSqBC(mSq12,mSq13);

		return coeff() * r.Centrifugal(mSqAB) * r.Zemach(mSqAB,mSqAC,mSqBC) * _propagator(mSqAB);
	}

};

} // namespace dafne
