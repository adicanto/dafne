#pragma once

#include <hydra/Function.h>
#include <hydra/Complex.h>

#include <physics/BaseAmplitude.h>
#include <physics/ThreeBodyPhaseSpace.h>

namespace dafne {

template<typename Arg1, typename Arg2, typename Signature=hydra::complex<double>(Arg1,Arg2)>
class NonResonant: public BaseAmplitude, public hydra::BaseFunctor<NonResonant<Arg1,Arg2>, Signature, 2>
{
	typedef hydra::BaseFunctor<NonResonant<Arg1,Arg2>, Signature, 2> super_type;
	using super_type::_par;

private:
	ThreeBodyPhaseSpace _phsp;
	
public:
	NonResonant() = delete;

	__hydra_dual__
	NonResonant(hydra::Parameter const& c_re, hydra::Parameter const& c_im, ThreeBodyPhaseSpace const &phsp, const char *name="", const char *label="") : BaseAmplitude(name,label), super_type({c_re, c_im}), _phsp(phsp)
	{}

	__hydra_dual__
	NonResonant(NonResonant<Arg1,Arg2> const& other) : BaseAmplitude(other), super_type(other), _phsp(other.PhaseSpace())
	{}

	__hydra_dual__
	NonResonant<Arg1,Arg2>& operator=(NonResonant<Arg1,Arg2> const& other)
	{
		if(this==&other) return *this;
		BaseAmplitude::operator=(other);
		super_type::operator=(other);
		_phsp = other.PhaseSpace();
		return *this;
	}

	__hydra_dual__ inline
	ThreeBodyPhaseSpace PhaseSpace() const { return _phsp; }

	__hydra_dual__ inline
	void FixParameters()
	{
		this->Parameter(0).SetFixed(true);
		this->Parameter(1).SetFixed(true);
	}

	__hydra_dual__ inline
	void ReleaseParameters()
	{
		this->Parameter(0).SetFixed(false);
		this->Parameter(1).SetFixed(false);
	}

	__hydra_dual__ inline
	hydra::complex<double> Evaluate(Arg1 mSq12, Arg1 mSq13) const
	{
		if (IsRemoved()) return hydra::complex<double>(0.,0.);

		if (_par[0]==0. && _par[1]==0.) return hydra::complex<double>(0.,0.);

		if (!_phsp.Contains<2,3>(mSq12,mSq13)) return hydra::complex<double>(0.,0.);

		return hydra::complex<double>(_par[0], _par[1]);
	}
};

} // namespace dafne
