#pragma once

#include <hydra/Function.h>
#include <hydra/Complex.h>

namespace dafne {

//========================================================================================
// GLASS amplitude as defined in BaBar:
//  - https://www.slac.stanford.edu/~wmd/kpi_swave/kpi_swave_fit.note
//  - https://inspirehep.net/literature/925859
//
template<typename R, size_t N=10, typename Signature=hydra::complex<double>(double)>
class GLassPropagator: public hydra::BaseFunctor<GLassPropagator<R>, Signature, N>
{
	typedef hydra::BaseFunctor<GLassPropagator<R>, Signature, N> super_type;
	using super_type::_par;

private:
	mutable R _resonance;

	__hydra_dual__ inline
	double lassR() const { return _par[4]; }

	__hydra_dual__ inline
	double lassB() const { return _par[5]; }

	__hydra_dual__ inline
	double phiR()  const { return _par[6]; }

	__hydra_dual__ inline
	double phiB()  const { return _par[7]; }

	__hydra_dual__ inline
	double lassr() const { return _par[8]; }

	__hydra_dual__ inline
	double lassa() const { return _par[9]; }

public:
	GLassPropagator() = delete;

	__hydra_dual__
	GLassPropagator(hydra::Parameter const &lassR, hydra::Parameter const &lassB, hydra::Parameter const &phiR, hydra::Parameter const &phiB, hydra::Parameter const &lassr, hydra::Parameter const &lassa, R const &resonance) : super_type({resonance.GetParameter(0), resonance.GetParameter(1), resonance.GetParameter(2), resonance.GetParameter(3), lassR, lassB, phiR, phiB, lassr, lassa}), _resonance(resonance)
	{
		HYDRA_STATIC_ASSERT( N==10, "Number of parameters must be 8.");
		if (_resonance.Spin()!=hydra::SWave)
			HYDRA_EXCEPTION("GLass propagator must be used only for spin-0 resonances");
	}

	__hydra_dual__ inline
	GLassPropagator<R>& operator=(GLassPropagator<R> const &other)
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
		double q = _resonance.Q(mSqAB);
		double rho0 = _resonance.Rho0();

		double qCotDeltaB = 1./lassa() + lassr()*q*q/2.;

		hydra::complex<double> rTerm = lassR() * exp( i*( phiR()+2.*phiB() ) );
		rTerm *= ( qCotDeltaB + i*q ) / ( qCotDeltaB - i*q );
		rTerm *= m*w / rho0 / (m*m - mSqAB - i*m*_resonance.RunningWidth(mSqAB) );

		hydra::complex<double> bTerm = lassB() * sqrt(mSqAB)/2. * exp(i*phiB());
		bTerm *= ( cos(phiB()) + sin(phiB()) * qCotDeltaB/q ) / ( qCotDeltaB - i*q);

		return rTerm + bTerm;
	}
};

} // namespace dafne

