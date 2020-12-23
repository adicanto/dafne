#pragma once

#include <hydra/Function.h>
#include <hydra/Complex.h>
#include <hydra/functions/Utils.h>
#include <hydra/functions/BlattWeisskopfFunctions.h>

#include <physics/Utils.h>
#include <physics/BaseAmplitude.h>
#include <physics/ThreeBodyPhaseSpace.h>
#include <physics/PDG.h>

#include <Eigen/Dense>

namespace dafne {

template<typename Arg1, typename Arg2, unsigned A, unsigned B, unsigned C=6-A-B, typename Signature=hydra::complex<double>(Arg1,Arg2)>
class Fvector: public BaseAmplitude, public hydra::BaseFunctor<Fvector<Arg1,Arg2,A,B,C>, Signature, 21>
{
	typedef hydra::BaseFunctor<Fvector<Arg1,Arg2,A,B,C>, Signature, 21> super_type;
	using super_type::_par;

private:
	ThreeBodyPhaseSpace _phsp;
	bool _usePvecSvp;// Decide if the slowly varying part of the P-vector should be used.
	
	Eigen::Matrix<double, 5, 1> _m0  = Eigen::Matrix<double, 5, 1>::Zero();
	Eigen::Matrix<double, 5, 5> _g0  = Eigen::Matrix<double, 5, 5>::Zero();
	Eigen::Matrix<double, 5, 5> _fSc = Eigen::Matrix<double, 5, 5>::Zero();

	double _s0sc = -3.92637;
	double _s0A  = -0.15;
	double _sA   = 1.0;

	// in principle the m0, _g0, fSc above could be static members, however, I am not sure how would device (like GPU) treat the static members. If you
	// are sure the statics members could be handled correctly in the device, you could use the following static members while removing the InitializeConstantMatrices()
	// in copy constructor and operator=, uncommenting the Eigen::Matrix initialization code in the end of this file
	// static bool hadInitializedConstantMatrices;
	// static Eigen::Matrix<double, 5, 1> _m0;
	// static Eigen::Matrix<double, 5, 5> _g0;
	// static Eigen::Matrix<double, 5, 5> _fSc;
	// static double _s0sc;
	// static double _s0A; 
	// static double _sA;

	__hydra_dual__ inline
	double mSq23(double const &mSq12, double const &mSq13) const { return _phsp.MSqJK(mSq12,mSq13); }

public:
	Fvector() = delete;
	
	__hydra_dual__
	Fvector(std::vector<hydra::Parameter> const& beta, std::vector<hydra::Parameter> const& fPr, hydra::Parameter const& s0pr, ThreeBodyPhaseSpace const &phsp, const char* name="", const char *label="") : BaseAmplitude(name,label), super_type(), _phsp(phsp), _usePvecSvp(true)
	{
		// judge the length of beta[] and fPr[], for the moment only 5*5 matrix is allowed
		if (beta.size() != 5*2 || fPr.size() != 5*2) { // real and image part
			HYDRA_EXCEPTION("For the moment, only 5x5 matrix is allowed!");
		}

		// map production parameters (beta) and production background parameters (fPr and s0pr) to _par[]
		for (size_t i = 0; i < 10; ++i) {
			this->SetParameter(i, beta[i]);
		}
		for (size_t i = 0; i < 10; ++i) {
			this->SetParameter(10+i, fPr[i]);
		}
		this->SetParameter(20, s0pr);

		// initialize the fixed KMatrix
		InitializeConstantMatrices();
	}

	__hydra_dual__
	Fvector(Fvector<Arg1,Arg2,A,B,C> const& other) : BaseAmplitude(other), super_type(other), _phsp(other.PhaseSpace()), _usePvecSvp(other.UsePvecSvp())
	{
		// when using m0, _g0, fSc as static members, the following line code be removed
		InitializeConstantMatrices();
	} 

	__hydra_dual__ inline
	Fvector<Arg1,Arg2,A,B,C>& operator=(Fvector<Arg1,Arg2,A,B,C> const& other)
	{
		if(this==&other) return *this;
		BaseAmplitude::operator=(other);
		super_type::operator=(other);
		_phsp = other.PhaseSpace();
		_usePvecSvp = other.UsePvecSvp();
		
		// when using m0, _g0, fSc as static members, the following line code be removed
		InitializeConstantMatrices();
		return *this;		
	}

	__hydra_dual__ inline
	void InitializeConstantMatrices()
	{
		// uncomment the following line when you use m0, _g0, fSc as static members
		// if (hadInitializedConstantMatrices == true) return;

		Eigen::Matrix<double, 5, 1> m0;
		m0  << 0.65100,
		       1.20360,
		       1.55817,
		       1.21000,
		       1.82206;
		_m0 = m0;

		Eigen::Matrix<double, 5, 5> g0;
		g0  <<  0.22889, -0.55377,  0.00000, -0.39899, -0.34639,
				0.94128,  0.55095,  0.00000,  0.39065,  0.31503,
				0.36856,  0.23888,  0.55639,  0.18340,  0.18681,
				0.33650,  0.40907,  0.85679,  0.19906, -0.00984,
				0.18171, -0.17558, -0.79658, -0.00355,  0.22358;
		_g0 = g0;

		_fSc = Eigen::Matrix<double, 5, 5>::Zero();
		_fSc(0, 0) = 0.23399,
	    _fSc(0, 1) = 0.15044,
	    _fSc(0, 2) = -0.20545,
	    _fSc(0, 3) = 0.32825,
	    _fSc(0, 4) = 0.35412;
		for (int r = 1; r < 5; ++r)
		    _fSc( r, 0 ) = _fSc( 0, r );

		// uncomment the following line when you use m0, _g0, fSc as static members
		// hadInitializedConstantMatrices = true; 
	}
	
	__hydra_dual__ inline
	bool UsePvecSvp() const { return _usePvecSvp; }

	__hydra_dual__ inline
	void FixParameters()
	{
		size_t n = this->GetNumberOfParameters();
		for (size_t i=0; i<n; ++i) this->Parameter(i).SetFixed(true);
	}

	__hydra_dual__ inline
	void ReleaseParameters()
	{
		size_t n = this->GetNumberOfParameters();
		for (size_t i=0; i<n; ++i) this->Parameter(i).SetFixed(false);
	}

	__hydra_dual__ inline
	ThreeBodyPhaseSpace PhaseSpace() const { return _phsp; }

	__hydra_dual__ inline
	hydra::complex<double> beta(size_t i) const
	{
		// judge the range of i, if needed
		return hydra::complex<double>(_par[i], _par[5+i]);
	}

	__hydra_dual__ inline
	hydra::complex<double> fPr(size_t i) const
	{
		// judge the range of i, if needed
		return hydra::complex<double>(_par[10+i], _par[15+i]);
	}

	__hydra_dual__ inline
	double s0pr() const
	{
		return _par[20];
	}

	__hydra_dual__ inline
	Eigen::Matrix<double, 5, 1> m0() { return _m0;}
	Eigen::Matrix<double, 5, 5> g0() { return _g0;}
	Eigen::Matrix<double, 5, 5> fSc() {return _fSc;};


	__hydra_dual__ inline
	double Rho(const double &mSqAB) const
	{
		return sqrt( kallen( mSqAB, MSqA(), MSqB() ) ) / mSqAB;
	}

	__hydra_dual__ inline
	hydra::complex<double> Rho( const double mCh, const double mSqAB ) const
	{
		double rhoSq = 1. - pow( mCh, 2 ) / mSqAB;

		if ( rhoSq >= 0. )
			return sqrt( rhoSq );

		const hydra::complex<double> i( 0., 1. );
		return i * sqrt( - rhoSq );
	}

	__hydra_dual__ inline
	hydra::complex<double> Rho4pi(const double mSqAB) const
	{
		const double mPi = Mass::Pi;

		if (mSqAB > 1.)
			return Rho(4.*mPi, mSqAB);

		const double m4AB = std::pow(mSqAB, 2);
		const double m6AB = std::pow(mSqAB, 3);
		const double m8AB = std::pow(mSqAB, 4);

		double term = 0.;
		term += 0.00370909 / m4AB;
		term -= 0.111203   / mSqAB;
		term += 1.2274;
		term -= 6.39017    * mSqAB;
		term += 16.8358    * m4AB;
		term -= 21.8845    * m6AB;
		term += 11.3153    * m8AB;

		return Rho(4.*mPi, 1.) * term;
	}

	__hydra_dual__ inline
	hydra::complex< double > Rho(const int index, const double mSqAB) const
	{
		const double mPi   = Mass::Pi;
		const double mK    = Mass::K; // Charged K mass.
		const double mEta  = Mass::Eta;
		const double mEtaP = Mass::EtaPrime;

		if ( index == 0 ) return Rho   (2. * mPi    , mSqAB);
		if ( index == 1 ) return Rho   (2. * mK     , mSqAB);
		if ( index == 2 ) return Rho4pi(              mSqAB);
		if ( index == 3 ) return Rho   (2. * mEta   , mSqAB);
		if ( index == 4 ) return Rho   (mEta + mEtaP, mSqAB);

		return hydra::complex<double>(1., 0.);
	}

	__hydra_dual__ inline
	double MSqA() const { return _phsp.MSq<A>(); }

	__hydra_dual__ inline
	double MSqB() const { return _phsp.MSq<B>(); }

	__hydra_dual__ inline
	double MSqC() const { return _phsp.MSq<C>(); }

	__hydra_dual__ inline
	double MSqAB(const double &mSq12, const double &mSq13) const
	{
		if ( C == 3 ) return mSq12;
		if ( C == 2 ) return mSq13;
		return mSq23(mSq12,mSq13);
	}

	__hydra_dual__ inline
	double MSqAC(const double &mSq12, const double &mSq13) const
	{
		if ( B == 3 ) return mSq12;
		if ( B == 2 ) return mSq13;
		return mSq23(mSq12,mSq13);
	}
    
	__hydra_dual__ inline
	double MSqBC(const double &mSq12, const double &mSq13) const
	{
		if ( A == 3 ) return mSq12;
		if ( A == 2 ) return mSq13;
		return mSq23(mSq12,mSq13);
	}
   
	// Momentum of a resonant particle in the rest frame of the resonant pair.
	__hydra_dual__ inline
	double Q(const double &mSqAB) const
	{
		return p_mother_frame( mSqAB, MSqA(), MSqB() );
	}

	__hydra_dual__ inline
	double  Q2(const double& mSqAB ) const
	{
 		return pow(Q(),2);
	}

	__hydra_dual__ inline
	hydra::complex<double> Evaluate(Arg1 mSq12, Arg2 mSq13) const
	{
		if (IsRemoved()) return hydra::complex<double>(0.,0.);

		if (!_phsp.Contains<2,3>(mSq12,mSq13)) return hydra::complex<double>(0.,0.);

		const double mSqAB = MSqAB(mSq12, mSq13);

		// Initialize the K matrix to zero.
		Eigen::Matrix<double, 5, 5> K = Eigen::Matrix<double, 5, 5>::Zero();

		// Resonant contribution.
		for (int row = 0; row < 5; ++row)
		for (int col = 0; col < 5; ++col)
		  for (int pole = 0; pole < 5; ++pole)
		    K(row, col) += _g0(pole, row) * _g0(pole, col) / (pow( _m0[pole], 2) - mSqAB ); // (pow( _m0[pole], 2) - mSqAB ) might have inf problem

		// Non-resonant contribution.
		K += _fSc * (1. - _s0sc) / (mSqAB - _s0sc);

		// Adler term.
		K *= (1. - _s0A) / (mSqAB - _s0A) * (mSqAB - _sA * _phsp.M<A>() * _phsp.M<B>() / 2. );

		Eigen::Matrix<std::complex<double>, 5, 5> M;
		const hydra::complex<double> i(0., 1.);

		// Build M = ( 1 - i K rho ).
		for (int row = 0; row < 5; ++row)
		for (int col = 0; col < 5; ++col)
	  		M(row, col) = 1.*(row == col) - i*K(row, col)*Rho(col, mSqAB);

		// Invert M: ( 1 - i K rho )^{ -1 }
		Eigen::Matrix<std::complex<double>, 5, 5> invM = M.inverse();

		// Decide if the slowly varying part should be used.
		double svp = 1.0;
		if (_usePvecSvp)
			svp = (1. - s0pr()) / (mSqAB - s0pr()); // Slowly varying part.

		// Build the P vector.
		std::vector<hydra::complex<double>> P;
		for (int row = 0; row < 5; ++row )
		{
			hydra::complex< double > poles = 0.;
			for ( int pole = 0; pole < 5; ++pole )
			  poles += beta(pole) * _g0(pole,row) / ( std::pow(_m0(pole), 2) - mSqAB);

			P.push_back(poles + fPr(row) * svp);
		}

		hydra::complex<double> result = 0.;

		// Compute the first component of the F vector: F_0 = ( 1 - i K rho )_{0j}^{ -1 } P_j.
		for (int line = 0; line < 5; ++line)
			result += hydra::complex<double>(invM(0, line)) * P[line];

		return result;
	}



	
};

// uncomment the following initialize lines when you use m0, _g0, fSc as static members
// template<typename Arg1, typename Arg2, unsigned A, unsigned B, unsigned C, typename Signature> bool Fvector<Arg1,Arg2,A,B,C,Signature>::hadInitializedConstantMatrices = false;
// template<typename Arg1, typename Arg2, unsigned A, unsigned B, unsigned C, typename Signature> Eigen::Matrix<double, 5, 1> Fvector<Arg1,Arg2,A,B,C,Signature>::_m0 = Eigen::Matrix<double, 5, 1>::Zero();
// template<typename Arg1, typename Arg2, unsigned A, unsigned B, unsigned C, typename Signature> Eigen::Matrix<double, 5, 5> Fvector<Arg1,Arg2,A,B,C,Signature>::_g0 = Eigen::Matrix<double, 5, 5>::Zero();
// template<typename Arg1, typename Arg2, unsigned A, unsigned B, unsigned C, typename Signature> Eigen::Matrix<double, 5, 5> Fvector<Arg1,Arg2,A,B,C,Signature>::_fSc = Eigen::Matrix<double, 5, 5>::Zero();
// template<typename Arg1, typename Arg2, unsigned A, unsigned B, unsigned C, typename Signature> double Fvector<Arg1,Arg2,A,B,C,Signature>::_s0sc = -3.92637;
// template<typename Arg1, typename Arg2, unsigned A, unsigned B, unsigned C, typename Signature> double Fvector<Arg1,Arg2,A,B,C,Signature>::_s0A = -0.15; 
// template<typename Arg1, typename Arg2, unsigned A, unsigned B, unsigned C, typename Signature> double Fvector<Arg1,Arg2,A,B,C,Signature>::_sA = 1.0;

	
} // dafne
