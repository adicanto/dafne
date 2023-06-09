#pragma once

#include <hydra/Function.h>
#include <hydra/functions/Utils.h>
#include <hydra/functions/BlattWeisskopfFunctions.h>

#include <physics/Utils.h>
#include <physics/ThreeBodyPhaseSpace.h>

namespace dafne {

//========================================================================================
// Functor that describes a resonance R in a 3-body decay:
//   M -> R ( -> A B) C
//
// The formalism follows CLEO (e.g., https://arxiv.org/abs/hep-ex/0011065) or,
// if BaBarKSpipi==true, BaBar's D0->KSpipi analysis (https://inspirehep.net/literature/925859)
//
// If Helicity==true the resonance pole mass is used at the denominator of the Zemach
// function in place of the two-body running mass
//
template<hydra::Wave L, unsigned A, unsigned B, bool BaBarKSpipi=false, bool Helicity=false, unsigned C=6-A-B, typename Signature=void(void)>
class Resonance: public hydra::BaseFunctor<Resonance<L,A,B,BaBarKSpipi,Helicity>, Signature, 4>
{
	typedef hydra::BaseFunctor<Resonance<L,A,B,BaBarKSpipi,Helicity>, Signature, 4> super_type;
	using super_type::_par;

private:
	ThreeBodyPhaseSpace _phsp;
	
	__hydra_dual__ inline
	double mSq23(double const &mSq12, double const &mSq13) const { return _phsp.MSqJK(mSq12,mSq13); }

public:
	Resonance() = delete;

	__hydra_dual__
	Resonance(hydra::Parameter const &mass, hydra::Parameter const &width, hydra::Parameter const &radius, hydra::Parameter const &radius_mother, ThreeBodyPhaseSpace const &phsp) : super_type({mass,width,radius,radius_mother}), _phsp(phsp)
	{
		HYDRA_STATIC_ASSERT(L<=hydra::DWave, "Resonance spin cannot exceed 2");
		HYDRA_STATIC_ASSERT(A>0 && A<4, "Invalid index of daughter A.");
		HYDRA_STATIC_ASSERT(B>0 && B<4, "Invalid index of daughter B.");
		HYDRA_STATIC_ASSERT(A!=B, "Indices of daughters cannot coincide.");
	}
	
	__hydra_dual__
	Resonance(Resonance<L,A,B,BaBarKSpipi,Helicity> const& other) : super_type(other), _phsp(other.PhaseSpace())
	{}
	
	__hydra_dual__ inline
	Resonance<L,A,B,BaBarKSpipi,Helicity>& operator=(Resonance<L,A,B,BaBarKSpipi,Helicity> const& other)
	{
		if(this==&other) return *this;
		super_type::operator=(other);
		_phsp = other.PhaseSpace();
		return *this;
	}
	
	__hydra_dual__ inline
	ThreeBodyPhaseSpace PhaseSpace() const { return _phsp; }
	
	__hydra_dual__ inline
	double Mass() const { return _par[0]; }

	__hydra_dual__ inline
	double Width() const { return _par[1]; }

	__hydra_dual__ inline
	double Radius() const { return _par[2]; }

	__hydra_dual__ inline
	double RadiusMother() const { return _par[3]; }

	__hydra_dual__ inline
	double RunningWidth(const double& mSqAB) const
	{
		double fR = BarrierFactor(mSqAB);
		return Width() * pow( Q(mSqAB)/Q0(), 2*L+1 ) * ( Mass()/sqrt(mSqAB) ) * fR * fR;
	}
	
	__hydra_dual__ inline
	hydra::Wave Spin() const { return L; }

	__hydra_dual__ inline
	double BarrierFactor(double const &mSqAB) const
	{
		return hydra::BarrierFactor<L>( Radius(), Q0(), Q(mSqAB) );
	}

	__hydra_dual__ inline
	double BarrierFactorMother(double const &mSqAB) const
	{
		return hydra::BarrierFactor<L>( RadiusMother(), P0(), P(mSqAB) );
	}

	__hydra_dual__ inline
	double Centrifugal(double const &mSqAB)
	{
		return (BaBarKSpipi) ? BarrierFactor(mSqAB) : BarrierFactor(mSqAB)*BarrierFactorMother(mSqAB);
	}
	
	__hydra_dual__ inline
	double Zemach(double const &mSqAB, double const& mSqAC, double const& mSqBC)
	{
		if (L == hydra::SWave) return 1.;
	
		auto diffSqMC = MSqMother() - MSqC();
		auto diffSqAB = MSqA() - MSqB();
		 
		auto mSq = ( Helicity ) ? Mass()*Mass() : mSqAB;
		 
		double zemach1 = mSqAC - mSqBC - diffSqMC*diffSqAB/mSq;
		 
		if (L == hydra::PWave) return zemach1;
		
		auto sumSqMC = MSqMother() + MSqC();
		auto sumSqAB = MSqA() + MSqB();

		double first  = mSqAB - 2.*sumSqMC + diffSqMC*diffSqMC/mSq;
		double second = mSqAB - 2.*sumSqAB + diffSqAB*diffSqAB/mSq;
		 
		return (zemach1*zemach1 - first*second/3.);
	}
	
	__hydra_dual__ inline
	double MSqA() const { return _phsp.MSq<A>(); }

	__hydra_dual__ inline
	double MSqB() const { return _phsp.MSq<B>(); }

	__hydra_dual__ inline
	double MSqC() const { return _phsp.MSq<C>(); }

	__hydra_dual__ inline
	double MSqMother() const { return _phsp.MSqMother(); }
	
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
   
	/*
	__hydra_dual__ inline
	double CosHelAB(const double &mSq12, const double &mSq13) const
	{
		 
	}
	*/
		
	__hydra_dual__ inline
	double P(const double &mSqAB) const
	{
		return (BaBarKSpipi) ? p_mother_frame( mSqAB, MSqMother(), MSqC() ) : p_mother_frame( MSqMother(), mSqAB, MSqC() );
	}
	
	__hydra_dual__ inline
	double Q(const double &mSqAB) const
	{
		return p_mother_frame( mSqAB, MSqA(), MSqB() );
	}

	__hydra_dual__ inline
	double Rho(const double &mSqAB) const
	{
		return sqrt( kallen( mSqAB, MSqA(), MSqB() ) ) / mSqAB;
	}

	__hydra_dual__ inline
	double P0() const
	{
		double mSqR = Mass()*Mass();
		return P(mSqR);
	}
	
	__hydra_dual__ inline
	double Q0() const
	{
		double mSqR = Mass()*Mass();
		return Q(mSqR);
	}

	__hydra_dual__ inline
	double Rho0() const
	{
		double mSqR = Mass()*Mass();
		return Rho(mSqR);
	}
	
	__hydra_dual__ inline
	bool InPhaseSpace(const double &mSq12, const double &mSq13) const
	{
		return _phsp.Contains<2,3>(mSq12,mSq13);
	}

	__hydra_dual__ inline
	void Evaluate(void ) const {}

};

} // namespace dafne
