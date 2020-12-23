#pragma once

#include <physics/Resonance.h>
#include <physics/BreitWignerPropagator.h>
#include <physics/GounarisSakuraiPropagator.h>
#include <physics/GLassPropagator.h>
#include <physics/LassPropagator.h>
#include <physics/FlattePropagator.h>
#include <physics/Fvector.h>
#include <physics/NonResonant.h>
#include <physics/Amplitude.h>

namespace dafne {

// Relativistic Breit-Wigner
template<typename Arg1, typename Arg2, hydra::Wave L, unsigned A, unsigned B, bool Jordi=false, bool Helicity=false>
auto BreitWignerAmplitude(hydra::Parameter const &c_re, hydra::Parameter const &c_im, hydra::Parameter const &mass, hydra::Parameter const &width, hydra::Parameter const &radius, hydra::Parameter const &radius_mother, ThreeBodyPhaseSpace const &phsp, const char *name="", const char *label="")
{
	auto r = Resonance<L,A,B,Jordi,Helicity>(mass,width,radius,radius_mother,phsp);
	auto p = BreitWignerPropagator<Resonance<L,A,B,Jordi,Helicity>>(r);
	return Amplitude<Arg1,Arg2,4,BreitWignerPropagator<Resonance<L,A,B,Jordi,Helicity>>>(c_re,c_im,p,name,label);
}

// Gounaris-Sakurai
template<typename Arg1, typename Arg2, hydra::Wave L, unsigned A, unsigned B, bool Jordi=false, bool Helicity=false>
auto GounarisSakuraiAmplitude(hydra::Parameter const &c_re, hydra::Parameter const &c_im, hydra::Parameter const &mass, hydra::Parameter const &width, hydra::Parameter const &radius, hydra::Parameter const &radius_mother, ThreeBodyPhaseSpace const &phsp, const char *name="", const char *label="")
{
	auto r = Resonance<L,A,B,Jordi,Helicity>(mass,width,radius,radius_mother,phsp);
	auto p = GounarisSakuraiPropagator<Resonance<L,A,B,Jordi,Helicity>>(r);
	return Amplitude<Arg1,Arg2,4,GounarisSakuraiPropagator<Resonance<L,A,B,Jordi,Helicity>>>(c_re,c_im,p,name,label);
}

// Generalized LASS
template<typename Arg1, typename Arg2, hydra::Wave L, unsigned A, unsigned B, bool Jordi=false, bool Helicity=false>
auto GLassAmplitude(hydra::Parameter const &c_re, hydra::Parameter const &c_im, hydra::Parameter const &mass, hydra::Parameter const &width, hydra::Parameter const &radius, hydra::Parameter const &radius_mother, hydra::Parameter const &lassR, hydra::Parameter const &lassB, hydra::Parameter const &phiR, hydra::Parameter const &phiB, hydra::Parameter const &lassr, hydra::Parameter const &lassa, ThreeBodyPhaseSpace const &phsp, const char *name="", const char *label="")
{
	auto r = Resonance<L,A,B,Jordi,Helicity>(mass,width,radius,radius_mother,phsp);
	auto p = GLassPropagator<Resonance<L,A,B,Jordi,Helicity>>(lassR,lassB,phiR,phiB,lassr,lassa,r);
	return Amplitude<Arg1,Arg2,10,GLassPropagator<Resonance<L,A,B,Jordi,Helicity>>>(c_re,c_im,p,name,label);
}

// LASS
template<typename Arg1, typename Arg2, hydra::Wave L, unsigned A, unsigned B, bool Jordi=false, bool Helicity=false>
auto LassAmplitude(hydra::Parameter const &c_re, hydra::Parameter const &c_im, hydra::Parameter const &mass, hydra::Parameter const &width, hydra::Parameter const &radius, hydra::Parameter const &radius_mother, hydra::Parameter const &b1, hydra::Parameter const &b2, hydra::Parameter const &b3, hydra::Parameter const &lassr, hydra::Parameter const &lassa, ThreeBodyPhaseSpace const &phsp, const char *name="", const char *label="")
{
	auto r = Resonance<L,A,B,Jordi,Helicity>(mass,width,radius,radius_mother,phsp);
	auto p = LassPropagator<Resonance<L,A,B,Jordi,Helicity>>(b1,b2,b3,lassr,lassa,r);
	return Amplitude<Arg1,Arg2,9,LassPropagator<Resonance<L,A,B,Jordi,Helicity>>>(c_re,c_im,p,name,label);
}

// Flatte
template<typename Arg1, typename Arg2, hydra::Wave L, unsigned A, unsigned B, bool Jordi=false, bool Helicity=false>
auto FlatteAmplitude(hydra::Parameter const &c_re, hydra::Parameter const &c_im, hydra::Parameter const &mass, hydra::Parameter const &width, hydra::Parameter const &radius, hydra::Parameter const &radius_mother, hydra::Parameter const &gamma1, hydra::Parameter const &gamma2, hydra::Parameter const &m02a, hydra::Parameter const &m02b, ThreeBodyPhaseSpace const &phsp, const char *name="", const char *label="")
{
	auto r = Resonance<L,A,B,Jordi,Helicity>(mass,width,radius,radius_mother,phsp);
	auto p = FlattePropagator<Resonance<L,A,B,Jordi,Helicity>>(gamma1,gamma2,m02a,m02b,r);
	return Amplitude<Arg1,Arg2,8,FlattePropagator<Resonance<L,A,B,Jordi,Helicity>>>(c_re,c_im,p,name,label);
}

// Fvector
template<typename Arg1, typename Arg2, unsigned A, unsigned B>
auto FvectorAmplitude(std::vector<hydra::Parameter> const& beta, std::vector<hydra::Parameter> const& fPr, hydra::Parameter const& s0pr, ThreeBodyPhaseSpace const &phsp, const char* name="", const char *label="")
{
	return Fvector<Arg1,Arg2,A,B>(beta,fPr,s0pr,phsp,name,label);
}

// Non resonant
template<typename Arg1, typename Arg2>
auto NonResonantAmplitude(hydra::Parameter const &c_re, hydra::Parameter const &c_im, ThreeBodyPhaseSpace const &phsp, const char *name="", const char *label="")
{
	return NonResonant<Arg1,Arg2>(c_re,c_im,phsp,name,label);
}

} // namespace dafne
