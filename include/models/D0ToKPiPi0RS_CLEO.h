#pragma once

#include <physics/Amplitudes.h>
#include <physics/Utils.h>
#include <physics/ThreeBodyPhaseSpace.h>

#include <TColor.h>

namespace dafne {

//======================================================
// Amplitude model for D0 -> K-pi+pi0
//======================================================
// Based on CLEO's PRD 63 (2001) 092001
// [https://arxiv.org/abs/hep-ex/0011065]
//
namespace D0ToKPiPi0RS_CLEO {

const double RHO_770_P_MASS  = 0.770;
const double RHO_770_P_WIDTH = 0.1507;
const double RHO_770_P_MAG   = 1.0;
const double RHO_770_P_PHI   = 0.0;
const double RHO_770_P_CRe   = RHO_770_P_MAG*cos(RHO_770_P_PHI);
const double RHO_770_P_CIm   = RHO_770_P_MAG*sin(RHO_770_P_PHI);

const double KST_892_M_MASS  = 0.8915;
const double KST_892_M_WIDTH = 0.0500;
const double KST_892_M_MAG   = 0.44;
const double KST_892_M_PHI   = 163.*Constants::DegToRad;
const double KST_892_M_CRe   = KST_892_M_MAG*cos(KST_892_M_PHI);
const double KST_892_M_CIm   = KST_892_M_MAG*sin(KST_892_M_PHI);

const double KST_892_0_MASS  = 0.8961;
const double KST_892_0_WIDTH = 0.0505;
const double KST_892_0_MAG   = 0.39;
const double KST_892_0_PHI   = -0.2*Constants::DegToRad;
const double KST_892_0_CRe   = KST_892_0_MAG*cos(KST_892_0_PHI);
const double KST_892_0_CIm   = KST_892_0_MAG*sin(KST_892_0_PHI);

const double K0_1430_M_MASS  = 1.412;
const double K0_1430_M_WIDTH = 0.294;
const double K0_1430_M_MAG   = 0.77;
const double K0_1430_M_PHI   = 55.5*Constants::DegToRad;
const double K0_1430_M_CRe   = K0_1430_M_MAG*cos(K0_1430_M_PHI);
const double K0_1430_M_CIm   = K0_1430_M_MAG*sin(K0_1430_M_PHI);

const double K0_1430_0_MASS  = 1.412;
const double K0_1430_0_WIDTH = 0.294;
const double K0_1430_0_MAG   = 0.85;
const double K0_1430_0_PHI   = 166.*Constants::DegToRad;
const double K0_1430_0_CRe   = K0_1430_0_MAG*cos(K0_1430_0_PHI);
const double K0_1430_0_CIm   = K0_1430_0_MAG*sin(K0_1430_0_PHI);

const double RHO_1700_P_MASS  = 1.700;
const double RHO_1700_P_WIDTH = 0.240;
const double RHO_1700_P_MAG   = 2.5;
const double RHO_1700_P_PHI   = 171.*Constants::DegToRad;
const double RHO_1700_P_CRe   = RHO_1700_P_MAG*cos(RHO_1700_P_PHI);
const double RHO_1700_P_CIm   = RHO_1700_P_MAG*sin(RHO_1700_P_PHI);

const double KST_1680_M_MASS  = 1.717;
const double KST_1680_M_WIDTH = 0.322;
const double KST_1680_M_MAG   = 2.5;
const double KST_1680_M_PHI   = 103.*Constants::DegToRad;
const double KST_1680_M_CRe   = KST_1680_M_MAG*cos(KST_1680_M_PHI);
const double KST_1680_M_CIm   = KST_1680_M_MAG*sin(KST_1680_M_PHI);

const double NR_MAG = 1.75;
const double NR_PHI = 31.2*Constants::DegToRad;
const double NR_CRe = NR_MAG*cos(NR_PHI);
const double NR_CIm = NR_MAG*sin(NR_PHI);

auto PhaseSpace()
{
	return ThreeBodyPhaseSpace( Mass::D0, {Mass::K, Mass::Pi, Mass::Pi0} );
}

template<typename MSq12, typename MSq13, bool J=false, bool H=true>
__hydra_dual__ inline
auto Amplitude(ThreeBodyPhaseSpace const &phsp)
{
	if (phsp.M<1>()!=Mass::K || phsp.M<2>()!=Mass::Pi || phsp.M<3>()!=Mass::Pi0)
		HYDRA_EXCEPTION("The particles order must be K- pi+ pi0");
	
	auto radiusReso   = hydra::Parameter::Create("radiusReso").Value(1.5).Error(0.1).Fixed();
	auto radiusMother = hydra::Parameter::Create("radiusMother").Value(5.0).Error(0.1).Fixed();

	//---------------------------------------------------------------------------------------
	// Amplitudes
	//---------------------------------------------------------------------------------------
	//rho(770)+ -> pi+pi0
	auto rho_770_p_m  = hydra::Parameter::Create("RHO_770_P_m").Value(RHO_770_P_MASS ).Error(0.001);
	auto rho_770_p_w  = hydra::Parameter::Create("RHO_770_P_w").Value(RHO_770_P_WIDTH).Error(0.001);
	auto rho_770_p_re = hydra::Parameter::Create("RHO_770_P_cRe").Value(RHO_770_P_CRe  ).Error(0.001).Fixed();
	auto rho_770_p_im = hydra::Parameter::Create("RHO_770_P_cIm").Value(RHO_770_P_CIm  ).Error(0.001).Fixed();

	auto RHO_770_P_Resonance  = BreitWignerAmplitude<MSq12,MSq13,hydra::PWave,2,3,J,H>(rho_770_p_re, rho_770_p_im, rho_770_p_m, rho_770_p_w, radiusReso, radiusMother, phsp, "RHO_770_P","#it{K}^{#minus}#it{#rho}(770)^{+}");
	RHO_770_P_Resonance.SetColor(kYellow+1);
	
	//K*(892)- -> K-pi0
	auto kst_892_m_m  = hydra::Parameter::Create("KST_892_M_m").Value(KST_892_M_MASS ).Error(0.001);
	auto kst_892_m_w  = hydra::Parameter::Create("KST_892_M_w").Value(KST_892_M_WIDTH).Error(0.001);
	auto kst_892_m_re = hydra::Parameter::Create("KST_892_M_cRe").Value(KST_892_M_CRe  ).Error(0.001);
	auto kst_892_m_im = hydra::Parameter::Create("KST_892_M_cIm").Value(KST_892_M_CIm  ).Error(0.001);

	auto KST_892_M_Resonance  = BreitWignerAmplitude<MSq12,MSq13,hydra::PWave,1,3,J,H>(kst_892_m_re, kst_892_m_im, kst_892_m_m, kst_892_m_w, radiusReso, radiusMother, phsp, "KST_892_M","#it{K}*(892)^{#minus}#it{#pi}^{+}");
	KST_892_M_Resonance.SetColor(kAzure-3);
	
	//K*(892)0 -> K-pi+
	auto kst_892_0_m  = hydra::Parameter::Create("KST_892_0_m").Value(KST_892_0_MASS ).Error(0.001);
	auto kst_892_0_w  = hydra::Parameter::Create("KST_892_0_w").Value(KST_892_0_WIDTH).Error(0.001);
	auto kst_892_0_re = hydra::Parameter::Create("KST_892_0_cRe").Value(KST_892_0_CRe  ).Error(0.001);
	auto kst_892_0_im = hydra::Parameter::Create("KST_892_0_cIm").Value(KST_892_0_CIm  ).Error(0.001);

	auto KST_892_0_Resonance  = BreitWignerAmplitude<MSq12,MSq13,hydra::PWave,1,2,J,H>(kst_892_0_re, kst_892_0_im, kst_892_0_m, kst_892_0_w, radiusReso, radiusMother, phsp, "KST_892_0","#it{K}*(892)^{0}#it{#pi}^{0}");
	KST_892_M_Resonance.SetColor(kGreen+2);
		
	//K*_0(1430)- -> K-pi0
	auto k0_1430_m_m  = hydra::Parameter::Create("K0_1430_M_m").Value(K0_1430_M_MASS ).Error(0.001);
	auto k0_1430_m_w  = hydra::Parameter::Create("K0_1430_M_w").Value(K0_1430_M_WIDTH).Error(0.001);
	auto k0_1430_m_re = hydra::Parameter::Create("K0_1430_M_cRe").Value(K0_1430_M_CRe  ).Error(0.001);
	auto k0_1430_m_im = hydra::Parameter::Create("K0_1430_M_cIm").Value(K0_1430_M_CIm  ).Error(0.001);

	auto K0_1430_M_Resonance  = BreitWignerAmplitude<MSq12,MSq13,hydra::SWave,1,3,J,H>(k0_1430_m_re, k0_1430_m_im, k0_1430_m_m, k0_1430_m_w, radiusReso, radiusMother, phsp, "K0_1430_M","#it{K}*_{0}(1430)^{#minus}#it{#pi}^{+}");
	K0_1430_M_Resonance.SetColor(kMagenta);
	
	//K*_0(1430)0 -> K-pi+
	auto k0_1430_0_m  = hydra::Parameter::Create("K0_1430_0_m").Value(K0_1430_0_MASS ).Error(0.001);
	auto k0_1430_0_w  = hydra::Parameter::Create("K0_1430_0_w").Value(K0_1430_0_WIDTH).Error(0.001);
	auto k0_1430_0_re = hydra::Parameter::Create("K0_1430_0_cRe").Value(K0_1430_0_CRe  ).Error(0.001);
	auto k0_1430_0_im = hydra::Parameter::Create("K0_1430_0_cIm").Value(K0_1430_0_CIm  ).Error(0.001);

	auto K0_1430_0_Resonance  = BreitWignerAmplitude<MSq12,MSq13,hydra::SWave,1,2,J,H>(k0_1430_0_re, k0_1430_0_im, k0_1430_0_m, k0_1430_0_w, radiusReso, radiusMother, phsp, "K0_1430_0","#it{K}*_{0}(1430)^{0}#it{#pi}^{0}");
	K0_1430_0_Resonance.SetColor(kBlue);
	
	//rho(1700)+ -> pi+pi0
	auto rho_1700_p_m  = hydra::Parameter::Create("RHO_1700_P_m").Value(RHO_1700_P_MASS ).Error(0.001);
	auto rho_1700_p_w  = hydra::Parameter::Create("RHO_1700_P_w").Value(RHO_1700_P_WIDTH).Error(0.001);
	auto rho_1700_p_re = hydra::Parameter::Create("RHO_1700_P_cRe").Value(RHO_1700_P_CRe  ).Error(0.001);
	auto rho_1700_p_im = hydra::Parameter::Create("RHO_1700_P_cIm").Value(RHO_1700_P_CIm  ).Error(0.001);

	auto RHO_1700_P_Resonance = BreitWignerAmplitude<MSq12,MSq13,hydra::PWave,2,3,J,H>(rho_1700_p_re, rho_1700_p_im, rho_1700_p_m, rho_1700_p_w, radiusReso, radiusMother, phsp, "RHO_1700_P","#it{K}^{#minus}#it{#rho}(1700)^{+}");
	RHO_1700_P_Resonance.SetColor(kOrange-3);
	
	//K*(1680)- -> K-pi0
	auto kst_1680_m_m  = hydra::Parameter::Create("KST_1680_M_m").Value(KST_1680_M_MASS ).Error(0.001);
	auto kst_1680_m_w  = hydra::Parameter::Create("KST_1680_M_w").Value(KST_1680_M_WIDTH).Error(0.001);
	auto kst_1680_m_re = hydra::Parameter::Create("KST_1680_M_cRe").Value(KST_1680_M_CRe  ).Error(0.001);
	auto kst_1680_m_im = hydra::Parameter::Create("KST_1680_M_cIm").Value(KST_1680_M_CIm  ).Error(0.001);

	auto KST_1680_M_Resonance = BreitWignerAmplitude<MSq12,MSq13,hydra::PWave,1,3,J,H>(kst_1680_m_re, kst_1680_m_im, kst_1680_m_m, kst_1680_m_w, radiusReso, radiusMother, phsp, "KST_1680_M","#it{K}*(1680)^{#minus}#it{#pi}^{+}");
	KST_1680_M_Resonance.SetColor(kGray+1);
	
	//Nonresonant
	auto nr_re = hydra::Parameter::Create("NR_cRe").Value(NR_CRe).Error(0.001);
	auto nr_im = hydra::Parameter::Create("NR_cIm").Value(NR_CIm).Error(0.001);

	auto NR = NonResonantAmplitude<MSq12,MSq13>(nr_re, nr_im, phsp, "Nonresonant");
	NR.SetStyle(7);

	//---------------------------------------------------------------------------------------
	// Total amplitude
	//---------------------------------------------------------------------------------------
	return hydra::sum(RHO_770_P_Resonance, KST_892_M_Resonance, KST_892_0_Resonance, K0_1430_M_Resonance, K0_1430_0_Resonance, RHO_1700_P_Resonance, KST_1680_M_Resonance, NR);
};

}// namespace D0ToKPiPi0RS_CLEO

}// namespace dafne
