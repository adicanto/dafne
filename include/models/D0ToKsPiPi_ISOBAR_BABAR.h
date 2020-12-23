#include <physics/ThreeBodyPhaseSpace.h>
#include <physics/Amplitudes.h>
#include <physics/Utils.h>

#include <TColor.h>

namespace dafne {

//======================================================
// Amplitude model for D0 -> Ks pi+ pi-
//======================================================
// Based on BaBar's PRL 105, 081803 (2010)
// [http://arxiv.org/abs/1004.5053] 

namespace D0ToKsPiPi_ISOBAR_BABAR {

auto PhaseSpace()
{
	return ThreeBodyPhaseSpace( 1.8645, {0.49767, 0.139570, 0.139570} );
}

template<typename MSq12, typename MSq13, bool Jordi=true, bool Helicity=false>
__hydra_dual__ inline
auto Amplitude(ThreeBodyPhaseSpace const &phsp)
{
	if (phsp.M<2>()!=phsp.M<3>())
		HYDRA_EXCEPTION("The particles order must be K0_S pi+ pi-");
	
	auto radius = hydra::Parameter::Create("rBW").Value(1.5).Error(0.1).Fixed();

	// Non resonant
	auto c_re  = hydra::Parameter::Create("reCtnt").Value(0.848984).Error(0.1);
	auto c_im  = hydra::Parameter::Create("imCtnt").Value(0.893618).Error(0.1);
	auto NR_amp = NonResonantAmplitude<MSq12,MSq13>(c_re,c_im,phsp, "noRes");
	NR_amp.SetStyle(2);
	
	// K*(892)- -> KS Pi-
	c_re  = hydra::Parameter::Create("reKstm").Value(-1.16356).Error(0.1);
	c_im  = hydra::Parameter::Create("imKstm").Value(1.19933).Error(0.1);
	auto mass  = hydra::Parameter::Create("mKst").Value(0.8936060).Error(0.1);
	auto width = hydra::Parameter::Create("wKst").Value(0.0463407).Error(0.1);
	
	auto KSTM_amp = BreitWignerAmplitude<MSq12,MSq13,hydra::PWave,1,3,Jordi,Helicity>(c_re,c_im,mass,width,radius,radius,phsp,"Kstm");
	KSTM_amp.SetColor(kYellow+1);

	// K*(892)+ -> KS pi+
	c_re  = hydra::Parameter::Create("reKstp").Value(0.106051).Error(0.1);
	c_im  = hydra::Parameter::Create("imKstp").Value(-0.118513).Error(0.1);
	auto KSTP_amp = BreitWignerAmplitude<MSq12,MSq13,hydra::PWave,1,2,Jordi,Helicity>(c_re,c_im,mass,width,radius,radius,phsp,"Kstp");
	KSTP_amp.SetColor(kAzure-3);

	// rho(770)0 -> pi+pi-
	c_re  = hydra::Parameter::Create("rerho0").Value(1.).Error(0.1);
	c_im  = hydra::Parameter::Create("imrho0").Value(0.).Error(0.1);
	mass  = hydra::Parameter::Create("mRho").Value(0.7758).Error(0.1);
	width = hydra::Parameter::Create("wRho").Value(0.1464).Error(0.1);
	auto RHO_amp = GounarisSakuraiAmplitude<MSq12,MSq13,hydra::PWave,2,3,Jordi,Helicity>(c_re,c_im,mass,width,radius,radius,phsp,"rho0");
	RHO_amp.SetColor(kGreen+2);
	
	// omega(782) -> pi+pi-
	c_re  = hydra::Parameter::Create("reomega").Value(-0.0249569).Error(0.1);
	c_im  = hydra::Parameter::Create("imomega").Value(0.0388072).Error(0.1);
	mass  = hydra::Parameter::Create("momega").Value(0.78259).Error(0.1);
	width = hydra::Parameter::Create("womega").Value(0.00849).Error(0.1);

	auto OMEGA_amp = BreitWignerAmplitude<MSq12,MSq13,hydra::PWave,2,3,Jordi,Helicity>(c_re,c_im,mass,width,radius,radius,phsp,"omega");
	OMEGA_amp.SetColor(kBlue);
	
	// f0(980) -> pi+pi-
	c_re  = hydra::Parameter::Create("ref0_980").Value(-0.423586).Error(0.1);
	c_im  = hydra::Parameter::Create("imf0_980").Value(-0.236099).Error(0.1);
	mass  = hydra::Parameter::Create("mf0_980").Value(0.975).Error(0.1);
	width = hydra::Parameter::Create("wf0_980").Value(0.044).Error(0.1);	

	auto F0_980_amp = BreitWignerAmplitude<MSq12,MSq13,hydra::SWave,2,3,Jordi,Helicity>(c_re,c_im,mass,width,radius,radius,phsp,"f0_980");
	F0_980_amp.SetColor(kYellow+1);
	F0_980_amp.SetStyle(2);
	
	// f0(1370) -> pi+pi-
	c_re  = hydra::Parameter::Create("ref0_1370").Value(-2.16486).Error(0.1);
	c_im  = hydra::Parameter::Create("imf0_1370").Value(3.62385).Error(0.1);
	mass  = hydra::Parameter::Create("mf0_1370").Value(1.434).Error(0.1);
	width = hydra::Parameter::Create("wf0_1370").Value(0.173).Error(0.1);	

	auto F0_1370_amp = BreitWignerAmplitude<MSq12,MSq13,hydra::SWave,2,3,Jordi,Helicity>(c_re,c_im,mass,width,radius,radius,phsp,"f0_1370");
	F0_1370_amp.SetColor(kMagenta);
	F0_1370_amp.SetStyle(2);

	// f2(1270) -> pi+pi-
	c_re  = hydra::Parameter::Create("ref2_1270").Value(0.217748).Error(0.1);
	c_im  = hydra::Parameter::Create("imf2_1270").Value(-0.133327).Error(0.1);
	mass  = hydra::Parameter::Create("mf2_1270").Value(1.2754).Error(0.1);
	width = hydra::Parameter::Create("wf2_1270").Value(0.1851).Error(0.1);	

	auto F2_1270_amp = BreitWignerAmplitude<MSq12,MSq13,hydra::DWave,2,3,Jordi,Helicity>(c_re,c_im,mass,width,radius,radius,phsp,"f2_1270");
	F2_1270_amp.SetColor(kMagenta);
	
	// K0*(1430)- -> Ks pi-
	c_re  = hydra::Parameter::Create("reK0stm_1430").Value(1.62128).Error(0.1);
	c_im  = hydra::Parameter::Create("imK0stm_1430").Value(1.06816).Error(0.1);
	mass  = hydra::Parameter::Create("mK0st_1430").Value(1.459).Error(0.1);
	width = hydra::Parameter::Create("wK0st_1430").Value(0.175).Error(0.1);	

	auto K0STM_1430_amp = BreitWignerAmplitude<MSq12,MSq13,hydra::SWave,1,3,Jordi,Helicity>(c_re,c_im,mass,width,radius,radius,phsp,"K0stm_1430");
	K0STM_1430_amp.SetColor(kOrange-3);
	
	// K0*(1430)+ -> Ks pi+
	c_re  = hydra::Parameter::Create("reK0stp_1430").Value(0.148802).Error(0.1);
	c_im  = hydra::Parameter::Create("imK0stp_1430").Value(-0.118513).Error(0.1);

	auto K0STP_1430_amp = BreitWignerAmplitude<MSq12,MSq13,hydra::SWave,1,2,Jordi,Helicity>(c_re,c_im,mass,width,radius,radius,phsp,"K0stp_1430");
	K0STP_1430_amp.SetColor(kGray+1);

	// K2*(1430)- -> Ks pi-
	c_re  = hydra::Parameter::Create("reK2stm_1430").Value(1.15489).Error(0.1);
	c_im  = hydra::Parameter::Create("imK2stm_1430").Value(-0.773363).Error(0.1);
	mass  = hydra::Parameter::Create("mK2st_1430").Value(1.4256).Error(0.1);
	width = hydra::Parameter::Create("wK2st_1430").Value(0.0985).Error(0.1);	

	auto K2STM_1430_amp = BreitWignerAmplitude<MSq12,MSq13,hydra::DWave,1,3,Jordi,Helicity>(c_re,c_im,mass,width,radius,radius,phsp,"K2stm_1430");
	K2STM_1430_amp.SetColor(kAzure-3);
	K2STM_1430_amp.SetStyle(2);
	
	// K2*(1430)+ -> Ks pi+
	c_re  = hydra::Parameter::Create("reK2stp_1430").Value(0.140865).Error(0.1);
	c_im  = hydra::Parameter::Create("imK2stp_1430").Value(-0.165378).Error(0.1);

	auto K2STP_1430_amp = BreitWignerAmplitude<MSq12,MSq13,hydra::DWave,1,2,Jordi,Helicity>(c_re,c_im,mass,width,radius,radius,phsp,"K2stp_1430");
	K2STP_1430_amp.SetColor(kGreen+2);
	K2STP_1430_amp.SetStyle(2);
	
	// sigma -> pi+pi-
	c_re  = hydra::Parameter::Create("resigma").Value(-1.55556).Error(0.1);
	c_im  = hydra::Parameter::Create("imsigma").Value(-0.931685).Error(0.1);
	mass  = hydra::Parameter::Create("msigma").Value(0.527699).Error(0.1);
	width = hydra::Parameter::Create("wsigma").Value(0.511861).Error(0.1);
	auto SIGMA_amp = GounarisSakuraiAmplitude<MSq12,MSq13,hydra::SWave,2,3,Jordi,Helicity>(c_re,c_im,mass,width,radius,radius,phsp,"sigma");
	SIGMA_amp.SetColor(kOrange-3);
	SIGMA_amp.SetStyle(2);

	// sigma2 -> pi+pi-
	c_re  = hydra::Parameter::Create("resigma2").Value(-0.273791).Error(0.1);
	c_im  = hydra::Parameter::Create("imsigma2").Value(-0.0535596).Error(0.1);
	mass  = hydra::Parameter::Create("msigma2").Value(1.03327).Error(0.1);
	width = hydra::Parameter::Create("wsigma2").Value(0.0987890).Error(0.1);
	auto SIGMA2_amp = GounarisSakuraiAmplitude<MSq12,MSq13,hydra::SWave,2,3,Jordi,Helicity>(c_re,c_im,mass,width,radius,radius,phsp,"sigma2");
	SIGMA2_amp.SetColor(kBlue);
	SIGMA2_amp.SetStyle(3);

	// K*(1680)- -> Ks pi-
	c_re  = hydra::Parameter::Create("reKstm_1680").Value(1.62128).Error(0.1);
	c_im  = hydra::Parameter::Create("imKstm_1680").Value(1.06816).Error(0.1);
	mass  = hydra::Parameter::Create("mKstm_1680").Value(1.459).Error(0.1);
	width = hydra::Parameter::Create("wKstm_1680").Value(0.175).Error(0.1);	

	auto KSTM_1680_amp = BreitWignerAmplitude<MSq12,MSq13,hydra::PWave,1,3,Jordi,Helicity>(c_re,c_im,mass,width,radius,radius,phsp,"Kstm_1680");
	KSTM_1680_amp.SetColor(kBlue);
	KSTM_1680_amp.SetStyle(2);
	
	return hydra::sum(NR_amp,KSTM_amp,KSTP_amp,RHO_amp,OMEGA_amp,F0_980_amp,F0_1370_amp,F2_1270_amp,K0STM_1430_amp,K0STP_1430_amp,K2STM_1430_amp,K2STP_1430_amp,SIGMA_amp,SIGMA2_amp,KSTM_1680_amp);
}

}// namespace D0ToKsPiPi_ISOBAR_BABAR

}// namespace dafne
