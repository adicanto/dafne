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

namespace D0ToKsPiPi_FVECTOR_BABAR {

__hydra_dual__ inline
auto PhaseSpace()
{
	return ThreeBodyPhaseSpace( 1.8645, {0.49767, 0.139570, 0.139570} );
}

__hydra_dual__ inline
auto PhaseSpaceWithTime()
{
	return ThreeBodyPhaseSpaceWithTime( 1.8645, {0.49767, 0.139570, 0.139570}, {0.,15.*Tau::D0} );
}

__hydra_dual__ inline
auto PhaseSpaceWithTimeAndTimeError()
{
	return ThreeBodyPhaseSpaceWithTimeAndTimeError( 1.8645, {0.49767, 0.139570, 0.139570}, {-2., 7.}, {0, 0.5});
}


template<typename MSq12, typename MSq13, bool BaBarKSpipi=true, bool Helicity=false>
__hydra_dual__ inline
auto Amplitude(ThreeBodyPhaseSpace const &phsp)
{
	if (phsp.M<2>()!=phsp.M<3>())
		HYDRA_EXCEPTION("The particles order must be K0_S pi+ pi-");
	
	auto radius = hydra::Parameter::Create("rBW").Value(1.5).Error(0.1).Fixed();
	
	// K*(892)- -> KS Pi-
	auto c_re  = hydra::Parameter::Create("reKstm").Value(-1.16356).Error(0.1);
	auto c_im  = hydra::Parameter::Create("imKstm").Value(1.19933).Error(0.1);
	auto mass  = hydra::Parameter::Create("mKst").Value(0.8936060).Error(0.1);
	auto width = hydra::Parameter::Create("wKst").Value(0.0463407).Error(0.1);
	
	auto KSTM_amp = BreitWignerAmplitude<MSq12,MSq13,hydra::PWave,1,3,BaBarKSpipi,Helicity>(c_re,c_im,mass,width,radius,radius,phsp,"Kstm","#it{K}*(892)^{#minus}#it{#pi}^{+}");
	KSTM_amp.SetColor(kYellow+1);
	
	// K*(892)+ -> KS pi+
	c_re  = hydra::Parameter::Create("reKstp").Value(0.106051).Error(0.1);
	c_im  = hydra::Parameter::Create("imKstp").Value(-0.118513).Error(0.1);
	auto KSTP_amp = BreitWignerAmplitude<MSq12,MSq13,hydra::PWave,1,2,BaBarKSpipi,Helicity>(c_re,c_im,mass,width,radius,radius,phsp,"Kstp","#it{K}*(892)^{+}#it{#pi}^{#minus}");
	KSTP_amp.SetColor(kAzure-3);
	
	// rho(770)0 -> pi+pi-
	c_re  = hydra::Parameter::Create("rerho0").Value(1.).Error(0.1);
	c_im  = hydra::Parameter::Create("imrho0").Value(0.).Error(0.1);
	mass  = hydra::Parameter::Create("mRho").Value(0.7758).Error(0.1);
	width = hydra::Parameter::Create("wRho").Value(0.1464).Error(0.1);
	auto RHO_amp = GounarisSakuraiAmplitude<MSq12,MSq13,hydra::PWave,2,3,BaBarKSpipi,Helicity>(c_re,c_im,mass,width,radius,radius,phsp,"rho0","#it{K}^{0}_{S}#it{#rho}(770)^{0}");
	RHO_amp.SetColor(kGreen+2);

	// omega(782) -> pi+pi-
	c_re  = hydra::Parameter::Create("reomega").Value(-0.0249569).Error(0.1);
	c_im  = hydra::Parameter::Create("imomega").Value(0.0388072).Error(0.1);
	mass  = hydra::Parameter::Create("momega").Value(0.78259).Error(0.1);
	width = hydra::Parameter::Create("womega").Value(0.00849).Error(0.1);

	auto OMEGA_amp = BreitWignerAmplitude<MSq12,MSq13,hydra::PWave,2,3,BaBarKSpipi,Helicity>(c_re,c_im,mass,width,radius,radius,phsp,"omega","#it{K}^{0}_{S}#it{#omega}^{0}");
	OMEGA_amp.SetColor(kBlue);
	
	// f2(1270) -> pi+pi-
	c_re  = hydra::Parameter::Create("ref2_1270").Value(0.217748).Error(0.1);
	c_im  = hydra::Parameter::Create("imf2_1270").Value(-0.133327).Error(0.1);
	mass  = hydra::Parameter::Create("mf2_1270").Value(1.2754).Error(0.1);
	width = hydra::Parameter::Create("wf2_1270").Value(0.1851).Error(0.1);	

	auto F2_1270_amp = BreitWignerAmplitude<MSq12,MSq13,hydra::DWave,2,3,BaBarKSpipi,Helicity>(c_re,c_im,mass,width,radius,radius,phsp,"f2_1270","#it{K}^{0}_{S}#it{f}_{2}(1270)");
	F2_1270_amp.SetColor(kMagenta);
	
	// K0*(1430)- -> Ks pi-
	c_re  = hydra::Parameter::Create("reK0stm_1430").Value(1.62128).Error(0.1);
	c_im  = hydra::Parameter::Create("imK0stm_1430").Value(1.06816).Error(0.1);
	mass  = hydra::Parameter::Create("mK0st_1430").Value(1.459).Error(0.1);
	width = hydra::Parameter::Create("wK0st_1430").Value(0.175).Error(0.1);
	auto lassR = hydra::Parameter::Create("lassR").Value(1.000000).Error(0.1);	
	auto lassB = hydra::Parameter::Create("lassB").Value(0.617734).Error(0.1);	
	auto lassPhiR = hydra::Parameter::Create("lassPhiR").Value(1.104390).Error(0.1);	
	auto lassPhiB = hydra::Parameter::Create("lassPhiB").Value(-0.099519).Error(0.1);	
	auto lassr = hydra::Parameter::Create("lassr").Value(-15.010300).Error(0.1);	
	auto lassa = hydra::Parameter::Create("lassa").Value(0.224004).Error(0.1);	

	auto K0STM_1430_amp = GLassAmplitude<MSq12,MSq13,hydra::SWave,1,3,BaBarKSpipi,Helicity>(c_re,c_im,mass,width,radius,radius,lassR,lassB,lassPhiR,lassPhiB,lassr,lassa,phsp,"K0stm_1430","#it{K}_{0}*(1430)^{#minus}#it{#pi}^{+}");
	K0STM_1430_amp.SetColor(kOrange-3);
	
	// K0*(1430)+ -> Ks pi+
	c_re  = hydra::Parameter::Create("reK0stp_1430").Value(0.148802).Error(0.1);
	c_im  = hydra::Parameter::Create("imK0stp_1430").Value(-0.118513).Error(0.1);

	auto K0STP_1430_amp = GLassAmplitude<MSq12,MSq13,hydra::SWave,1,2,BaBarKSpipi,Helicity>(c_re,c_im,mass,width,radius,radius,lassR,lassB,lassPhiR,lassPhiB,lassr,lassa,phsp,"K0stp_1430","#it{K}_{0}*(1430)^{+}#it{#pi}^{#minus}");
	K0STP_1430_amp.SetColor(kGray+1);
	
	// K2*(1430)- -> Ks pi-
	c_re  = hydra::Parameter::Create("reK2stm_1430").Value(1.15489).Error(0.1);
	c_im  = hydra::Parameter::Create("imK2stm_1430").Value(-0.773363).Error(0.1);
	mass  = hydra::Parameter::Create("mK2st_1430").Value(1.4256).Error(0.1);
	width = hydra::Parameter::Create("wK2st_1430").Value(0.0985).Error(0.1);	

	auto K2STM_1430_amp = BreitWignerAmplitude<MSq12,MSq13,hydra::DWave,1,3,BaBarKSpipi,Helicity>(c_re,c_im,mass,width,radius,radius,phsp,"K2stm_1430","#it{K}_{2}*(1430)^{#minus}#it{#pi}^{+}");
	K2STM_1430_amp.SetColor(kAzure-3);
	K2STM_1430_amp.SetStyle(2);
	
	// K2*(1430)+ -> Ks pi+
	c_re  = hydra::Parameter::Create("reK2stp_1430").Value(0.140865).Error(0.1);
	c_im  = hydra::Parameter::Create("imK2stp_1430").Value(-0.165378).Error(0.1);

	auto K2STP_1430_amp = BreitWignerAmplitude<MSq12,MSq13,hydra::DWave,1,2,BaBarKSpipi,Helicity>(c_re,c_im,mass,width,radius,radius,phsp,"K2stp_1430","#it{K}_{0}*(1430)^{+}#it{#pi}^{#minus}");
	K2STP_1430_amp.SetColor(kGreen+2);
	K2STP_1430_amp.SetStyle(2);

	// K*(1680)- -> Ks pi-
	c_re  = hydra::Parameter::Create("reKstm_1680").Value(1.62128).Error(0.1);
	c_im  = hydra::Parameter::Create("imKstm_1680").Value(1.06816).Error(0.1);
	mass  = hydra::Parameter::Create("mKstm_1680").Value(1.459).Error(0.1);
	width = hydra::Parameter::Create("wKstm_1680").Value(0.175).Error(0.1);	

	auto KSTM_1680_amp = BreitWignerAmplitude<MSq12,MSq13,hydra::PWave,1,3,BaBarKSpipi,Helicity>(c_re,c_im,mass,width,radius,radius,phsp,"Kstm_1680","#it{K}*(1680)^{#minus}#it{#pi}^{+}");
	KSTM_1680_amp.SetColor(kBlue);
	KSTM_1680_amp.SetStyle(2);
	
	// Fvector for pi pi S-Wave
	std::vector<hydra::Parameter> beta;
	beta.push_back(hydra::Parameter::Create("reBeta0").Value(- 5.532630).Error(0.059820));
	beta.push_back(hydra::Parameter::Create("reBeta1").Value( 15.634400).Error(0.057220));
	beta.push_back(hydra::Parameter::Create("reBeta2").Value( 40.865600).Error(1.118000));
	beta.push_back(hydra::Parameter::Create("reBeta3").Value(  6.139390).Error(0.249900));
	beta.push_back(hydra::Parameter::Create("reBeta4").Value(  0.0     ).Error(0.1     ));
	beta.push_back(hydra::Parameter::Create("imBeta0").Value(  0.298584).Error(0.040560));
	beta.push_back(hydra::Parameter::Create("imBeta1").Value(  0.262657).Error(0.070300));
	beta.push_back(hydra::Parameter::Create("imBeta2").Value(-17.780900).Error(0.754000));
	beta.push_back(hydra::Parameter::Create("imBeta3").Value(- 6.929000).Error(0.178100));
	beta.push_back(hydra::Parameter::Create("imBeta4").Value(  0.0     ).Error(0.1     ));

	std::vector<hydra::Parameter> fPr;
	fPr.push_back(hydra::Parameter::Create("refPr_pipi").Value(-11.424400).Error(0.111700));
	fPr.push_back(hydra::Parameter::Create("refPr_kk"  ).Value(- 6.602050).Error(0.367000));
	fPr.push_back(hydra::Parameter::Create("refPr_4pi" ).Value(- 3.795820).Error(0.731100));
	fPr.push_back(hydra::Parameter::Create("refPr_ee"  ).Value(  0.000000).Error(0.000000));
	fPr.push_back(hydra::Parameter::Create("refPr_eep" ).Value(  0.000000).Error(0.000000));
	fPr.push_back(hydra::Parameter::Create("imfPr_pipi").Value(  0.056421).Error(0.104400));
	fPr.push_back(hydra::Parameter::Create("imfPr_kk"  ).Value( 14.002300).Error(0.407000));
	fPr.push_back(hydra::Parameter::Create("imfPr_4pi" ).Value(- 5.828540).Error(0.750500));
	fPr.push_back(hydra::Parameter::Create("imfPr_ee"  ).Value(  0.000000).Error(0.000000));
	fPr.push_back(hydra::Parameter::Create("imfPr_eep" ).Value(  0.000000).Error(0.000000));

	auto s0pr = hydra::Parameter::Create("s0pr").Value( -3.92637).Error(0.1);

	auto fvec_amp = FvectorAmplitude<MSq12,MSq13,2,3>(beta,fPr,s0pr,phsp,"fvec","#it{#pi}^{+}#it{#pi}^{#minus} #it{S}-wave");
	fvec_amp.SetColor(kMagenta);
	fvec_amp.SetStyle(2);
	
	return hydra::sum(KSTM_amp,KSTP_amp,RHO_amp,OMEGA_amp,F2_1270_amp,K0STM_1430_amp,K0STP_1430_amp,K2STM_1430_amp,K2STP_1430_amp,KSTM_1680_amp,fvec_amp);
}

}// namespace D0ToKsPiPi_FVECTOR_BABAR

}// namespace dafne
