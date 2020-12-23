#include <physics/ThreeBodyPhaseSpace.h>
#include <physics/Amplitudes.h>
#include <physics/Utils.h>

#include <TColor.h>

namespace dafne {

//======================================================
// Amplitude model for D0 -> Ks K+ K-
//======================================================
// Based on BaBar's PRL 105, 081803 (2010)
// [http://arxiv.org/abs/1004.5053] 

namespace D0ToKsKK_ISOBAR_BABAR {

auto PhaseSpace()
{
	return ThreeBodyPhaseSpace( 1.8645, {0.49767, 0.49368, 0.49368} );
}

template<typename MSq12, typename MSq13, bool Jordi=true, bool Helicity=false>
__hydra_dual__ inline
auto Amplitude(ThreeBodyPhaseSpace const &phsp)
{
	if (phsp.M<2>()!=phsp.M<3>())
		HYDRA_EXCEPTION("The particles order must be K0_S K+ K-");
	
	auto radius = hydra::Parameter::Create("rBW").Value(1.5).Error(0.1).Fixed();

	// phi -> K+K-
	auto c_re  = hydra::Parameter::Create("re_phi").Value(-1.29497e-01).Error(0.1);
	auto c_im  = hydra::Parameter::Create("im_phi").Value(1.91649e-01).Error(0.1);
	auto mass  = hydra::Parameter::Create("m_phi").Value(1.01955).Error(0.1);
	auto width = hydra::Parameter::Create("w_phi").Value(0.00459787).Error(0.1);

	auto PHI_amp = BreitWignerAmplitude<MSq12,MSq13,hydra::PWave,2,3,Jordi,Helicity>(c_re,c_im,mass,width,radius,radius,phsp,"phi","#it{#phi}(1020)");
	PHI_amp.SetColor(kYellow+1);
	
	// a_0(980)^0 -> K+K-
	mass  = hydra::Parameter::Create("m_a0_980").Value(0.999).Error(0.1);
	width = hydra::Parameter::Create("w_a0_980").Value(0.1135943621).Error(0.1);
	
	auto g1_a0_980 = hydra::Parameter::Create("g1_a0_980").Value(0.6219180253).Error(0.1);
	auto g2_a0_980 = hydra::Parameter::Create("g2_a0_980").Value(0.7830823519).Error(0.1);
	
	auto mEta = hydra::Parameter::Create("m_eta").Value(0.54730).Fixed();
	auto mPi  = hydra::Parameter::Create("w_pi").Value(0.139570).Fixed();
	
	c_re  = hydra::Parameter::Create("re_a00_980").Value(12.2344861).Error(0.1).Fixed();
	c_im  = hydra::Parameter::Create("im_a00_980").Value(0.0).Error(0.1).Fixed();
	
	auto A00_980_amp = FlatteAmplitude<MSq12,MSq13,hydra::SWave,2,3,Jordi,Helicity>(c_re,c_im,mass,width,radius,radius,g1_a0_980,g2_a0_980,mEta,mPi,phsp,"a00_980","#it{a}_{0}(980)^{0}");
	A00_980_amp.SetColor(kBlue);
	
	// a_0(980)^+  -> K+KS
	c_re  = hydra::Parameter::Create("re_a0p_980").Value(-7.56633713).Error(0.1);
	c_im  = hydra::Parameter::Create("im_a0p_980").Value(-1.73609903).Error(0.1);
	
	auto A0p_980_amp = FlatteAmplitude<MSq12,MSq13,hydra::SWave,1,2,Jordi,Helicity>(c_re,c_im,mass,width,radius,radius,g1_a0_980,g2_a0_980,mEta,mPi,phsp,"a0p_980","#it{a}_{0}(980)^{+}");
	A0p_980_amp.SetColor(kMagenta);
	
	// a_0(980)^-  -> KSK-
	c_re  = hydra::Parameter::Create("re_a0m_980").Value(-1.19249603).Error(0.1);
	c_im  = hydra::Parameter::Create("im_a0m_980").Value(0.94854397).Error(0.1);
	
	auto A0m_980_amp = FlatteAmplitude<MSq12,MSq13,hydra::SWave,1,3,Jordi,Helicity>(c_re,c_im,mass,width,radius,radius,g1_a0_980,g2_a0_980,mEta,mPi,phsp,"a0m_980","#it{a}_{0}(980)^{#minus}");
	A0m_980_amp.SetColor(kGreen+2);
	
	// a_0(1450)^0 -> K+K-
	mass  = hydra::Parameter::Create("m_a00_1450").Value(1.474).Error(0.1);
	width = hydra::Parameter::Create("w_a00_1450").Value(0.265).Error(0.1);

	c_re  = hydra::Parameter::Create("re_a00_1450").Value(-2.91593e-01).Error(0.1);
	c_im  = hydra::Parameter::Create("im_a00_1450").Value(-7.72051e-01).Error(0.1);
		
	auto A00_1450_amp = BreitWignerAmplitude<MSq12,MSq13,hydra::SWave,2,3,Jordi,Helicity>(c_re,c_im,mass,width,radius,radius,phsp,"a00_1450","#it{a}_{0}(1450)^{0}");
	A00_1450_amp.SetColor(kBlue);
	A00_1450_amp.SetStyle(2);

	// a_0(1450)^+ -> K+KS
	c_re  = hydra::Parameter::Create("re_a0p_1450").Value(-8.72480e-02).Error(0.1);
	c_im  = hydra::Parameter::Create("im_a0p_1450").Value(9.26874e-01).Error(0.1);
	
	auto A0p_1450_amp = BreitWignerAmplitude<MSq12,MSq13,hydra::SWave,1,2,Jordi,Helicity>(c_re,c_im,mass,width,radius,radius,phsp,"a0p_1450","#it{a}_{0}(1450)^{+}");
	A0p_1450_amp.SetColor(kMagenta);
	A0p_1450_amp.SetStyle(2);
	
	// a_0(1450)^- -> KSK-
	c_re  = hydra::Parameter::Create("re_a0m_1450").Value(0.0).Error(0.1).Fixed();
	c_im  = hydra::Parameter::Create("im_a0m_1450").Value(0.0).Error(0.1).Fixed();
	
	auto A0m_1450_amp = BreitWignerAmplitude<MSq12,MSq13,hydra::SWave,1,3,Jordi,Helicity>(c_re,c_im,mass,width,radius,radius,phsp,"a0m_1450","#it{a}_{0}(1450)^{#minus}");
	A0m_1450_amp.SetColor(kGreen+2);
	A0m_1450_amp.SetStyle(2);
	
	A0m_1450_amp.Remove();

	// f0(1370) -> K+K-
	c_re  = hydra::Parameter::Create("re_f0_1370").Value(1.55304e-01).Error(0.1);
	c_im  = hydra::Parameter::Create("im_f0_1370").Value(2.37891e-02).Error(0.1);
	mass  = hydra::Parameter::Create("m_f0_1370").Value(1.434).Error(0.1);
	width = hydra::Parameter::Create("w_f0_1370").Value(0.265).Error(0.1);

	auto F0_1370_amp = BreitWignerAmplitude<MSq12,MSq13,hydra::SWave,2,3,Jordi,Helicity>(c_re,c_im,mass,width,phsp,"f0_1370","#if{f}_{0}(1370)");
	F0_1370_amp.SetColor(kOrange-3);
	
	// f2(1270) -> K+K-
	c_re  = hydra::Parameter::Create("re_f2_1270").Value(3.84159e-01).Error(0.1);
	c_im  = hydra::Parameter::Create("im_f2_1270").Value(2.36121e-02).Error(0.1);
	mass  = hydra::Parameter::Create("m_f2_1270").Value(1.2754).Error(0.1);
	width = hydra::Parameter::Create("w_f2_1270").Value(0.1851).Error(0.1);

	auto F2_1270_amp = BreitWignerAmplitude<MSq12,MSq13,hydra::DWave,2,3,Jordi,Helicity>(c_re,c_im,mass,width,phsp,"f2_1270","#if{f}_{2}(1370)");
	F2_1270_amp.SetColor(kGray+1);
	
	return hydra::sum(PHI_amp,A00_980_amp,A0p_980_amp,A0m_980_amp,A00_1450_amp,A0p_1450_amp,A0m_1450_amp,F0_1370_amp,F2_1270_amp);
}

}// namespace D0ToKsKK_ISOBAR_BABAR

}// namespace dafne
