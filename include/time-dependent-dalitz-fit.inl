#include <hydra/host/System.h>
#include <hydra/device/System.h>
#include <hydra/FunctorArithmetic.h>
#include <hydra/Algorithm.h>
#include <hydra/Pdf.h>
#include <hydra/LogLikelihoodFCN.h>

#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnUserParameterState.h>
#include <Minuit2/MnPrint.h>
#include <Minuit2/MnMigrad.h>

#include <TFile.h>
#include <TTree.h>

#include <physics/PDG.h>
#include <physics/Amplitudes.h>
#include <physics/Rate.h>
using namespace dafne;

// Configuration parameters
const unsigned nevents(100000);
std::string outdir("./");

// Define the arguments of the amplitude
declarg(MSqZero  , double)
declarg(MSqPlus  , double)
declarg(MSqMinus , double)
declarg(DecayTime, double)
using namespace hydra::arguments;

// Main
int main()
{
	// Phase space for D0 -> KS pi+ pi- decay with decay-time in [0,10tau]
	ThreeBodyPhaseSpaceWithTime phsp(Mass::D0,{Mass::K0,Mass::Pi,Mass::Pi},{0.,10.*Tau::D0});

	// Radii of the Blatt-Weisskopf barrier factors
	auto radiusReso   = hydra::Parameter::Create("radiusReso").Value(1.5).Error(0.1).Fixed();
	auto radiusMother = hydra::Parameter::Create("radiusMother").Value(5.0).Error(0.1).Fixed();

	// rho(770)0 -> pi+pi-
	auto c_re  = hydra::Parameter::Create("RHO_C_Re ").Value(1.).Fixed();
	auto c_im  = hydra::Parameter::Create("RHO_C_Im ").Value(0.).Fixed();
	auto mass  = hydra::Parameter::Create("RHO_Mass ").Value(0.7758).Fixed();
	auto width = hydra::Parameter::Create("RHO_Width").Value(0.1464).Fixed();

	auto RHO_amp = GounarisSakuraiAmplitude<MSqPlus,MSqMinus,hydra::PWave,2,3>(c_re,c_im,mass,width,radiusReso,radiusMother,phsp);

	// K*(892)- -> KS pi-
	c_re  = hydra::Parameter::Create("KSTM_C_Re ").Value(-1.1961).Error(0.0001);
	c_im  = hydra::Parameter::Create("KSTM_C_Im ").Value( 1.1247).Error(0.0001);
	mass  = hydra::Parameter::Create("KSTM_Mass ").Value(0.8936060).Fixed();
	width = hydra::Parameter::Create("KSTM_Width").Value(0.0463407).Fixed();
	
	auto KSTM_amp = BreitWignerAmplitude<MSqPlus,MSqMinus,hydra::PWave,1,3>(c_re,c_im,mass,width,radiusReso,radiusMother,phsp);

	// Non resonant
	c_re  = hydra::Parameter::Create("NR_C_Re").Value(0.5).Error(0.0001);
	c_im  = hydra::Parameter::Create("NR_C_Im").Value(0.1).Error(0.0001);
	auto NR_amp = NonResonantAmplitude<MSqPlus,MSqMinus>(c_re,c_im,phsp);
	
	// Total amplitude for the direct decay
	auto Adir = hydra::sum(RHO_amp,KSTM_amp,NR_amp);
	
	// Total amplitude for the decay that undergoes mixing
	auto Abar = hydra::wrap_lambda( [&Adir] __hydra_dual__ (MSqPlus a, MSqMinus b) {
		MSqPlus switched_a = b;
		MSqMinus switched_b = a;
		return Adir(switched_a,switched_b);
	});
	
	// Time-dependent parameters
	auto x   = hydra::Parameter::Create("x").Value(0.01).Error(0.0001).Limits(-1.,1.);
	auto y   = hydra::Parameter::Create("y").Value(0.06).Error(0.0001).Limits(-1.,1.);
	auto tau = hydra::Parameter::Create("tau").Value(Tau::D0).Fixed();
	
	// CP-violation parameters
	auto qop = hydra::Parameter::Create("qop").Value(1.0).Error(0.0001).Limits(0.,10.).Fixed();
	auto phi = hydra::Parameter::Create("phi").Value(0.0).Error(0.0001).Limits(-1.,1.).Fixed();
	
	// D0 rate
	auto model_dz = time_dependent_rate<Flavor::Positive,DecayTime>(tau,x,y,qop,phi,Adir,Abar);
	
	// D0bar rate
	auto model_db = time_dependent_rate<Flavor::Negative,DecayTime>(tau,x,y,qop,phi,Adir,Abar);

	//---------------------------------------------------------------------------------------
	// Generate input data
	//---------------------------------------------------------------------------------------
	std::cout << "***** Generation" << std::endl;
	auto data_dz = phsp.GenerateDataWithTime<MSqPlus,MSqMinus,MSqZero,DecayTime>(model_dz,nevents,12345);
	std::cout << "Generated " << data_dz.size() << " D0 candidates." << std::endl;
	auto data_db = phsp.GenerateDataWithTime<MSqPlus,MSqMinus,MSqZero,DecayTime>(model_db,nevents,67890);
	std::cout << "Generated " << data_db.size() << " D0bar candidates." << std::endl;

	//---------------------------------------------------------------------------------------
	// Save to ROOT file
	//---------------------------------------------------------------------------------------
	{
		std::string outfilename = outdir + "time-dependent-dalitz-data.root";
		std::cout << "Saving output to " << outfilename << std::endl;
		TFile *outfile = new TFile(outfilename.c_str(),"recreate");

		double m2_12, m2_13, m2_23, t;
		int flavor(+1);
		TTree *ntp = new TTree("ntp","ntp");
		ntp->Branch("mSq12",&m2_12);
		ntp->Branch("mSq13",&m2_13);
		ntp->Branch("mSq23",&m2_23);
		ntp->Branch("t",&t);
		ntp->Branch("flavor",&flavor);
		
		for( auto event : data_dz )
		{
			MSqPlus   a = hydra::get<0>(event);
			MSqMinus  b = hydra::get<1>(event);
			MSqZero   c = hydra::get<2>(event);
			DecayTime d = hydra::get<3>(event);
			
			m2_12 = a.Value();
			m2_13 = b.Value();
			m2_23 = c.Value();
			t     = d.Value();
			
			ntp->Fill();
		}

		flavor = -1;
		for( auto event : data_db )
		{
			MSqPlus   a = hydra::get<0>(event);
			MSqMinus  b = hydra::get<1>(event);
			MSqZero   c = hydra::get<2>(event);
			DecayTime d = hydra::get<3>(event);
			
			m2_12 = a.Value();
			m2_13 = b.Value();
			m2_23 = c.Value();
			t     = d.Value();
			
			ntp->Fill();
		}

		outfile->Write();
		outfile->Close();
	}
	
	//---------------------------------------------------------------------------------------
	// Build pdf and log-likelihood function from model
	//---------------------------------------------------------------------------------------
	std::cout << "***** Build likelihood" << std::endl;
	auto pdf_dz = hydra::make_pdf( model_dz, phsp.IntegratorWithTime() );
	std::cout << "Initial normalization for D0 PDF: "<< pdf_dz.GetNorm() << " +/- " << pdf_dz.GetNormError() << std::endl;
	
	auto pdf_db = hydra::make_pdf( model_db, phsp.IntegratorWithTime() );
	std::cout << "Initial normalization for D0bar PDF: "<< pdf_db.GetNorm() << " +/- " << pdf_db.GetNormError() << std::endl;
	
	auto fcn_dz = hydra::make_loglikehood_fcn( pdf_dz, data_dz.begin(), data_dz.end() );
	auto fcn_db = hydra::make_loglikehood_fcn( pdf_db, data_db.begin(), data_db.end() );
	
	auto fcn = hydra::make_simultaneous_fcn(fcn_dz, fcn_db);
	
	//---------------------------------------------------------------------------------------
	// Configure and run MINUIT
	//---------------------------------------------------------------------------------------
	std::cout << "***** Fit" << std::endl;
	ROOT::Minuit2::MnPrint::SetLevel(2);
	ROOT::Minuit2::MnStrategy strategy(2);
	
	// print starting values of parameters
	auto pstart = fcn.GetParameters().GetMnState();
	std::cout << "Starting parameters:" << pstart << std::endl;
		
	ROOT::Minuit2::MnMigrad migrad(fcn, pstart, strategy);
	auto minimum = ROOT::Minuit2::FunctionMinimum( migrad(5000., 1.) );

	std::cout << minimum << std::endl;

	if ( !minimum.IsValid() || !minimum.HasAccurateCovar() ) {
		std::cout << "Fit did not converge or covariance matrix is not accurate." << std::endl;
		return -2;
	}

	auto parameters = minimum.UserParameters();
	auto covariance = minimum.UserCovariance();
	std::cout << "***** Fit results:\n" << parameters << covariance << std::endl;
	
	return 0;
}

