#include <hydra/host/System.h>
#include <hydra/device/System.h>
#include <hydra/Pdf.h>
#include <hydra/LogLikelihoodFCN.h>
#include <hydra/Placeholders.h>

#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnUserParameterState.h>
#include <Minuit2/MnPrint.h>
#include <Minuit2/MnMigrad.h>

#include <tools/Plotter.h>
#include <tools/Arguments.h>
#include <physics/PDG.h>
#include <physics/Amplitudes.h>
#include <physics/Rate.h>
using namespace dafne;

// Define the arguments of the amplitude
declarg(MSqZero , double)
declarg(MSqPlus , double)
declarg(MSqMinus, double)
using namespace hydra::arguments;

// Main
int main( int argc, char** argv  )
{
	//---------------------------------------------------------------------------------------
	// Read command line arguments
	//---------------------------------------------------------------------------------------
	Arguments args("Dalitz-fit example","0.1","");
	try {
		args.Read( argc, argv );
	}
	catch ( TCLAP::ArgException& e ) {
		std::cerr << "Error: " << e.error() << " for argument " << e.argId() << "." << std::endl;
		return -1;
	}
	args.Print();
	
	//---------------------------------------------------------------------------------------
	// Build simplistic model for D0->pi+pi-pi0 with only rho(770) resonances and nonresonant
	//---------------------------------------------------------------------------------------
	ThreeBodyPhaseSpace phsp( Mass::D0, {Mass::Pi, Mass::Pi, Mass::Pi0} );
	
	auto radiusReso   = hydra::Parameter::Create("radiusReso").Value(1.5).Error(0.1).Fixed();
	auto radiusMother = hydra::Parameter::Create("radiusMother").Value(5.0).Error(0.1).Fixed();
	
	auto mass  = hydra::Parameter::Create("RHO_770_mass  ").Value(0.7758).Fixed();
	auto width = hydra::Parameter::Create("RHO_770_width ").Value(0.1503).Fixed();

	// rho(770)+ -> pi+pi0
	auto cRe = hydra::Parameter::Create("RHO_770_P_cRe ").Value(1.).Fixed();
	auto cIm = hydra::Parameter::Create("RHO_770_P_cIm ").Value(0.).Fixed();
	
	auto rho_770_p_amp = BreitWignerAmplitude<MSqZero,MSqPlus,hydra::PWave,1,3>(cRe,cIm,mass,width,radiusReso,radiusMother,phsp,"rho_770_p","#rho(770)^{+}");
	rho_770_p_amp.SetColor(kBlue);

	// rho(770)0 -> pi+pi-
	cRe = hydra::Parameter::Create("RHO_770_0_cRe ").Value(+0.53).Error(0.01).Limits(-1.,1.);
	cIm = hydra::Parameter::Create("RHO_770_0_cIm ").Value(-0.15).Error(0.01).Limits(-1.,1.);
	
	auto rho_770_0_amp = BreitWignerAmplitude<MSqZero,MSqPlus,hydra::PWave,1,2>(cRe,cIm,mass,width,radiusReso,radiusMother,phsp,"rho_770_0","#rho(770)^{0}");
	rho_770_0_amp.SetColor(kMagenta);
	
	// rho(770)- -> pi-pi0
	cRe = hydra::Parameter::Create("RHO_770_M_cRe ").Value(0.73).Error(0.01).Limits(-1.,1.);
	cIm = hydra::Parameter::Create("RHO_770_M_cIm ").Value(0.02).Error(0.01).Limits(-1.,1.);
	
	auto rho_770_m_amp = BreitWignerAmplitude<MSqZero,MSqPlus,hydra::PWave,2,3>(cRe,cIm,mass,width,radiusReso,radiusMother,phsp,"rho_770_m","#rho(770)^{#minus}");
	rho_770_m_amp.SetColor(kGreen+2);
	
	// NR
	cRe = hydra::Parameter::Create("NR_cRe        ").Value(+0.40).Error(0.01).Limits(-1.,1.);
	cIm = hydra::Parameter::Create("NR_cIm        ").Value(-0.16).Error(0.01).Limits(-1.,1.);
	auto nr_amp = NonResonantAmplitude<MSqZero,MSqPlus>(cRe,cIm,phsp,"Nonresonant");
	nr_amp.SetStyle(2);
	
	// Total model
	auto model = rate(rho_770_p_amp, rho_770_0_amp, rho_770_m_amp, nr_amp);
	
	//---------------------------------------------------------------------------------------
	// Generate input data
	//---------------------------------------------------------------------------------------
	std::cout << "***** Generation" << std::endl;
	auto data = phsp.GenerateData<MSqZero,MSqPlus,MSqMinus>(model,args.nevents,args.seed);
	auto ncands = data.size();
	std::cout << "Generated " << ncands << " data events." << std::endl;
	
	//---------------------------------------------------------------------------------------
	// Build pdf and log-likelihood function from model
	//---------------------------------------------------------------------------------------
	std::cout << "***** Build likelihood" << std::endl;
	auto model_pdf = hydra::make_pdf( model, phsp.Integrator() );
	std::cout << "Initial PDF normalization: "<< model_pdf.GetNorm() << " +/- " << model_pdf.GetNormError() << std::endl;
	
	auto fcn = hydra::make_loglikehood_fcn( model_pdf, data.begin(), data.end() );
	
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
	
	//---------------------------------------------------------------------------------------
	// Plot results
	//---------------------------------------------------------------------------------------
	if (args.plot) {
		TApplication myapp("myapp",0,0);

		auto plotter = DalitzPlotter<MSqZero,MSqPlus,MSqMinus>(phsp,"#it{#pi}^{+}","#it{#pi}^{#minus}","#it{#pi}^{0}",(args.prlevel>3));
		
		// make sure parameters are updated to their value at the minimum
		fcn.GetParameters().UpdateParameters(minimum);
		
		// get the pdf, then model functor, and finally the amplitude functor
		auto fit_model = fcn.GetPDF().GetFunctor().GetFunctor(hydra::placeholders::_1);
		
		plotter.FillHistograms(data,fit_model);

		// 2D data projections
		TCanvas c1("c1","c1",1600,500);
		c1.Divide(3,1);
		
		c1.cd(1);
		plotter.Plot2DProjectionData(0,1);
		
		c1.cd(2);
		plotter.Plot2DProjectionData(0,2);
		
		c1.cd(3);
		plotter.Plot2DProjectionData(2,1);
		
		std::string outfilename = args.outdir + "Kpipi0-2d-projection.pdf";
		c1.Print(outfilename.c_str());
		
		// 1D projections
		TCanvas c2("c2","c2",1600,500);
		c2.Divide(3,1);
		
		c2.cd(1);
		plotter.Plot1DProjections(0, 0);
		c2.cd(1)->SetLogy();

		c2.cd(2);
		plotter.Plot1DProjections(1, 0);
		c2.cd(2)->SetLogy();

		c2.cd(3);
		plotter.Plot1DProjections(2, 1); // plot legend in this pad
		c2.cd(3)->SetLogy();
		
		outfilename = args.outdir + "Kpipi0-1d-projection.pdf";
		c2.Print(outfilename.c_str());

		if (args.interactive) {
			std::cout << "Press Crtl+C to terminate" << std::endl;
			myapp.Run();
		}
	}
	
	return 0;
}

