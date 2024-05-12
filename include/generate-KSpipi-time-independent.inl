#include <hydra/host/System.h>
#include <hydra/device/System.h>
#include <hydra/Placeholders.h>

#include <iostream>
#include <iomanip>
#include <chrono>

#include <tools/Printer.h>
#include <tools/Plotter.h>
#include <tools/ConfigFile.h>
#include <tools/Arguments.h>
#include <physics/PDG.h>
#include <physics/Rate.h>
#include <physics/Resolution.h>
#include <models/D0ToKsPiPi_FVECTOR_BABAR.h>
using namespace dafne;

// Define the arguments of the amplitude
declarg(MSqPlus , double)
declarg(MSqMinus, double)
declarg(MSqZero , double)
using namespace hydra::arguments;

// output prefix
std::string outprefix = "generate-KSpipi-time-independent";

// Main
int main( int argc, char** argv  )
{
	//---------------------------------------------------------------------------------------
	// Read command line arguments
	//---------------------------------------------------------------------------------------
	StandardArguments args("Generation of D0 -> K0_S pi+ pi- decays according to BABAR 2010 amplitude model","0.1","");
	try {
		args.Read( argc, argv );
	}
	catch ( TCLAP::ArgException& e ) {
		std::cerr << "Error: " << e.error() << " for argument " << e.argId() << "." << std::endl;
		return -1;
	}
	args.Print();

	//---------------------------------------------------------------------------------------
	// Build model for D0->KS pi+ pi- decays
	//---------------------------------------------------------------------------------------
	auto phsp = D0ToKsPiPi_FVECTOR_BABAR::PhaseSpace();

	// build baseline model for the direct decay amplitude
	auto amp = D0ToKsPiPi_FVECTOR_BABAR::Amplitude<MSqPlus,MSqMinus>(phsp);

	// config the model according to configuration file
	ConfigFile config(args.config_file.c_str(), (args.prlevel>3) );
	config.ConfigureModel(amp);

	// efficiency plane
	
	// efficiency plane described by irregular binning 2D histogram
	ArbitraryBinningHistogram2D efficiency_hist = config.ConfigureEfficiencyHistogram();

	// build and check the efficiency plane
	TCanvas cefficiency("cefficiency", "cefficiency", 800, 600);
	gStyle->SetOptStat(0);
	gPad->SetRightMargin(0.15);
	efficiency_hist.GetTH2D((outprefix + "_efficiency_hist").c_str(), 
		                    (outprefix + "_efficiency_hist").c_str(),
		                    "m^{2}_{#it{#pi#pi}} [GeV^{2}/#it{c}^{4}]", "cos(#theta_{#it{#pi#pi}})")->Draw("COLZ");
	Print::Canvas(cefficiency,  args.outdir + outprefix + "_efficiency_hist");
	gStyle->SetOptStat(1);

	auto efficiency = hydra::wrap_lambda(
		[phsp, efficiency_hist] __hydra_dual__ (MSqPlus m2p, MSqMinus m2m) {

		// judge whether in phase space or not
		if (!phsp.Contains<2,3>(m2p, m2m)) return 0.0;

		double m2z = phsp.MSqJK(m2p, m2m);
		double helicity_z = phsp.HelicityJK<2,3>(m2p, m2m);

		double x = m2z;
		double y = helicity_z;
		return efficiency_hist.GetValue(x, y);

	}); 

	// efficiency plane could also be described by an analytic function
	/*auto a0  = hydra::Parameter::Create("eff_a0").Value(0.0).Error(0.1);
	auto a1  = hydra::Parameter::Create("eff_a1").Value(0.0).Error(0.1);
	auto a2  = hydra::Parameter::Create("eff_a2").Value(0.0).Error(0.1);
	auto a3  = hydra::Parameter::Create("eff_a3").Value(0.0).Error(0.1);
	auto center  = hydra::Parameter::Create("eff_center").Value(1.9).Error(0.1);

	auto efficiency = hydra::wrap_lambda(
		[phsp] __hydra_dual__ (unsigned int npar, const hydra::Parameter* params, MSqPlus m2p, MSqMinus m2m) {

		// judge whether in phase space or not
		if (!phsp.Contains<2,3>(m2p, m2m)) return 0.0;

		double center = params[4];
		double x = m2m - center;
		double y = m2p - center;

		double a0 = params[0];
		double a1 = params[1];
		double a2 = params[2];
		double a3 = params[3];

		return 1.0 + a0 * (      x      +      y      ) +
              		 a1 *        x      *      y        +
              		 a2 * ( pow( x, 2 ) + pow( y, 2 ) ) +
            		 a3 * (      x      -      y      );

	}, a0, a1, a2, a3, center); 

	// configure efficiency according to configuration file
	config.Debug();
	config.ConfigureEfficiency(efficiency);*/

	// compute the decay rate
	auto model = rate(amp)*efficiency;

	// prepare the resolution for smearing
	Resolution2D resolution2D;
	config.ConfigureResolution2D(resolution2D);
	TCanvas cresolution("cresolution", "cresolution", 800, 600);
	resolution2D.Draw("m^{2}_{#it{#pi#pi}} [GeV^{2}/#it{c}^{4}]", "cos(#theta_{#it{#pi#pi}})");
	Print::Canvas(cresolution,  args.outdir + outprefix + "_resolution");

	//---------------------------------------------------------------------------------------
	// Generate data
	//---------------------------------------------------------------------------------------
	std::cout << "***** Generation" << std::endl;

	auto start = std::chrono::high_resolution_clock::now();

	auto data = phsp.GenerateData<MSqPlus,MSqMinus,MSqZero>(model,args.nevents,args.seed);
	auto ngenerated = data.size();
	
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> elapsed = end - start;

	std::cout << "Generated " << ngenerated << " data events." << std::endl;
	std::cout << "Time elapsed (ms): "<< elapsed.count() << std::endl;

	//---------------------------------------------------------------------------------------
	// Save to ROOT file
	//---------------------------------------------------------------------------------------
	{
		std::string outfilename = args.outdir + outprefix + "-data.root";
		std::cout << "Saving output to " << outfilename << std::endl;
		TFile *outfile = new TFile(outfilename.c_str(),"recreate");

		double m2p, m2m, m2z, helicity_z; // truth
		double m2p_s, m2m_s, m2z_s, helicity_z_s; // smeared
		int flavor(+1);
		TTree *ntp = new TTree("ntp","ntp");
		ntp->Branch("mSqP",&m2p);
		ntp->Branch("mSqM",&m2m);
		ntp->Branch("mSqZ",&m2z);
		ntp->Branch("HelicityZ",&helicity_z);
		ntp->Branch("mSqP_S",&m2p_s);
		ntp->Branch("mSqM_S",&m2m_s);
		ntp->Branch("mSqZ_S",&m2z_s);
		ntp->Branch("HelicityZ_S",&helicity_z_s);
		ntp->Branch("flavor",&flavor);
		
		for( auto event : data )
		{
			MSqPlus  a = hydra::get<0>(event);
			MSqMinus b = hydra::get<1>(event);
			MSqZero  c = hydra::get<2>(event);
			m2p = a.Value();
			m2m = b.Value();
			m2z = c.Value();
			helicity_z = phsp.HelicityJK<2,3>(m2p, m2m);
			
			do {

				auto smear_result = resolution2D.Smear(m2z, helicity_z);
				m2z_s = std::get<0>(smear_result);
				helicity_z_s = std::get<1>(smear_result);

				auto inv_helicity_z = phsp.invHelicityJK<2,3>(helicity_z_s, m2z_s);
				m2p_s = std::get<0>(inv_helicity_z); 
				m2m_s = std::get<1>(inv_helicity_z); 

			} while (!phsp.Contains<2,3>(m2p_s, m2m_s));

			ntp->Fill();
		}
		
		outfile->Write();
		outfile->Close();
	}
	
	//---------------------------------------------------------------------------------------
	// Plot results
	//---------------------------------------------------------------------------------------
	if (args.plot) {
		std::cout << "***** Plot data and model" << std::endl;

		TApplication* myapp = NULL;
		if (args.interactive) myapp = new TApplication("myapp",0,0);
		
		
		auto plotter = DalitzPlotter<MSqPlus, MSqMinus, MSqZero>(phsp,"#it{K}^{0}_{S}","#it{#pi}^{+}","#it{#pi}^{#minus}",(args.prlevel>3));
		
		std::string outfilename = args.outdir + outprefix + "-HIST.root";
		plotter.FillDataHistogram(data);
		plotter.FillModelHistogram(model);
		plotter.FillComponentHistograms(amp, efficiency, 1);

		// 1D projections
		TCanvas c1("c1","c1",1600,500);
		c1.Divide(3,1);

		c1.cd(1);
		plotter.Plot1DProjections(0, 0);

		c1.cd(2);
		plotter.Plot1DProjections(1, 0);

		c1.cd(3);
		plotter.Plot1DProjections(2, 1); // plot legend in this pad

		outfilename = args.outdir + outprefix + "-1d-projection";
		Print::Canvas(c1,outfilename);

		// 2D projections
		TCanvas c2("c2","c2",1200,500);
		c2.Divide(2,1);

		c2.cd(1);
		plotter.Plot2DProjectionData(0,1);

		c2.cd(2);
		plotter.Plot2DProjectionData(1,2);

		outfilename = args.outdir + outprefix + "-2d-projection";
		Print::Canvas(c2,outfilename);

		// phase difference
		TCanvas c3("c3","c3",600,500);
		c3.SetRightMargin(.14);
		
		auto abar = hydra::wrap_lambda( [&amp] __hydra_dual__ (MSqPlus a, MSqMinus b) {
			MSqPlus switched_a = b;
			MSqMinus switched_b = a;
			return amp(switched_a,switched_b);
		});
		
		plotter.PlotPhaseDifference(amp,abar);
		
		outfilename = args.outdir + outprefix + "-phase-difference";
		Print::Canvas(c3,outfilename);

		if (args.interactive) {
			std::cout << "Press Crtl+C to terminate" << std::endl;
			myapp->Run();
		}
	}
	
	return 0;
}

