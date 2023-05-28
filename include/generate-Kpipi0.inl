#include <hydra/host/System.h>
#include <hydra/device/System.h>
#include <hydra/Placeholders.h>

#include <iostream>
#include <iomanip>

#include <tools/Plotter.h>
#include <tools/Printer.h>
#include <tools/ConfigFile.h>
#include <tools/Arguments.h>
#include <physics/PDG.h>
#include <physics/Rate.h>
#include <models/D0ToKPiPi0RS_CLEO.h>
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
	StandardArguments args("Generation of D0 -> K- pi+ pi0 decays according to CLEO model","0.1","");
	try {
		args.Read( argc, argv );
	}
	catch ( TCLAP::ArgException& e ) {
		std::cerr << "Error: " << e.error() << " for argument " << e.argId() << "." << std::endl;
		return -1;
	}
	args.Print();

	//---------------------------------------------------------------------------------------
	// Build model for D0->K-pi+pi0 decays
	//---------------------------------------------------------------------------------------
	auto phsp = D0ToKPiPi0RS_CLEO::PhaseSpace();
	auto amp  = D0ToKPiPi0RS_CLEO::Amplitude<MSqZero,MSqMinus>(phsp);
	
	// config the model according to configuration file
	ConfigFile config(args.config_file.c_str(), (args.prlevel>3) );
	config.ConfigureModel(amp);

	// compute the decay rate
	auto model = rate(amp);

	//---------------------------------------------------------------------------------------
	// Obtain individual components
	//---------------------------------------------------------------------------------------
	auto  rho_770_p_amp = amp.GetFunctor(hydra::placeholders::_0);
	auto  kst_892_m_amp = amp.GetFunctor(hydra::placeholders::_1);
	auto  kst_892_0_amp = amp.GetFunctor(hydra::placeholders::_2);
	auto  k0_1430_m_amp = amp.GetFunctor(hydra::placeholders::_3);
	auto  k0_1430_0_amp = amp.GetFunctor(hydra::placeholders::_4);
	auto rho_1700_p_amp = amp.GetFunctor(hydra::placeholders::_5);
	auto kst_1680_m_amp = amp.GetFunctor(hydra::placeholders::_6);
	auto         nr_amp = amp.GetFunctor(hydra::placeholders::_7);
	
	//---------------------------------------------------------------------------------------
	// Compute fit fractions
	//---------------------------------------------------------------------------------------
	std::cout << "***** Fit fractions:" << std::endl;
	auto  rho_770_p_FF = phsp.FitFraction( model,  rate(rho_770_p_amp) );
	auto  kst_892_m_FF = phsp.FitFraction( model,  rate(kst_892_m_amp) );
	auto  kst_892_0_FF = phsp.FitFraction( model,  rate(kst_892_0_amp) );
	auto  k0_1430_m_FF = phsp.FitFraction( model,  rate(k0_1430_m_amp) );
	auto  k0_1430_0_FF = phsp.FitFraction( model,  rate(k0_1430_0_amp) );
	auto rho_1700_p_FF = phsp.FitFraction( model, rate(rho_1700_p_amp) );
	auto kst_1680_m_FF = phsp.FitFraction( model, rate(kst_1680_m_amp) );
	auto         nr_FF = phsp.FitFraction( model,         rate(nr_amp) );
	auto  rho_770_p_kst_892_0_Interference = phsp.Interference(model, rho_770_p_amp, kst_892_0_amp);
	auto  rho_770_p_k0_1430_0_Interference = phsp.Interference(model, rho_770_p_amp, k0_1430_0_amp);
	std::cout << std::setprecision(6) << std::fixed;
	std::cout <<   rho_770_p_amp.Name() << " = " <<  rho_770_p_FF.first << " +- " <<  rho_770_p_FF.second << std::endl;
	std::cout <<   kst_892_m_amp.Name() << " = " <<  kst_892_m_FF.first << " +- " <<  kst_892_m_FF.second << std::endl;
	std::cout <<   kst_892_0_amp.Name() << " = " <<  kst_892_0_FF.first << " +- " <<  kst_892_0_FF.second << std::endl;
	std::cout <<   k0_1430_m_amp.Name() << " = " <<  k0_1430_m_FF.first << " +- " <<  k0_1430_m_FF.second << std::endl;
	std::cout <<   k0_1430_0_amp.Name() << " = " <<  k0_1430_0_FF.first << " +- " <<  k0_1430_0_FF.second << std::endl;
	std::cout <<  rho_1700_p_amp.Name() << " = " << rho_1700_p_FF.first << " +- " << rho_1700_p_FF.second << std::endl;
	std::cout <<  kst_1680_m_amp.Name() << " = " << kst_1680_m_FF.first << " +- " << kst_1680_m_FF.second << std::endl;
	std::cout <<          nr_amp.Name() << " = " <<         nr_FF.first << " +- " <<         nr_FF.second << std::endl;

	std::cout << rho_770_p_amp.Name() << ", " << kst_892_0_amp.Name() << " interference ";
	std::cout << " = " <<         rho_770_p_kst_892_0_Interference.first << " +- ";
	std::cout <<     rho_770_p_kst_892_0_Interference.second << std::endl;
	std::cout << rho_770_p_amp.Name() << ", " << k0_1430_0_amp.Name() << " interference ";
	std::cout << " = " <<         rho_770_p_k0_1430_0_Interference.first << " +- ";
	std::cout <<     rho_770_p_k0_1430_0_Interference.second << std::endl;
	//---------------------------------------------------------------------------------------
	// Generate data
	//---------------------------------------------------------------------------------------
	
	std::cout << "***** Generation" << std::endl;
	auto start = std::chrono::high_resolution_clock::now();	

	auto data = phsp.GenerateData<MSqZero,MSqMinus,MSqPlus>(model,args.nevents,args.seed);
	auto ngenerated = data.size();
	std::cout << "Generated " << ngenerated << " data events." << std::endl;

	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> elapsed = end - start;
	std::cout << "Time elapsed (ms):"<< elapsed.count() << std::endl;

	
	//---------------------------------------------------------------------------------------
	// Save to ROOT file
	//---------------------------------------------------------------------------------------
	{
		std::string outfilename = args.outdir + "Kpipi0-data.root";
		std::cout << "Saving output to " << outfilename << std::endl;
		TFile *outfile = new TFile(outfilename.c_str(),"recreate");

		double m2_12, m2_13, m2_23;
		TTree *ntp = new TTree("ntp","ntp");
		ntp->Branch("mSqZ",&m2_12);
		ntp->Branch("mSqM",&m2_13);
		ntp->Branch("mSqP",&m2_23);
		
		for( auto event : data )
		{
			MSqZero  a = hydra::get<0>(event);
			MSqMinus b = hydra::get<1>(event);
			MSqPlus  c = hydra::get<2>(event);
			m2_12 = a.Value();
			m2_13 = b.Value();
			m2_23 = c.Value();

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

		auto plotter = DalitzPlotter<MSqZero,MSqMinus,MSqPlus>(phsp,"#it{K}^{#minus}","#it{#pi}^{+}","#it{#pi}^{0}",(args.prlevel>3));
		plotter.FillDataHistogram(data, 140, 140, 100, 0.3, 3.1, 0.3, 3.1, 0.0, 1.9);
		plotter.FillModelHistogram(model, 140, 140, 100, 0.3, 3.1, 0.3, 3.1, 0.0, 1.9);
		auto efficiency = ConstantFunctor(1);
		double signal_fraction = 1; // pure signal
		plotter.FillComponentHistograms(amp, efficiency, signal_fraction, 140, 140, 100, 0.3, 3.1, 0.3, 3.1, 0.0, 1.9);
		
		// 2D data projections
		TCanvas c1("c1","c1",1600,500);
		c1.Divide(3,1);
		
		c1.cd(1);
		plotter.Plot2DProjectionData(0,1);
		
		c1.cd(2);
		plotter.Plot2DProjectionData(0,2);
		
		c1.cd(3);
		plotter.Plot2DProjectionData(2,1);
		
		std::string outfilename = args.outdir + "Kpipi0-2d-projection";
		Print::Canvas(c1, outfilename);
		
		// 1D projections
		TCanvas c2("c2","c2",1600,500);
		c2.Divide(3,1);
		
		c2.cd(1);
		gPad->SetRightMargin(0.07);
		gPad->SetLeftMargin(0.15);
		plotter.Plot1DProjections(0, 0);

		c2.cd(2);
		gPad->SetRightMargin(0.07);
		gPad->SetLeftMargin(0.15);
		plotter.Plot1DProjections(1, 0);

		c2.cd(3);
		gPad->SetRightMargin(0.07);
		gPad->SetLeftMargin(0.15);
		plotter.Plot1DProjections(2, 1); // plot legend in this pad
		
		outfilename = args.outdir + "Kpipi0-1d-projection";
		Print::Canvas(c2, outfilename);

		if (args.interactive) {
			std::cout << "Press Crtl+C to terminate" << std::endl;
			myapp->Run();
		}
	}
	
	return 0;
}

