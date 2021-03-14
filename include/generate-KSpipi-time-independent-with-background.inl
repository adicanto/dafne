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
#include <tools/FunctorTools.h>
#include <models/D0ToKsPiPi_FVECTOR_BABAR.h>

#include <hydra/Function.h>
#include <hydra/FunctorArithmetic.h>
#include <hydra/Lambda.h>

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
	Print::Canvas(cefficiency,  args.outdir + outprefix + "efficiency_hist");
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
	auto model_truth = rate(amp);
	auto model = model_truth*efficiency;
	auto plainIntegrator = phsp.Integrator(10*args.nevents);
	auto normalized_model = make_numerical_normalized_functor(model, plainIntegrator);

	// add background

	// auto f_rnd = hydra::Parameter::Create("f_rnd").Value(0.1).Error(0.0001);
	auto f_cmb = hydra::Parameter::Create("f_cmb").Value(0.5).Error(0.0001).Limits(-1.,1.);

	// auto random_background = [phsp](MSqPlus m2p, MSqMinus m2m)->double {
	// 		// judge whether in phase space or not
	// 		if (!phsp.Contains<2,3>(m2p, m2m)) return 0.0;
	// 		return 1;
	// };

	auto combinatorial_background = hydra::wrap_lambda( 
		[phsp] __hydra_dual__ (MSqPlus m2p, MSqMinus m2m) {
			// judge whether in phase space or not
			if (!phsp.Contains<2,3>(m2p, m2m)) return 0.0;
			return 1.;
		}
	);
	
	auto combinatorial_background_norm = plainIntegrator(combinatorial_background).first;
	auto combinatorial_background_pdf =  divideBy<double>(combinatorial_background, combinatorial_background_norm); 



	// build averaged sum pdf for the plotting

	auto f_cmb_functor = SingleValue(f_cmb);

	auto _build_averaged_sum_pdf = hydra::wrap_lambda(
			  [] __hydra_dual__ (hydra::tuple< double, double, double> input_functors){
			  		auto _pdf_sig = hydra::get<0>(input_functors);
			  		auto _f_cmb = hydra::get<1>(input_functors);
			  		auto _pdf_cmb = hydra::get<2>(input_functors);

			  		return  (1-_f_cmb)*_pdf_sig + _f_cmb*_pdf_cmb;
			  }
	);

	auto averaged_sum_pdf = hydra::compose(_build_averaged_sum_pdf, normalized_model, f_cmb_functor, combinatorial_background_pdf);

	std::cout << "combinatorial_background_pdf(1., 1.): " << combinatorial_background_pdf(1., 1.) << std::endl ; // for debug
	std::cout << "normalized_model(1., 1.): " << normalized_model(1., 1.) << std::endl ; // for debug
	std::cout << "(1-f_cmb) * normalized_model(1., 1.) + f_cmb * combinatorial_background_pdf(1., 1.): ";
	std::cout << (1-f_cmb.GetValue()) * normalized_model(MSqPlus(1.), MSqMinus(1.)) 
				 + f_cmb.GetValue() * combinatorial_background_pdf(MSqPlus(1.), MSqMinus(1.)) << std::endl; // for debug
	std::cout << "averaged_sum_pdf(1., 1.): " << averaged_sum_pdf(MSqPlus(1.), MSqMinus(1.)) << std::endl ; // for debug
	// for debugging
	// typename decltype(averaged_sum_pdf)::argument_type  test{};
	// std::cout << test.dummy << '\n';


	//---------------------------------------------------------------------------------------
	// Generate data
	//---------------------------------------------------------------------------------------
	std::cout << "***** Generation" << std::endl;

	hydra::multivector<hydra::tuple<MSqPlus,MSqMinus,MSqZero>, hydra::device::sys_t> data;

	std::cout << "Generating signal ... ... " << std::endl;
	auto start = std::chrono::high_resolution_clock::now();

	auto data_sig = phsp.GenerateData<MSqPlus,MSqMinus,MSqZero>(model,(1-f_cmb.GetValue())*args.nevents,args.seed);
	auto ngenerated = data_sig.size();
	
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> elapsed = end - start;
	std::cout << "Generated " << ngenerated << " signal events." << std::endl;
	std::cout << "Time elapsed (ms):"<< elapsed.count() << std::endl;


	std::cout << "Generating combinatorial background ... ... " << std::endl;
	start = std::chrono::high_resolution_clock::now();

	auto data_cmb = phsp.GenerateData<MSqPlus,MSqMinus,MSqZero>(combinatorial_background,f_cmb.GetValue()*args.nevents,args.seed);
	ngenerated = data_cmb.size();
	
	end = std::chrono::high_resolution_clock::now();
	elapsed = end - start;
	std::cout << "Generated " << ngenerated << " combinatorial background events." << std::endl;
	std::cout << "Time elapsed (ms):"<< elapsed.count() << std::endl;

	// merge and shuffle
	data.insert( data.end(), data_sig.begin(), data_sig.end());
	data.insert( data.end(), data_cmb.begin(), data_cmb.end());
	std::random_shuffle(data.begin(), data.end());
	ngenerated = data.size();
	std::cout << "Generated " << ngenerated << " events in total." << std::endl;



	//---------------------------------------------------------------------------------------
	// Save to ROOT file
	//---------------------------------------------------------------------------------------

	// prepare the resolution for smearing
	Resolution2D resolution2D;
	config.ConfigureResolution2D(resolution2D);

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

		// plotting procedure
		TApplication myapp("myapp",0,0);
		
		auto plotter = DalitzPlotter<MSqPlus, MSqMinus, MSqZero>(phsp,"#it{K}^{0}_{S}","#it{#pi}^{+}","#it{#pi}^{#minus}",(args.prlevel>3));
		
		std::string outfilename = args.outdir + outprefix + "-HIST.root";
		// plotter.FillHistograms(data, amp, outfilename); // currently, the efficiency plane is not included in the plotting, and it would be included after the new plotting funtion is ready
		plotter.FillDataHistogram(data);
		plotter.FillModelHistogram(averaged_sum_pdf);
		plotter.FillOtherHistogram("cmb_bkg", "background", combinatorial_background_pdf, f_cmb(), 16, 7, 38);
		if (outfilename != "") plotter.SaveHistograms(outfilename);
		plotter.SetCustomAxesTitles("#it{m}^{2}_{+} [GeV^{2}/#it{c}^{4}]","#it{m}^{2}_{#minus} [GeV^{2}/#it{c}^{4}]","#it{m}^{2}_{#it{#pi#pi}} [GeV^{2}/#it{c}^{4}]");

		// 1D Projection
		TCanvas c1("c1","c1",1800,700);
		TPad *pad1 = new TPad("pad1","pad1",0.01,0.25,0.33,0.99);
		TPad *pad2 = new TPad("pad2","pad2",0.01,0.01,0.33,0.25);
		TPad *pad3 = new TPad("pad3","pad3",0.33,0.25,0.65,0.99);
		TPad *pad4 = new TPad("pad4","pad4",0.33,0.01,0.65,0.25);
		TPad *pad5 = new TPad("pad5","pad5",0.65,0.25,0.97,0.99);
		TPad *pad6 = new TPad("pad6","pad6",0.65,0.01,0.97,0.25);
		pad1->Draw();
		pad2->Draw();
		pad3->Draw();
		pad4->Draw();
		pad5->Draw();
		pad6->Draw();
		pad1->SetLeftMargin(0.15);
		pad2->SetLeftMargin(0.15);
		pad3->SetLeftMargin(0.15);
		pad4->SetLeftMargin(0.15);
		pad5->SetLeftMargin(0.15);
		pad6->SetLeftMargin(0.15);

		pad1->cd();
		plotter.Plot1DProjections(0, 0);
		pad2->cd();
		plotter.Plot1DPull(0);

		pad3->cd();
		plotter.Plot1DProjections(1, 1); // plot legend in this pad
		pad4->cd();
		plotter.Plot1DPull(1);

		pad5->cd();
		plotter.Plot1DProjections(2, 0); 
		pad6->cd();
		plotter.Plot1DPull(2);

		outfilename = args.outdir + outprefix + "-1d-projection";
		Print::Canvas(c1, outfilename);

		// 2D Projection
		TCanvas c2("c2","c2",1500,1000);
		c2.Divide(3,2);

		c2.cd(1);
		plotter.Plot2DProjectionData(0,1);

		c2.cd(2);
		plotter.Plot2DProjectionModel(0,1);

		c2.cd(3);
		plotter.Plot2DPull(0,1);

		c2.cd(4);
		plotter.Plot2DProjectionData(1,2);

		c2.cd(5);
		plotter.Plot2DProjectionModel(1,2);

		c2.cd(6);
		plotter.Plot2DPull(1,2);

		outfilename = args.outdir + outprefix + "-2d-projection";
		Print::Canvas(c2,outfilename);

		// 2D Projection
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
			myapp.Run();
		}
	}
	
	return 0;
}

