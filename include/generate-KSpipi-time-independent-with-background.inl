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
declarg(SignalFaction, double)
declarg(CombinatorialBackgroundFaction, double)

using namespace hydra::arguments;

// output prefix
std::string outprefix = "generate-KSpipi-time-independent-with-background";

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
	auto model_truth = rate(amp);
	auto model = model_truth*efficiency;
	auto plainIntegrator = phsp.Integrator(10*args.nevents);
	auto normalized_model = make_numerical_normalized_functor(model, plainIntegrator);

	// add combinatorial background

	// define combinatorial background
	auto combinatorial_background = hydra::wrap_lambda( 
		[phsp] __hydra_dual__ (MSqPlus m2p, MSqMinus m2m) {
			// judge whether in phase space or not
			if (!phsp.Contains<2,3>(m2p, m2m)) return 0.0;
			return 1.;
		}
	);
	// normalize the background and build a background pdf
	auto combinatorial_background_norm = plainIntegrator(combinatorial_background).first;
	auto combinatorial_background_pdf =  divideBy<double>(combinatorial_background, combinatorial_background_norm); 

	// build a sum pdf for the plotting. The f_cmb in the pdf below is the average of the 
	// the f_cmb of all events.
	auto _build_sum_pdf = hydra::wrap_lambda(
			  [] __hydra_dual__ (hydra::tuple< double, double, double> input_functors){
			  		auto _pdf_sig = hydra::get<0>(input_functors);
			  		auto _f_cmb = hydra::get<1>(input_functors);
			  		auto _pdf_cmb = hydra::get<2>(input_functors);

			  		return  (1-_f_cmb)*_pdf_sig + _f_cmb*_pdf_cmb;
			  }
	);
	auto f_cmb = hydra::Parameter::Create("f_cmb").Value(0.15).Error(0.01).Limits(0,1.);
	auto f_cmb_functor = PassParameter(f_cmb);
	auto averaged_sum_pdf = hydra::compose(_build_sum_pdf, normalized_model, f_cmb_functor, combinatorial_background_pdf);

	// validate the pdf combination, by checking the values in a certain point
	if (args.prlevel > 3) {
		std::cout << "combinatorial_background_pdf(1., 1.): ";
		std::cout << combinatorial_background_pdf(1., 1.) << std::endl ; 

		std::cout << "normalized_model(1., 1.): ";
		std::cout << normalized_model(1., 1.) << std::endl ;

		std::cout << "(1-f_cmb) * normalized_model(1., 1.) + f_cmb * combinatorial_background_pdf(1., 1.): ";
		std::cout << (1-f_cmb.GetValue()) * normalized_model(MSqPlus(1.), MSqMinus(1.)) 
					 + f_cmb.GetValue() * combinatorial_background_pdf(MSqPlus(1.), MSqMinus(1.)) << std::endl;

		std::cout << "averaged_sum_pdf(1., 1.): ";
		std::cout << averaged_sum_pdf(MSqPlus(1.), MSqMinus(1.)) << std::endl ;
	}


	//---------------------------------------------------------------------------------------
	// Generate data
	//---------------------------------------------------------------------------------------
	std::cout << "***** Generation" << std::endl;

	hydra::multivector<hydra::tuple<MSqPlus,MSqMinus,MSqZero>, hydra::device::sys_t> data;

	// define the expected average combinatorial fraction
	double f_cmb_expected = 0.15;

	std::cout << "Generating signal ... ... " << std::endl;
	auto start = std::chrono::high_resolution_clock::now();

	auto data_sig = phsp.GenerateData<MSqPlus,MSqMinus,MSqZero>(model, (1-f_cmb_expected)*args.nevents, args.seed);
	auto ngenerated = data_sig.size();
	
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> elapsed = end - start;
	std::cout << "Generated " << ngenerated << " signal events." << std::endl;
	std::cout << "Time elapsed (ms): "<< elapsed.count() << std::endl;


	std::cout << "Generating combinatorial background ... ... " << std::endl;
	start = std::chrono::high_resolution_clock::now();

	auto data_cmb = phsp.GenerateData<MSqPlus,MSqMinus,MSqZero>(combinatorial_background, f_cmb_expected*args.nevents, args.seed);
	ngenerated = data_cmb.size();
	
	end = std::chrono::high_resolution_clock::now();
	elapsed = end - start;
	std::cout << "Generated " << ngenerated << " combinatorial background events." << std::endl;
	std::cout << "Time elapsed (ms): "<< elapsed.count() << std::endl;

	// merge and shuffle
	data.insert( data.end(), data_sig.begin(), data_sig.end());
	data.insert( data.end(), data_cmb.begin(), data_cmb.end());
	std::random_shuffle(data.begin(), data.end());
	ngenerated = data.size();

	std::cout << "Generated " << ngenerated << " events in total." << std::endl;


	// events by events signal and background fractions
	double f_sig_expected = 1.0 - f_cmb_expected;
	hydra::multivector< hydra::tuple<SignalFaction,CombinatorialBackgroundFaction>, hydra::device::sys_t> fractions;
	for (int i = 0; i < ngenerated; ++i) {
		SignalFaction _f_sig = f_sig_expected;
		CombinatorialBackgroundFaction _f_cmb = f_cmb_expected;
		fractions.push_back(hydra::make_tuple(_f_sig, _f_cmb));
	}

	// The following code samples background fraction from a Gaussian distribution
	// double N_gaussian = 2*f_sig_expected/erf(1./sqrt(2.)); // suppose the fractions are from a gaussian distribution with
	//                                                    // x in [-1sigma, 1sigma]
	// auto gaussian_for_f_cmb = hydra::wrap_lambda( [N_gaussian] __hydra_dual__ (double _x)
	// {
	// 		double sigma = 1;
	// 		double mean = 0;
	// 		double x = (_x-mean)/sigma;
	// 		double x2 = x*x;

	// 		const double one_div_sqrt_2_pi = 0.39894228;

	// 		double _f_sig = N_gaussian * one_div_sqrt_2_pi * std::exp(-x2/2.);
	// 		double _f_cmb = 1.0 - _f_sig;

	// 		SignalFaction _f_sig_output = _f_sig;
	// 		CombinatorialBackgroundFaction _f_cmb_output = _f_cmb;

	// 		return hydra::make_tuple(_f_sig, _f_cmb);
	// });
	// 
	// auto uniform = hydra::UniformShape<double>(-1, 1);
	// hydra::device::vector<double> xs_for_fractions(ngenerated);
	// hydra::fill_random(xs_for_fractions, uniform, args.seed);
	// hydra::multivector< hydra::tuple<SignalFaction,CombinatorialBackgroundFaction>, hydra::device::sys_t> fractions = xs_for_fractions | gaussian_for_f_cmb;
	// 
	// some sampling process to match the background fraction value and the background population
	// under this value ... ... 

	auto data_with_fraction = data.meld(fractions);




	//---------------------------------------------------------------------------------------
	// Save to ROOT file
	//---------------------------------------------------------------------------------------

	// prepare the resolution for smearing
	Resolution2D resolution2D;
	config.ConfigureResolution2D(resolution2D);
	TCanvas cresolution("cresolution", "cresolution", 800, 600);
	resolution2D.Draw("m^{2}_{#it{#pi#pi}} [GeV^{2}/#it{c}^{4}]", "cos(#theta_{#it{#pi#pi}})");
	Print::Canvas(cresolution,  args.outdir + outprefix + "_resolution");

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

		// branches for event by event fractions
		double f_sig, f_cmb;
		ntp->Branch("f_sig",&f_sig);
		ntp->Branch("f_cmb",&f_cmb);		

		// save to file
		for( auto event : data_with_fraction )
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

			// add event by event fraction
			SignalFaction e = hydra::get<3>(event);
			CombinatorialBackgroundFaction f = hydra::get<4>(event);
			f_sig = e.Value();
			f_cmb = f.Value();	

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
		
		// calculate the f_cmb for the event by event case
		double f_cmb_plotting = 0;
		for( auto event : data_with_fraction ) {
			CombinatorialBackgroundFaction a = hydra::get<4>(event);
			f_cmb_plotting += a;
		}		
		f_cmb_plotting = f_cmb_plotting / data.size();
		std::cout << "The calculated f_cmb_plotting = " << f_cmb_plotting << std::endl;
		averaged_sum_pdf.SetParameter("f_cmb", f_cmb_plotting);

		// plot the Dalitz-plot distributions		
		auto plotter = DalitzPlotter<MSqPlus, MSqMinus, MSqZero>(phsp,"#it{K}^{0}_{S}","#it{#pi}^{+}","#it{#pi}^{#minus}",(args.prlevel>3));
		
		plotter.FillDataHistogram(data);
		plotter.FillModelHistogram(averaged_sum_pdf);
		plotter.FillHistogram("cmb_bkg", "background", combinatorial_background_pdf, f_cmb(), 16, 7, 38);
		plotter.FillComponentHistograms(amp, efficiency, 1. - f_cmb());
		plotter.SetCustomAxesTitles("#it{m}^{2}_{+} [GeV^{2}/#it{c}^{4}]","#it{m}^{2}_{#minus} [GeV^{2}/#it{c}^{4}]","#it{m}^{2}_{#it{#pi#pi}} [GeV^{2}/#it{c}^{4}]");

		std::string outfilename = args.outdir + outprefix + "-HIST.root";
		plotter.SaveHistograms(outfilename);

		// 1D projections
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
		auto h1_pull = plotter.Plot1DPull(0);
		plotter.PlotPullLines(h1_pull->GetXaxis()->GetXmin(), h1_pull->GetXaxis()->GetXmax());


		pad3->cd();
		plotter.Plot1DProjections(1, 1); // plot legend in this pad
		pad4->cd();
		h1_pull = plotter.Plot1DPull(1);
		plotter.PlotPullLines(h1_pull->GetXaxis()->GetXmin(), h1_pull->GetXaxis()->GetXmax());

		pad5->cd();
		plotter.Plot1DProjections(2, 0); 
		pad6->cd();
		h1_pull = plotter.Plot1DPull(2);
		plotter.PlotPullLines(h1_pull->GetXaxis()->GetXmin(), h1_pull->GetXaxis()->GetXmax());
		outfilename = args.outdir + outprefix + "-1d-projection";
		Print::Canvas(c1, outfilename);

		// 2D projections
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

