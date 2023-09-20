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
#include <physics/Rate.h>
#include <physics/ThreeBodyPhaseSpace.h>
#include <models/D0ToKsPiPi_FVECTOR_BABAR.h>
#include <tools/Plotter.h>
#include <tools/Printer.h>
#include <tools/MinuitTools.h>
#include <tools/ConfigFile.h>
#include <tools/Arguments.h>
#include <tools/ArbitraryBinningHistogram2D.h>
#include <tools/ArbitraryBinningHistogram2DFunctor.h>
#include <physics/Amplitudes.h>
using namespace dafne;

// Output-files prefix
std::string outprefix("fit-KSpipi-time-independent-with-background");

// Define the arguments of the amplitude
declarg(MSqZero  , double)
declarg(MSqPlus  , double)
declarg(MSqMinus , double)
declarg(SignalFraction, double)
declarg(CombinatorialBackgroundFraction, double)

using namespace hydra::arguments;

// Main
int main(int argc, char** argv)
{
	//---------------------------------------------------------------------------------------
	// Read command line arguments
	//---------------------------------------------------------------------------------------
	StandardArguments args(
		"Generation of D0 -> K0_S pi+ pi- decays according to BABAR 2010 amplitude model","0.1","");

	try {
		args.Read( argc, argv );
	}
	catch ( TCLAP::ArgException& e ) {
		std::cerr << "Error: " << e.error() << " for argument " << e.argId() << "." << std::endl;
		return -1;
	}
	args.Print();



	//---------------------------------------------------------------------------------------
	// Get input data from the ROOT file
	//---------------------------------------------------------------------------------------
	std::cout << "***** Input data" << std::endl;
	std::cout << "Creating data containers ... ...  " << std::endl;

	// build a data container with event by event fraction
	hydra::multivector<hydra::tuple<MSqPlus,MSqMinus,MSqZero,SignalFraction,CombinatorialBackgroundFraction>, hydra::device::sys_t> data; // create data container for D0 events

	auto file = TFile::Open(args.input.c_str());
	if (!file) {
		std::cout << "Failed to open input file: " << args.input << std::endl;
		return -1;
	}

	std::cout << "Reading data from input file: " << args.input << std::endl;
	auto ntp = (TTree*) file->Get("ntp");

	auto nentries = ntp->GetEntriesFast();
	if (nentries<1) {
		std::cout << "Empty input file!" << std::endl;
		return -1;
	}

	double m2p, m2m, m2z;
	ntp->SetBranchAddress("mSqP",&m2p);
	ntp->SetBranchAddress("mSqM",&m2m);
	ntp->SetBranchAddress("mSqZ",&m2z);

	// add event by event fraction
	double f_sig_i, f_cmb_i;
	if (!ntp->GetBranch("f_sig") && !ntp->GetBranch("f_cmb")) {
		std::cout << "It seems f_sig and f_cmb are not in input file. Quit!" << std::endl;
		exit(-1);
	}
	ntp->SetBranchAddress("f_sig",&f_sig_i);
	ntp->SetBranchAddress("f_cmb",&f_cmb_i);		

	for (auto i=0; i<nentries; ++i) {
		ntp->GetEntry(i);

		if (i % 100000 == 0) std::cout << "Reading " << i << " events" << std::endl;

		// add event by event fraction to the container
		// in the fit, we only use f_cmb
		data.push_back(hydra::make_tuple(MSqPlus(m2p),MSqMinus(m2m),MSqZero(m2z),SignalFraction(f_sig_i),CombinatorialBackgroundFraction(f_cmb_i))); 
	}
	
	file->Close();

	auto ncands = data.size();

	if (ncands < 1) {
		std::cout << "The events number in hydra container is less than 1, something must be wrong!" << std::endl;
		return -1;
	}
	std::cout << "Read " << ncands << " data candidates ... ... " << std::endl;

	if (args.prlevel > 3) {
		std::cout << "Dataset test: " << std::endl;
		std::cout << data[0] << std::endl;
		std::cout << data[1] << std::endl;
		std::cout << data[2] << std::endl;
	}



	//---------------------------------------------------------------------------------------
	// Build model for D0->KS pi+ pi- decays
	//---------------------------------------------------------------------------------------
	auto phsp = D0ToKsPiPi_FVECTOR_BABAR::PhaseSpace();

	// build baseline amplitudes model
	auto amp = D0ToKsPiPi_FVECTOR_BABAR::Amplitude<MSqPlus,MSqMinus>(phsp);
	
	// config the model according to configuration file
	ConfigFile config(args.config_file.c_str(), (args.prlevel>3) );
	config.ConfigureModel(amp); // config the amplitude

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
	auto plainIntegrator = phsp.Integrator(10*ncands);
	auto normalized_model = make_numerical_normalized_functor(model, plainIntegrator);


	// add background
	auto combinatorial_background = hydra::wrap_lambda( 
		[phsp] __hydra_dual__ (MSqPlus m2p, MSqMinus m2m) {
			// judge whether in phase space or not
			if (!phsp.Contains<2,3>(m2p, m2m)) return 0.0;
			return 1.;
		}
	);
	
	auto combinatorial_background_norm = plainIntegrator(combinatorial_background).first;
	auto combinatorial_background_pdf = divideBy<double>(combinatorial_background, combinatorial_background_norm); 


	// load combinatorial background from a histrogram
	// TH2D* th2_input_cmb = new TH2D("th2_input_cmb", "th2_input_cmb", 200, 0, 3.2, 200, 0, 3.2);
	// for (int i_x = 0; i_x < th2_input_cmb->GetXaxis()->GetNbins(); ++i_x)
	// for (int i_y = 0; i_y < th2_input_cmb->GetYaxis()->GetNbins(); ++i_y) {
	// 	double m2m = th2_input_cmb->GetXaxis()->GetBinCenter(i_x);
	// 	double m2p = th2_input_cmb->GetYaxis()->GetBinCenter(i_y);
	// 	if (phsp.Contains<2,3>(m2p, m2m)) th2_input_cmb->SetBinContent(i_x, i_y, 1.0);
	// 	else th2_input_cmb->SetBinContent(i_x, i_y, 0.0);
	// }

	// // normalize the background
	// double combinatorial_background_norm = 0;
	// for (int i_x = 0; i_x < th2_input_cmb->GetXaxis()->GetNbins(); ++i_x)
	// for (int i_y = 0; i_y < th2_input_cmb->GetYaxis()->GetNbins(); ++i_y) {
	// 	combinatorial_background_norm += th2_input_cmb->GetBinContent(i_x, i_y);
	// }

	// for (int i_x = 0; i_x < th2_input_cmb->GetXaxis()->GetNbins(); ++i_x)
	// for (int i_y = 0; i_y < th2_input_cmb->GetYaxis()->GetNbins(); ++i_y) {
	// 	double value = th2_input_cmb->GetBinContent(i_x, i_y);
	// 	double widthx = th2_input_cmb->GetXaxis()->GetBinWidth(1);
	// 	double widthy = th2_input_cmb->GetYaxis()->GetBinWidth(1);
	// 	th2_input_cmb->SetBinContent(i_x, i_y, value/combinatorial_background_norm/widthx/widthy);
	// }
	
	// ArbitraryBinningHistogram2D hist_input_cmb(*th2_input_cmb);
	// th2_input_cmb->Draw("COLZ");
	// Print::Canvas(cefficiency,  args.outdir + outprefix + "_cmb_input_th2");

	// TCanvas chist_input_cmb("chist_input_cmb", "chist_input_cmb", 800, 600);
	// gStyle->SetOptStat(0);
	// gPad->SetRightMargin(0.15);
	// hist_input_cmb.GetTH2D((outprefix + "_cmb_input_hist").c_str(), 
	//                     (outprefix + "_cmb_input_hist").c_str(),
	//                     "m^{2}_{#it{-}} [GeV^{2}/#it{c}^{4}]", "m^{2}_{#it{+}} [GeV^{2}/#it{c}^{4}]")->Draw("COLZ");
	// Print::Canvas(chist_input_cmb,  args.outdir + outprefix + "_cmb_input_hist");
	// gStyle->SetOptStat(1);



	// hydra::device::vector<double> combinatorial_background_device_buffer_xs;
	// hydra::device::vector<double> combinatorial_background_device_buffer_ys;
	// hydra::device::vector<double> combinatorial_background_device_buffer_zs;
	// hydra::device::vector<double> combinatorial_background_device_buffer_zerrors;
	// auto combinatorial_background_pdf_test = hydra::make_arbitrary_binning_histogram_2d_functor<MSqPlus, MSqMinus, double>(
	// 	combinatorial_background_device_buffer_xs,
	// 	combinatorial_background_device_buffer_ys,
	// 	combinatorial_background_device_buffer_zs,
	// 	combinatorial_background_device_buffer_zerrors,
	// 	hist_input_cmb
	// 	);



	// build pdf with event by event fraction 
	auto f_cmb_passor = PassVariable<CombinatorialBackgroundFraction>();
	auto evt_frac_pdf = make_pdf_with_background_functor( model, plainIntegrator,
														  f_cmb_passor, combinatorial_background_pdf);



	auto pdf = hydra::make_pdf_from_normalized_functor( evt_frac_pdf );
	std::cout << "Initial normalization for D0 PDF: "<< pdf.GetNorm() << " +/- " << pdf.GetNormError() << std::endl;


	//---------------------------------------------------------------------------------------
	// Build pdf and log-likelihood function from model
	//---------------------------------------------------------------------------------------
	std::cout << "***** Build likelihood" << std::endl;
	
	auto fcn = hydra::make_loglikehood_fcn(pdf, data.begin(), data.end() );
	MinuitTools::CheckFCN(fcn);
	
	//---------------------------------------------------------------------------------------
	// Configure and run MINUIT
	//---------------------------------------------------------------------------------------
	std::cout << "***** Fit" << std::endl;
	// ROOT::Minuit2::MnPrint::SetLevel(2); // for earlier versions of CERN ROOT (before 6.24)
	ROOT::Minuit2::MnPrint::SetGlobalLevel(2); 
	ROOT::Minuit2::MnStrategy strategy(2);
	
	// print starting values of parameters
	auto pstart = fcn.GetParameters().GetMnState();
	std::cout << "Starting parameters:" << pstart << std::endl;
		
	ROOT::Minuit2::MnMigrad migrad(fcn, pstart, strategy);

	auto start = std::chrono::high_resolution_clock::now();
	auto minimum = ROOT::Minuit2::FunctionMinimum( migrad(5000., 1.) );
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> elapsed = end - start;
	
	std::cout << "-----------------------------------------" << std::endl;
	std::cout << "| [Migrad] Time (ms) ="<< elapsed.count() << std::endl;
	std::cout << "-----------------------------------------" << std::endl;


	std::cout << minimum << std::endl;
	MinuitTools::MinimizerStatus(minimum);

	if ( !minimum.IsValid() || !minimum.HasAccurateCovar() ) {
		std::cout << "Fit did not converge or covariance matrix is not accurate." << std::endl;
		return -2;
	}

	std::cout << "***** Fit results:\n" << minimum.UserState() << std::endl;
	MinuitTools::CovarianceMatrixStatus(minimum);

	//---------------------------------------------------------------------------------------
	// Plot fit result
	//---------------------------------------------------------------------------------------
	if (args.plot) {
		TApplication* myapp = NULL;
		if (args.interactive) myapp = new TApplication("myapp",0,0);

		std::cout << "***** Build the fitted model" << std::endl;

		// Build a fitted model according to the fit result
		auto normalized_model_fitted = dafne::MinuitTools::UpdateParameters<MSqPlus, MSqMinus, MSqZero>(phsp, normalized_model, minimum);

		// preparing a fitted sum amplitude, in order to plot the component
		auto amp_rate_fitted = dafne::MinuitTools::UpdateParameters<MSqPlus, MSqMinus, MSqZero>(phsp, rate(amp), minimum);
		auto amp_fitted = amp_rate_fitted.GetFunctor(hydra::placeholders::_1);

		// calculate the average f_cmb, for adding background in to the fitted model
		double f_cmb_average = 0;
		for( auto event : data ) {
			CombinatorialBackgroundFraction a = hydra::get<4>(event);
			f_cmb_average += a;
		}		
		f_cmb_average = f_cmb_average / data.size();
		std::cout << "calculated f_cmb_average = " << f_cmb_average << std::endl;

		auto f_cmb_for_plotting = hydra::Parameter::Create("f_cmb").Value(f_cmb_average).Error(0.01).Limits(0,1.).Fixed();
		auto uniform_f_cmb_functor_for_plotting = PassParameter(f_cmb_for_plotting);

		// Add the background to the fitted model, the background fraction is set to the average value
		// auto averaged_sum_pdf_fitted = hydra::compose(_build_sum_pdf, normalized_model_fitted, uniform_f_cmb_functor_for_plotting, combinatorial_background_pdf);
		auto averaged_sum_pdf_fitted = make_pdf_with_background_functor( model, plainIntegrator,
														  uniform_f_cmb_functor_for_plotting, combinatorial_background_pdf);

		std::cout << "***** Plot data and the fitted model" << std::endl;	

		auto plotter = DalitzPlotter<MSqPlus, MSqMinus, MSqZero>(phsp,"#it{K}^{0}_{S}","#it{#pi}^{+}","#it{#pi}^{#minus}",(args.prlevel>3));

		// prepare the dataset to plot

		// event by event fraction case, remove the dimensions of the fractions or weights
		auto data_for_plotting = data.column(hydra::placeholders::_0, hydra::placeholders::_1, hydra::placeholders::_2);

		// prepare the data and model histograms to plot
		plotter.FillDataHistogram(data_for_plotting);
		plotter.FillModelHistogram(averaged_sum_pdf_fitted);
		plotter.FillHistogram("cmb_bkg", "background", combinatorial_background_pdf, f_cmb_for_plotting.GetValue(), 16, 7, 38);
		plotter.FillComponentHistograms(amp_fitted, efficiency, 1. - f_cmb_for_plotting.GetValue());
		plotter.SetCustomAxesTitles("#it{m}^{2}_{+} [GeV^{2}/#it{c}^{4}]","#it{m}^{2}_{#minus} [GeV^{2}/#it{c}^{4}]","#it{m}^{2}_{#it{#pi#pi}} [GeV^{2}/#it{c}^{4}]");

		std::string outfilename = args.outdir + outprefix + "-HIST.root";	
		plotter.SaveHistograms(outfilename);
		

		// plot 
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

		// 1D projections
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
		TCanvas c2("c2","c2",1500,500);
		c2.Divide(3,1);

		c2.cd(1);
		plotter.Plot2DProjectionData(0,1);

		c2.cd(2);
		plotter.Plot2DProjectionData(1,2);

		c2.cd(3);
		plotter.Plot2DPull(0,1);

		outfilename = args.outdir + outprefix + "-2d-projection";
		Print::Canvas(c2, outfilename);

		// Phase difference
		TCanvas c3("c3","c3",600,500);
		c3.SetRightMargin(.14);
		
		auto abar = hydra::wrap_lambda( [&amp] __hydra_dual__ (MSqPlus a, MSqMinus b) {
			MSqPlus switched_a = b;
			MSqMinus switched_b = a;
			return amp(switched_a,switched_b);
		});
		
		plotter.PlotPhaseDifference(amp,abar);
		
		outfilename = args.outdir + outprefix + "-phase-difference";
		Print::Canvas(c3, outfilename);

		if (args.interactive) {
			std::cout << "Press Crtl+C to terminate" << std::endl;
			myapp->Run();
		}
	}

	return 0;
}

