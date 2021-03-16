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
#include <physics/Amplitudes.h>
using namespace dafne;

// Output-files prefix
std::string outprefix("fit-KSpipi-time-independent-with-background");

// Define the arguments of the amplitude
declarg(MSqZero  , double)
declarg(MSqPlus  , double)
declarg(MSqMinus , double)
declarg(SignalFaction, double)
// declarg(RandomBackfroundFaction, double)
declarg(CombinatorialBackfroundFaction, double)

using namespace hydra::arguments;

// Main
int main(int argc, char** argv)
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
	auto f_cmb = hydra::Parameter::Create("f_cmb").Value(0.5).Error(0.0001).Limits(-1.,1.).Fixed(0);

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

	// build averaged sum pdf for fitting with uniform background fraction

	auto uniform_f_cmb_functor = PassParameter(f_cmb);

	auto _build_sum_pdf = hydra::wrap_lambda(
			  [] __hydra_dual__ (hydra::tuple< double, double, double> input_functors){
			  		auto _pdf_sig = hydra::get<0>(input_functors);
			  		auto _f_cmb = hydra::get<1>(input_functors);
			  		auto _pdf_cmb = hydra::get<2>(input_functors);

			  		return  (1-_f_cmb)*_pdf_sig + _f_cmb*_pdf_cmb;
			  }
	);

	auto averaged_sum_pdf = hydra::compose(_build_sum_pdf, normalized_model, uniform_f_cmb_functor, combinatorial_background_pdf);

	// build pdf with event by event fraction
	auto f_cmb_passor = PassVariable<CombinatorialBackfroundFaction>();

	auto evt_frac_pdf = hydra::compose(_build_sum_pdf, normalized_model, f_cmb_passor, combinatorial_background_pdf);




	// Build hydra Pdf<FUNCTOR, INTEGRATOR>. Since the averaged_sum_pdf is already normalized, we only need to
	// add constant "1" as INTEGRATOR, to adapt to hydra framework's convention.
	// auto pdf = hydra::make_pdf( averaged_sum_pdf, ConstantIntegrator<hydra::device::sys_t>(1.0));
	auto pdf = hydra::make_pdf( evt_frac_pdf, ConstantIntegrator<hydra::device::sys_t>(1.0));
	std::cout << "Initial normalization for D0 PDF: "<< pdf.GetNorm() << " +/- " << pdf.GetNormError() << std::endl;


	//---------------------------------------------------------------------------------------
	// get input data from ROOT file
	//---------------------------------------------------------------------------------------
	std::cout << "***** Input data" << std::endl;
	std::cout << "Creating data containers ... ...  " << std::endl;

	// with uniform fraction
	//hydra::multivector<hydra::tuple<MSqPlus,MSqMinus,MSqZero>, hydra::device::sys_t> data; // create data container for D0 events

	// with event by event fraction
	hydra::multivector<hydra::tuple<MSqPlus,MSqMinus,MSqZero,SignalFaction,CombinatorialBackfroundFaction>, hydra::device::sys_t> data; // create data container for D0 events

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

		if (i % 100000 == 0) std::cout << "Reading " << i << "events" << std::endl;

		// data.push_back(hydra::make_tuple(MSqPlus(m2p),MSqMinus(m2m),MSqZero(m2z)));

		// with event by event fraction
		data.push_back(hydra::make_tuple(MSqPlus(m2p),MSqMinus(m2m),MSqZero(m2z),SignalFaction(f_sig_i),CombinatorialBackfroundFaction(f_cmb_i)));
	}
	
	file->Close();

	auto ncands = data.size();

	if (ncands < 1) {
		std::cout << "The events number in hydra container is less than 1, something must be wrong!" << std::endl;
		return -1;
	}
	std::cout << "Read " << ncands << " data candidates ... ... " << std::endl;

	std::cout << "Dataset test: " << std::endl;
	std::cout << data[0] << std::endl;
	std::cout << data[1] << std::endl;
	std::cout << data[2] << std::endl;

	//---------------------------------------------------------------------------------------
	// Build pdf and log-likelihood function from model
	//---------------------------------------------------------------------------------------
	std::cout << "***** Build likelihood" << std::endl;
	
	auto fcn = hydra::make_loglikehood_fcn( pdf, data.begin(), data.end() );
	MinuitTools::CheckFCN(fcn);
	
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

	// auto parameters = minimum.UserParameters();
	// auto covariance = minimum.UserCovariance();
	// std::cout << "***** Fit results:\n" << parameters << covariance << std::endl;
	std::cout << "***** Fit results:\n" << minimum.UserState() << std::endl;
	MinuitTools::CovarianceMatrixStatus(minimum);

	//---------------------------------------------------------------------------------------
	// Plot fit result
	//---------------------------------------------------------------------------------------
	if (args.plot) {
		// set and check the parameters
		std::cout << "***** Set model according to fit result" << std::endl;

		// parameters after fitting
		auto fitted_parameters = fcn.GetParameters().GetVariables();

		// build a dummy fcn to easily synchronize the paramters in plotting functor and fitting functor
		hydra::multivector<hydra::tuple<MSqPlus,MSqMinus>, hydra::device::sys_t> data_dummy;
		data_dummy.push_back(hydra::make_tuple(MSqPlus(1.0),MSqMinus(1.0)));
		auto averaged_sum_pdf_for_plotting = dafne::MinuitTools::UpdateParametersByBuildingDummyFCN(data_dummy, averaged_sum_pdf, minimum);
		auto amp_for_plotting_temp = dafne::MinuitTools::UpdateParametersByBuildingDummyFCN(data_dummy, rate(amp), minimum);
		auto amp_for_plotting = amp_for_plotting_temp.GetFunctor(hydra::placeholders::_1);

		// directly get the f_cmb for the uniform fraction case
		// double f_cmb_plotting = MinuitTools::GetParameterPointer(fitted_parameters, "f_cmb")->GetValue();

		// calculate the f_cmb for the event by event case
		double f_cmb_sum = 0;
		for( auto event : data ) {
			CombinatorialBackfroundFaction a = hydra::get<4>(event);
			f_cmb_sum += a;
		}		
		double f_cmb_plotting = f_cmb_sum / data.size();
		averaged_sum_pdf_for_plotting.SetParameter("f_cmb", f_cmb_plotting);

		// plot the model and its components
		std::cout << "***** Plot data and model" << std::endl;

		TApplication myapp("myapp",0,0);
		
		auto plotter = DalitzPlotter<MSqPlus, MSqMinus, MSqZero>(phsp,"#it{K}^{0}_{S}","#it{#pi}^{+}","#it{#pi}^{#minus}",(args.prlevel>3));
		
		std::string outfilename = args.outdir + outprefix + "-HIST.root";
		// uniform fraction case
		// auto data_for_plotting = data;
		// event by event fraction case
		auto data_for_plotting = data.column(hydra::placeholders::_0, hydra::placeholders::_1, hydra::placeholders::_2);
		plotter.FillDataHistogram(data_for_plotting);
		plotter.FillModelHistogram(averaged_sum_pdf_for_plotting);
		plotter.FillOtherHistogram("cmb_bkg", "background", combinatorial_background_pdf, f_cmb_plotting, 16, 7, 38);
		plotter.FillComponentsHistogramsWithEfficiency(amp_for_plotting, efficiency, 1. - f_cmb_plotting);
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

		// 2D Projection
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
			myapp.Run();
		}
	}

	return 0;
}

