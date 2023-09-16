#include <hydra/host/System.h>
#include <hydra/device/System.h>
#include <hydra/FunctorArithmetic.h>
#include <hydra/Algorithm.h>
#include <hydra/Pdf.h>
#include <hydra/LogLikelihoodFCN.h>
#include <hydra/functions/JohnsonSUShape.h>

#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnUserParameterState.h>
#include <Minuit2/MnPrint.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnContours.h>
#include <Minuit2/MnPlot.h>

#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <Math/QuantFuncMathCore.h>

#include <physics/PDG.h>
#include <physics/Rate.h>
#include <physics/ThreeBodyPhaseSpace.h>
#include <models/D0ToKsPiPi_FVECTOR_BABAR.h>
#include <tools/Plotter.h>
#include <tools/Printer.h>
#include <tools/MinuitTools.h>
#include <tools/ConfigFile.h>
#include <tools/Arguments.h>
#include <physics/Amplitudes.h>
#include <tools/ArbitraryBinningHistogram2D.h>
#include <tools/ArbitraryBinningHistogram2DFunctor.h>
using namespace dafne;

std::string outprefix("fit-KSpipi-time-dependent-time-resolution-with-background");

// define the arguments of the amplitude
declarg(MSqZero  , double)
declarg(MSqPlus  , double)
declarg(MSqMinus , double)
declarg(DecayTime, double)
declarg(DecayTimeError, double)
using namespace hydra::arguments;

// Main
int main(int argc, char** argv)
{
	//---------------------------------------------------------------------------------------
	// Read command line arguments
	//---------------------------------------------------------------------------------------
	D0ToKsPiPiArguments args("Fit to D0 -> K0_S pi+ pi- dataset with BABAR 2010 amplitude model","0.1","");
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
	hydra::multivector<hydra::tuple<MSqPlus,MSqMinus,MSqZero,DecayTime,DecayTimeError>, hydra::device::sys_t> data_dz; // create data container for D0 events
	hydra::multivector<hydra::tuple<MSqPlus,MSqMinus,MSqZero,DecayTime,DecayTimeError>, hydra::device::sys_t> data_db; // create data container for D0-bar events
	auto file = TFile::Open(args.input.c_str());

	if (!file) {
		std::cout << "Failed to open input file: " << args.input << std::endl;
		return -1;
	}

	std::cout << "Reading data from input file: " << args.input << std::endl;
	auto ntp = (TTree*) file->Get("ntp");

	auto nentries = ntp->GetEntries();
	if (nentries<1) {
		std::cout << "Empty input file!" << std::endl; 
		return -1;
	}

	double m2p, m2m, m2z, t, sigmat;
	int flavor;
	ntp->SetBranchAddress("mSqP",&m2p);
	ntp->SetBranchAddress("mSqM",&m2m);
	ntp->SetBranchAddress("mSqZ",&m2z);
	ntp->SetBranchAddress("t",&t);
	ntp->SetBranchAddress("sigmat",&sigmat);
	ntp->SetBranchAddress("flavor",&flavor);

	for (auto i=0; i<nentries; ++i) {
		ntp->GetEntry(i);

		if (i % 100000 == 0) std::cout << "Reading " << i << " events" << std::endl;

		if (flavor > 0) {
			data_dz.push_back(hydra::make_tuple(MSqPlus(m2p),MSqMinus(m2m),MSqZero(m2z),DecayTime(t),DecayTimeError(sigmat)));
		} else if (flavor < 0) {
			data_db.push_back(hydra::make_tuple(MSqPlus(m2p),MSqMinus(m2m),MSqZero(m2z),DecayTime(t),DecayTimeError(sigmat)));
		} else {
			std::cout << "entry " << i << " : wrong flavor value!" << std::endl; 
		}

	}
	
	file->Close();

	auto ncands_dz = data_dz.size();
	auto ncands_db = data_db.size();
	auto ncands = ncands_dz + ncands_db;
	if (ncands < 1) {
		std::cout << "The events number in hydra container is less than 1, some thing must be wrong!" << std::endl;
		return -1;
	}
	std::cout << "Read " << ncands_dz << " D0 data candidates and " << ncands_db << " D0-bar events " << std::endl;

	if (args.prlevel > 3) {
		std::cout << "D0 dataset test: " << std::endl;
		std::cout << data_dz[0] << std::endl;
		std::cout << data_dz[1] << std::endl;
		std::cout << data_dz[2] << std::endl;
		std::cout << "D0-bar data test: " << std::endl;
		std::cout << data_db[0] << std::endl;
		std::cout << data_db[1] << std::endl;
		std::cout << data_db[2] << std::endl;
	}



	//---------------------------------------------------------------------------------------
	// Build model for D0->KS pi+ pi- decays
	//---------------------------------------------------------------------------------------
	auto phsp = D0ToKsPiPi_FVECTOR_BABAR::PhaseSpaceWithTimeAndTimeError();

	// build baseline amplitudes model
	auto Adir = D0ToKsPiPi_FVECTOR_BABAR::Amplitude<MSqPlus,MSqMinus>(phsp);
	
	// config the model according to configuration file
	ConfigFile config(args.config_file.c_str(), (args.prlevel>3) );
	config.ConfigureModel(Adir);

	auto Abar = hydra::wrap_lambda( [&Adir] __hydra_dual__ (MSqPlus a, MSqMinus b) {
		MSqPlus switched_a = b;
		MSqMinus switched_b = a;
		return Adir(switched_a,switched_b);
	});

	
	// Time-dependent parameters
	auto tau = hydra::Parameter::Create("tau").Value(Tau::D0).Error(0.0001);
	auto x   = hydra::Parameter::Create("x").Value(0.003).Error(0.0001).Limits(-1.,1.);
	auto y   = hydra::Parameter::Create("y").Value(0.006).Error(0.0001).Limits(-1.,1.);
	auto qop = hydra::Parameter::Create("qop").Value(1.0).Error(0.0001).Limits(0.,10.);
	auto phi = hydra::Parameter::Create("phi").Value(0.0).Error(0.0001).Limits(-1.,1.);
	config.ConfigureMixingParameters(tau, x, y, qop, phi);
	

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

	// the time dependence of the efficiency is ignored 
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

	// time resolution
	auto b   = hydra::Parameter::Create("b").Value(0.0).Error(0.0001).Limits(-1.,1.);
	auto s   = hydra::Parameter::Create("s").Value(1.0).Error(0.0001).Limits(0.5,1.5);

	auto johnson_delta  = hydra::Parameter::Create().Name("johnson_delta").Value(1.65335e+00).Error(0.01);
	auto johnson_lambda = hydra::Parameter::Create().Name("johnson_lambda").Value(1.87922e-02).Error(0.001);
	auto johnson_gamma  = hydra::Parameter::Create().Name("johnson_gamma" ).Value(-2.57429e+00).Error(0.01);
	auto johnson_xi     = hydra::Parameter::Create().Name("johnson_xi").Value(4.27580e-02).Error(0.001);

	config.ConfigureTimeResolutionParameters({&b, &s, &johnson_delta, &johnson_lambda, &johnson_gamma, &johnson_xi});

	auto johnson_su = NormalizedJohnsonSU<DecayTimeError>(johnson_gamma, johnson_delta, johnson_xi, johnson_lambda, phsp.TimeErrorMin(), phsp.TimeErrorMax());
	
	auto model_dz = time_dependent_rate_with_time_resolution_pdf<Flavor::Positive,MSqPlus,MSqMinus,DecayTime,DecayTimeError>(
						tau,x,y,qop,phi,b,s,efficiency,Adir,Abar,johnson_su,
						{phsp.MSqMin<1,2>(),phsp.MSqMax<1,2>()},
						{phsp.MSqMin<1,3>(),phsp.MSqMax<1,3>()},
						phsp.TimeRange(),
						phsp.TimeErrorRange(), 100*ncands); 
	auto model_db = time_dependent_rate_with_time_resolution_pdf<Flavor::Negative,MSqPlus,MSqMinus,DecayTime,DecayTimeError>(
						tau,x,y,qop,phi,b,s,efficiency,Adir,Abar,johnson_su,
						{phsp.MSqMin<1,2>(),phsp.MSqMax<1,2>()},
						{phsp.MSqMin<1,3>(),phsp.MSqMax<1,3>()},
						phsp.TimeRange(),
						phsp.TimeErrorRange(), 100*ncands); 


	// add background

	auto f_rnd = hydra::Parameter::Create("f_rnd").Value(0.0035).Fixed();
	auto f_cmb = hydra::Parameter::Create("f_cmb").Value(0.0117).Fixed();
	auto tau_cmb = hydra::Parameter::Create("tau_cmb").Value(Tau::D0*0.9).Fixed();
	auto b_cmb   = hydra::Parameter::Create("b_cmb").Value(0.0).Fixed();
	auto s_cmb   = hydra::Parameter::Create("s_cmb").Value(1.0).Fixed();
	config.ConfigureBackgroundParameters({&f_rnd, &f_cmb, &tau_cmb, &b_cmb, &s_cmb});

	auto f_rnd_functor = PassParameter(f_rnd);
	auto random_background_pdf = hydra::wrap_lambda( 
		[model_dz, model_db] __hydra_dual__ (MSqPlus a, MSqMinus b, DecayTime t, DecayTimeError sigma_t) {
			MSqPlus switched_a = b;
			MSqMinus switched_b = a;
			return 0.5*model_dz(a,b,t,sigma_t) + 0.5*model_db(switched_a,switched_b,t,sigma_t);
	});

	
	auto f_cmb_functor = PassParameter(f_cmb);

	// load combinatorial background from a histrogram
	TH2D* th2_input_cmb = new TH2D("th2_input_cmb", "th2_input_cmb", 200, 0, 3.2, 200, 0, 3.2);
	for (int i_x = 0; i_x < th2_input_cmb->GetXaxis()->GetNbins(); ++i_x)
	for (int i_y = 0; i_y < th2_input_cmb->GetYaxis()->GetNbins(); ++i_y) {
		double m2m = th2_input_cmb->GetXaxis()->GetBinCenter(i_x);
		double m2p = th2_input_cmb->GetYaxis()->GetBinCenter(i_y);
		if (phsp.Contains<2,3>(m2p, m2m)) th2_input_cmb->SetBinContent(i_x, i_y, 1.0);
		else th2_input_cmb->SetBinContent(i_x, i_y, 0.0);
	}

	// normalize the background
	double combinatorial_background_without_time_norm = 0;
	for (int i_x = 0; i_x < th2_input_cmb->GetXaxis()->GetNbins(); ++i_x)
	for (int i_y = 0; i_y < th2_input_cmb->GetYaxis()->GetNbins(); ++i_y) {
		combinatorial_background_without_time_norm += th2_input_cmb->GetBinContent(i_x, i_y);
	}

	for (int i_x = 0; i_x < th2_input_cmb->GetXaxis()->GetNbins(); ++i_x)
	for (int i_y = 0; i_y < th2_input_cmb->GetYaxis()->GetNbins(); ++i_y) {
		double value = th2_input_cmb->GetBinContent(i_x, i_y);
		th2_input_cmb->SetBinContent(i_x, i_y, value/combinatorial_background_without_time_norm);
	}
	
	ArbitraryBinningHistogram2D hist_input_cmb(*th2_input_cmb);
	th2_input_cmb->Draw("COLZ");
	Print::Canvas(cefficiency,  args.outdir + outprefix + "_cmb_input_th2");

	TCanvas chist_input_cmb("chist_input_cmb", "chist_input_cmb", 800, 600);
	gStyle->SetOptStat(0);
	gPad->SetRightMargin(0.15);
	hist_input_cmb.GetTH2D((outprefix + "_cmb_input_hist").c_str(), 
	                    (outprefix + "_cmb_input_hist").c_str(),
	                    "m^{2}_{#it{-}} [GeV^{2}/#it{c}^{4}]", "m^{2}_{#it{+}} [GeV^{2}/#it{c}^{4}]")->Draw("COLZ");
	Print::Canvas(chist_input_cmb,  args.outdir + outprefix + "_cmb_input_hist");
	gStyle->SetOptStat(1);

	// auto combinatorial_background_without_time = hydra::wrap_lambda( 
	// 	[phsp, hist_input_cmb] __hydra_dual__ (MSqPlus m2p, MSqMinus m2m) {
	// 		// judge whether in phase space or not
	// 		if (!phsp.Contains<2,3>(m2p, m2m)) return 0.0;

	// 		double x = m2m;
	// 		double y = m2p;
	// 		return hist_input_cmb.GetValue(x, y);
	// 	}
	// );

	hydra::device::vector<double> combinatorial_background_device_buffer_xs;
	hydra::device::vector<double> combinatorial_background_device_buffer_ys;
	hydra::device::vector<double> combinatorial_background_device_buffer_zs;
	hydra::device::vector<double> combinatorial_background_device_buffer_zerrors;
	auto combinatorial_background_without_time = hydra::make_arbitrary_binning_histogram_2d_functor<MSqPlus, MSqMinus, double>(
		combinatorial_background_device_buffer_xs,
		combinatorial_background_device_buffer_ys,
		combinatorial_background_device_buffer_zs,
		combinatorial_background_device_buffer_zerrors,
		hist_input_cmb
		);


    double time_min = phsp.TimeMin();
    double time_max = phsp.TimeMax();
	auto combinatorial_psi = hydra::wrap_lambda( 
		[time_min, time_max] __hydra_dual__ (unsigned int npar, const hydra::Parameter* params, DecayTime t, DecayTimeError sigmat) {
			double _tau_cmb = params[0];
			double _s_cmb = params[1];
			double _b_cmb = params[2];

			double _t = t;
			double _sigmat = sigmat;
			double _sigma = _s_cmb*_sigmat;
			double chi = (_t-_b_cmb) / _sigma;
			double kappa_0 = _sigma / _tau_cmb;
			double chi0 = (time_min-_b_cmb) / _sigma;
			double chi1 = (time_max-_b_cmb) / _sigma;
			double time_dependent_component = _psi(chi, kappa_0); 
			double time_dependent_component_norm = _int_psi_dt(_sigma, chi1, kappa_0) - _int_psi_dt(_sigma, chi0, kappa_0);
			return time_dependent_component / time_dependent_component_norm;
	}, tau_cmb, s_cmb, b_cmb);


	auto combinatorial_background_pdf = combinatorial_background_without_time * combinatorial_psi * johnson_su;
	// auto combinatorial_background_pdf = hydra::wrap_lambda( 
	// 	[phsp, hist_input_cmb] __hydra_dual__ (unsigned int npar, const hydra::Parameter* params, MSqPlus m2p, MSqMinus m2m, DecayTime t, DecayTimeError sigmat) {
	// 		if (!phsp.Contains<2,3>(m2p, m2m)) return 0.0;

	// 		double x = m2m;
	// 		double y = m2p;
	// 		double time_independent_component = hist_input_cmb.GetValue(x, y);

	// 		double _tau_cmb = params[0];
	// 		double _s_cmb = params[1];
	// 		double _b_cmb = params[2];

	// 		double _t = t;
	// 		double _sigmat = sigmat;
	// 		double _sigma = _s_cmb*_sigmat;
	// 		double chi = (_t-_b_cmb) / _sigma;
	// 		double kappa_0 = _sigma / _tau_cmb;
	// 		double chi0 = (phsp.TimeMin()-_b_cmb) / _sigma;
	// 		double chi1 = (phsp.TimeMax()-_b_cmb) / _sigma;
	// 		double time_dependent_component = _psi(chi, kappa_0); 
	// 		double time_dependent_component_norm = _int_psi_dt(_sigma, chi1, kappa_0) - _int_psi_dt(_sigma, chi0, kappa_0);
	// 		time_dependent_component = time_dependent_component / time_dependent_component_norm;

	// 		// JohnsonSU
	// 		double _johnson_delta  = 1.65335e+00;
	// 		double _johnson_lambda = 1.87922e-02;
	// 		double _johnson_gamma  = -2.57429e+00;
	// 		double _johnson_xi     = 4.27580e-02;
	// 		double inverse_lambda = 1.0/_johnson_lambda;
	// 		double z0 = (phsp.TimeMin()-_johnson_xi)*inverse_lambda;
	// 		double C0 = _johnson_gamma + _johnson_delta*::asinh(z0);
	// 		double johnson_su_int_0 = 0.5*(1.0 + ::erf(C0*hydra::math_constants::inverse_sqrt2));
	// 		double z1 = (phsp.TimeMax()-_johnson_xi)*inverse_lambda;
	// 		double C1 = _johnson_gamma + _johnson_delta*::asinh(z1);
	// 		double johnson_su_int_1 = 0.5*(1.0 + ::erf(C1*hydra::math_constants::inverse_sqrt2));
	// 		double r = johnson_su_int_1 - johnson_su_int_0;

	// 		double z = (t-_johnson_xi)*inverse_lambda;
	// 		double A = inverse_lambda*_johnson_delta*hydra::math_constants::inverse_sqrt2Pi;
	// 		double B = 1.0/::sqrt( 1 + z*z);
	// 		double C = _johnson_gamma + _johnson_delta*::asinh(z); C *=C;

	// 		double sigmat_component = A*B*::exp(-0.5*C) / r;
	// 		return time_independent_component * time_dependent_component * sigmat_component;


	// }, tau_cmb, s_cmb, b_cmb);

	auto _build_averaged_sum_pdf = hydra::wrap_lambda(
			  [] __hydra_dual__ (hydra::tuple< double, double, double, double, double> input_functors){
			  		auto _pdf_sig = hydra::get<0>(input_functors);
			  		auto _f_rnd = hydra::get<1>(input_functors);
			  		auto _pdf_rnd = hydra::get<2>(input_functors);
			  		auto _f_cmb = hydra::get<3>(input_functors);
			  		auto _pdf_cmb = hydra::get<4>(input_functors);


			  		return  (1-_f_rnd-_f_cmb)*_pdf_sig + _f_rnd*_pdf_rnd + _f_cmb*_pdf_cmb;
			  }
	);

	auto averaged_sum_pdf_dz = hydra::compose(_build_averaged_sum_pdf, model_dz, 
											f_rnd_functor, random_background_pdf,
											f_cmb_functor, combinatorial_background_pdf);

	auto averaged_sum_pdf_db = hydra::compose(_build_averaged_sum_pdf, model_db, 
											f_rnd_functor, random_background_pdf,
											f_cmb_functor, combinatorial_background_pdf);



	// build pdf
	
	// the time_dependent_rate_with_time_resolution_pdf_type1 functor is a pdf itself, but the FCN needs 
	// Pdf<functor, integrator> as input, therefore we add a dummy constant integrator, always returning 
	// 1.0, here
	auto pdf_dz = hydra::make_pdf( averaged_sum_pdf_dz, ConstantIntegrator<hydra::device::sys_t>(1.0) ); 
	std::cout << "Initial normalization for D0 PDF: "<< pdf_dz.GetNorm() << " +/- " << pdf_dz.GetNormError() << std::endl;

	auto pdf_db = hydra::make_pdf( averaged_sum_pdf_db, ConstantIntegrator<hydra::device::sys_t>(1.0) );
	std::cout << "Initial normalization for D0bar PDF: "<< pdf_db.GetNorm() << " +/- " << pdf_db.GetNormError() << std::endl;





	//---------------------------------------------------------------------------------------
	// Build pdf and log-likelihood function from model
	//---------------------------------------------------------------------------------------
	std::cout << "***** Build likelihood" << std::endl;
	
	auto fcn_dz = hydra::make_loglikehood_fcn( pdf_dz, data_dz.begin(), data_dz.end() );
	auto fcn_db = hydra::make_loglikehood_fcn( pdf_db, data_db.begin(), data_db.end() );
	
	auto fcn = hydra::make_simultaneous_fcn(fcn_dz, fcn_db);
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
		std::cout << "***** Plot data and fit result" << std::endl;

		// plotting procedure
		TApplication* myapp = NULL;
		if (args.interactive) myapp = new TApplication("myapp",0,0);


		// data_dz + data_db will be plotted together with model_dz ignoring the CPV
		auto data = data_dz;
		data.insert(data.end(), data_db.begin(), data_db.end());


		// model for plotting
		std::cout << "build a plotting model for signal ....." << std::endl;

		// synchronize the paramters from the fitter to the plotting model
		auto model_truth_dz = time_dependent_rate<Flavor::Positive,DecayTime>(tau,x,y,qop,phi,Adir,Abar); 
		auto phsp_for_plotting = D0ToKsPiPi_FVECTOR_BABAR::PhaseSpaceWithTime();
		auto model_truth_dz_fitted = dafne::MinuitTools::UpdateParameters<MSqPlus, MSqMinus, MSqZero, DecayTime>(phsp_for_plotting, model_truth_dz, minimum);


		std::cout << "build a plotting model for random background ....." << std::endl;
		auto model_truth_db =  time_dependent_rate<Flavor::Negative,DecayTime>(tau,x,y,qop,phi,Adir,Abar);
		auto random_background_truth = hydra::wrap_lambda( 
		[&model_truth_dz, &model_truth_db] __hydra_dual__ (MSqPlus a, MSqMinus b, DecayTime t) {
			MSqPlus switched_a = b;
			MSqMinus switched_b = a;
			return 0.5*model_truth_dz(a,b,t) + 0.5*model_truth_db(switched_a,switched_b,t);
		});
		auto random_background_truth_fitted = dafne::MinuitTools::UpdateParameters<MSqPlus, MSqMinus, MSqZero, DecayTime>(phsp_for_plotting, random_background_truth, minimum);

		std::cout << "build a plotting model for combinatorial background ....." << std::endl;
		// get the fitted parameters
		fcn_dz.GetParameters().UpdateParameters(minimum);
		auto fitted_parameters = fcn_dz.GetParameters().GetVariables();
		double Gamma_cmb_fitted = 1. / MinuitTools::GetParameterPointer(fitted_parameters, "tau_cmb")->GetValue();
		auto exp_cmb_fitted = hydra::wrap_lambda( [Gamma_cmb_fitted] __hydra_dual__ (DecayTime t)
		{
			return std::exp(-Gamma_cmb_fitted * t)  ;
		});
		auto combinatorial_background_truth_fitted = combinatorial_background_without_time * exp_cmb_fitted;


		// plot Dalitz-plot distributions 
		auto plotter = DalitzPlotterWithTimeAndTimeError<MSqPlus, MSqMinus, MSqZero, DecayTime, DecayTimeError>(phsp,"#it{K}^{0}_{S}","#it{#pi}^{+}","#it{#pi}^{#minus}",(args.prlevel>3));

		std::string outfilename = args.outdir + outprefix + "-HIST.root";
		double tau_for_plot = MinuitTools::GetParameterPointer(fitted_parameters, "tau")->GetValue();
		double y_for_plot = MinuitTools::GetParameterPointer(fitted_parameters, "y")->GetValue();
		double b_for_plot = MinuitTools::GetParameterPointer(fitted_parameters, "b")->GetValue();
		double s_for_plot = MinuitTools::GetParameterPointer(fitted_parameters, "s")->GetValue();
		double f_rnd_for_plot = MinuitTools::GetParameterPointer(fitted_parameters, "f_rnd")->GetValue();
		double f_cmb_for_plot = MinuitTools::GetParameterPointer(fitted_parameters, "f_cmb")->GetValue();

		double tau_cmb_for_plot = MinuitTools::GetParameterPointer(fitted_parameters, "tau_cmb")->GetValue();
		double b_cmb_for_plot = MinuitTools::GetParameterPointer(fitted_parameters, "b_cmb")->GetValue();
		double s_cmb_for_plot = MinuitTools::GetParameterPointer(fitted_parameters, "s_cmb")->GetValue();

		THnSparseD* h_data =  plotter.FillDataHistogram(data);
		THnSparseD* h_sig = plotter.FillHistogram("signal", "Signal", model_truth_dz_fitted, efficiency, 
			tau_for_plot, y_for_plot, b_for_plot, s_for_plot, 
			johnson_su, 1-f_rnd_for_plot-f_cmb_for_plot,
			2, 1, 0);
		THnSparseD* h_rnd = plotter.FillHistogram("rnd_bkg", "Random #pi^{s}", random_background_truth_fitted, efficiency, 
			tau_for_plot, y_for_plot, b_for_plot, s_for_plot, 
			johnson_su, f_rnd_for_plot,
			28, 3, 6);

		auto dummy_efficiency = ConstantFunctor(1);
		THnSparseD* h_cmb = plotter.FillHistogram("cmb_bkg", "Combinatorial", combinatorial_background_truth_fitted, dummy_efficiency, 
			tau_cmb_for_plot, 0, b_cmb_for_plot, s_cmb_for_plot,  
			johnson_su, f_cmb_for_plot,  
			16, 7, 41);
		plotter.FillModelHistogramFromExternalHistograms({"signal", "rnd_bkg", "cmb_bkg"});
		plotter.SetCustomAxesTitles("#it{m}^{2}_{+} [GeV^{2}/#it{c}^{4}]","#it{m}^{2}_{#minus} [GeV^{2}/#it{c}^{4}]","#it{m}^{2}_{#it{#pi#pi}} [GeV^{2}/#it{c}^{4}]");
		if (outfilename != "") plotter.SaveHistograms(outfilename);
		


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
		TH1D * h1_sig = plotter.Plot1DProjection("signal", 0, "histo");
		h1_sig->SetLineWidth(1);
		TH1D * h1_cmb = plotter.Plot1DProjection("cmb_bkg", 0, "histo same");
		TH1D * h1_rnd = plotter.Plot1DProjection("rnd_bkg", 0, "histo same");

		THStack * h1_all = new THStack("h1_all", "h1_all");
		h1_all->Add(h1_cmb);
		h1_all->Add(h1_rnd);
		h1_all->Add(h1_sig);

		h1_all->Draw();
		h1_all->SetMaximum(h1_all->GetMaximum()*1.1);
		plotter.Plot1DProjectionData(0, "e1 same");


		pad2->cd();
		TH1D* h1_pull = plotter.Plot1DPull(0);
		plotter.PlotPullLines(h1_pull->GetXaxis()->GetXmin(), h1_pull->GetXaxis()->GetXmax());


		pad3->cd();
		h1_sig = plotter.Plot1DProjection("signal", 1, "histo");
		h1_sig->SetLineWidth(1);
		h1_cmb = plotter.Plot1DProjection("cmb_bkg", 1, "histo same");
		h1_rnd = plotter.Plot1DProjection("rnd_bkg", 1, "histo same");

		h1_all = new THStack("h1_all", "h1_all");
		h1_all->Add(h1_cmb);
		h1_all->Add(h1_rnd);
		h1_all->Add(h1_sig);

		h1_all->Draw();
		h1_all->SetMaximum(h1_all->GetMaximum()*1.1);
		plotter.Plot1DProjectionData(1, "e1 same");


		pad4->cd();
		h1_pull = plotter.Plot1DPull(1);
		plotter.PlotPullLines(h1_pull->GetXaxis()->GetXmin(), h1_pull->GetXaxis()->GetXmax());


		pad5->cd();
		h1_sig = plotter.Plot1DProjection("signal", 2, "histo");
		h1_sig->SetLineWidth(1);
		h1_cmb = plotter.Plot1DProjection("cmb_bkg", 2, "histo same");
		h1_rnd = plotter.Plot1DProjection("rnd_bkg", 2, "histo same");

		h1_all = new THStack("h1_all", "h1_all");
		h1_all->Add(h1_cmb);
		h1_all->Add(h1_rnd);
		h1_all->Add(h1_sig);

		h1_all->Draw();
		h1_all->SetMaximum(h1_all->GetMaximum()*1.1);
		TH1D* h1_data = plotter.Plot1DProjectionData(2, "e1 same");

		TLegend* leg = new TLegend(0.55,0.67,0.82,0.88);
		leg->SetBorderSize(0);
		leg->SetFillStyle(0);
		leg->AddEntry(h1_data, h_data->GetTitle(), "pe");
		leg->AddEntry(h1_sig, h_sig->GetTitle(), "l");
		leg->AddEntry(h1_rnd, h_rnd->GetTitle(), "lf");
		leg->AddEntry(h1_cmb, h_cmb->GetTitle(), "lf");
		leg->Draw();


		pad6->cd();
		h1_pull = plotter.Plot1DPull(2);
		plotter.PlotPullLines(h1_pull->GetXaxis()->GetXmin(), h1_pull->GetXaxis()->GetXmax());

		outfilename = args.outdir + outprefix + "-1d-projection";
		Print::Canvas(c1,outfilename);

		// 2D projections
		TCanvas c2("c2","c2",1500,500);
		c2.Divide(3,1);

		c2.cd(1);
		plotter.Plot2DProjectionData(0,1); 

		c2.cd(2);
		plotter.Plot2DProjectionData(1,2);

		c2.cd(3);
		plotter.Plot2DPull(0,1);

		c2.Update();
		outfilename = args.outdir + outprefix + "-2d-projection";
		Print::Canvas(c2,outfilename);
		
		// Plot phase difference
		TCanvas c3("c3","c3",600,500);
		c3.SetRightMargin(.14);

		auto phspWithoutTime = D0ToKsPiPi_FVECTOR_BABAR::PhaseSpace();
		auto plotterWithoutTime = DalitzPlotter<MSqPlus, MSqMinus, MSqZero>(phsp,"#it{K}^{0}_{S}","#it{#pi}^{+}","#it{#pi}^{#minus}",(args.prlevel>3));
		plotterWithoutTime.PlotPhaseDifference(Adir,Abar);
		
		outfilename = args.outdir + outprefix + "-phase-difference";
		Print::Canvas(c3,outfilename);

		// time distribution
		TCanvas c4("c4","c4",1200,500);
		TPad *pad7 = new TPad("pad7","pad7",0.01,0.25,0.49,0.99);
		TPad *pad8 = new TPad("pad8","pad8",0.01,0.01,0.49,0.25);
		TPad *pad9 = new TPad("pad9","pad9",0.5,0.25,0.99,0.99);
		TPad *pad10 = new TPad("pad10","pad10",0.5,0.01,0.99,0.25);	
		pad7->Draw();
		pad8->Draw();
		pad9->Draw();
		pad10->Draw();
		pad7->SetLeftMargin(0.15);
		pad8->SetLeftMargin(0.15);
		pad9->SetLeftMargin(0.15);
		pad10->SetLeftMargin(0.15);

		pad7->cd();
		h1_sig = plotter.Plot1DProjection("signal", 3, "histo");
		h1_sig->SetLineWidth(1);
		h1_cmb = plotter.Plot1DProjection("cmb_bkg", 3, "histo same");
		h1_rnd = plotter.Plot1DProjection("rnd_bkg", 3, "histo same");

		h1_all = new THStack("h1_all", "h1_all");
		h1_all->Add(h1_cmb);
		h1_all->Add(h1_rnd);
		h1_all->Add(h1_sig);

		h1_all->Draw();
		h1_all->SetMaximum(h1_all->GetMaximum()*1.1);
		h1_data = plotter.Plot1DProjectionData(3, "e1 same");


		pad8->cd();
		h1_pull = plotter.Plot1DPull(3);
		plotter.PlotPullLines(h1_pull->GetXaxis()->GetXmin(), h1_pull->GetXaxis()->GetXmax());


		pad9->cd();
		h1_all->Draw();
		h1_all->SetMaximum(h1_all->GetMaximum()*1.1);
		h1_data = plotter.Plot1DProjectionData(3, "e1 same");
		pad9->SetLogy();


		leg = new TLegend(0.6,0.67,0.8,0.88);
		leg->SetBorderSize(0);
		leg->SetFillStyle(0);
		leg->AddEntry(h1_data, h_data->GetTitle(), "pe");
		leg->AddEntry(h1_sig, h_sig->GetTitle(), "l");
		leg->AddEntry(h1_rnd, h_rnd->GetTitle(), "lf");
		leg->AddEntry(h1_cmb, h_cmb->GetTitle(), "lf");
		leg->Draw();

		pad10->cd();
		h1_pull->Draw();
		plotter.PlotPullLines(h1_pull->GetXaxis()->GetXmin(), h1_pull->GetXaxis()->GetXmax());

		outfilename = args.outdir + outprefix + "-decay-time";
		Print::Canvas(c4,outfilename);


		// sigmat distribution
		TCanvas c5("c5","c5",1200,500);
		pad7 = new TPad("pad7","pad7",0.01,0.25,0.49,0.99);
		pad8 = new TPad("pad8","pad8",0.01,0.01,0.49,0.25);
		pad9 = new TPad("pad9","pad9",0.5,0.25,0.99,0.99);
		pad10 = new TPad("pad10","pad10",0.5,0.01,0.99,0.25);	
		pad7->Draw();
		pad8->Draw();
		pad9->Draw();
		pad10->Draw();
		pad7->SetLeftMargin(0.15);
		pad8->SetLeftMargin(0.15);
		pad9->SetLeftMargin(0.15);
		pad10->SetLeftMargin(0.15);


		pad7->cd();
		h1_sig = plotter.Plot1DProjection("signal", 4, "histo");
		h1_sig->SetLineWidth(1);
		h1_cmb = plotter.Plot1DProjection("cmb_bkg", 4, "histo same");
		h1_rnd = plotter.Plot1DProjection("rnd_bkg", 4, "histo same");

		h1_all = new THStack("h1_all", "h1_all");
		h1_all->Add(h1_cmb);
		h1_all->Add(h1_rnd);
		h1_all->Add(h1_sig);

		h1_all->Draw();
		h1_all->SetMaximum(h1_all->GetMaximum()*1.1);
		h1_data = plotter.Plot1DProjectionData(4, "e1 same");


		pad8->cd();
		h1_pull = plotter.Plot1DPull(4);
		plotter.PlotPullLines(h1_pull->GetXaxis()->GetXmin(), h1_pull->GetXaxis()->GetXmax());


		pad9->cd();
		h1_all->Draw();
		h1_all->SetMaximum(h1_all->GetMaximum()*1.1);
		h1_data = plotter.Plot1DProjectionData(4, "e1 same");
		pad9->SetLogy();


		leg = new TLegend(0.6,0.67,0.8,0.88);
		leg->SetBorderSize(0);
		leg->SetFillStyle(0);
		leg->AddEntry(h1_data, h_data->GetTitle(), "pe");
		leg->AddEntry(h1_sig, h_sig->GetTitle(), "l");
		leg->AddEntry(h1_rnd, h_rnd->GetTitle(), "lf");
		leg->AddEntry(h1_cmb, h_cmb->GetTitle(), "lf");
		leg->Draw();

		pad10->cd();
		h1_pull->Draw();
		plotter.PlotPullLines(h1_pull->GetXaxis()->GetXmin(), h1_pull->GetXaxis()->GetXmax());

		outfilename = args.outdir + outprefix + "-sigmat";
		Print::Canvas(c5,outfilename);

		// plot contours
		if (args.contour_option != "") {
			TCanvas c5("c5","c5",600,500);
			c5.SetRightMargin(.14);

			// a function to cope with root plotting issue, the "CFA" in root does not work some time
			auto plotCFA = [](TGraph * tgin, const bool drawAxis=1, const char* xtitle="x [%]", const char* ytitle="y [%]") {
				TGraph* fgraph = new TGraph(*tgin);
				TGraph* cgraph = new TGraph(*tgin);

				if (drawAxis) {
					fgraph->Draw("FA");
					fgraph->GetXaxis()->SetTitle(xtitle);
					fgraph->GetYaxis()->SetTitle(ytitle);
				}
				else fgraph->Draw("F");
				cgraph->Draw("C");

			};
			
			if (args.contour_option == "test") {
				TGraph* contGraphSigma1 = MinuitTools::Contour(fcn, minimum, "cont_sigma1", "x", "y", 1.1478745, 15, 1, "silence"); // 1 sigma

				contGraphSigma1->SetLineColor(45);
				contGraphSigma1->SetLineWidth(3);
				contGraphSigma1->SetFillStyle(3244);
				contGraphSigma1->SetFillColor(46);

				plotCFA(contGraphSigma1);

				leg = new TLegend(0.75,0.9-0.06,0.9,0.9);
				leg->AddEntry(contGraphSigma1, " 1 #sigma ", "f");
			} else if (args.contour_option == "sigma1") {
				TGraph* contGraphSigma1 = MinuitTools::Contour(fcn, minimum, "cont_sigma1", "x", "y", 1.1478745, 30, 1, "silence"); // 1 sigma

				contGraphSigma1->SetLineColor(45);
				contGraphSigma1->SetLineWidth(3);
				contGraphSigma1->SetFillStyle(3244);
				contGraphSigma1->SetFillColor(46);

				plotCFA(contGraphSigma1);

				leg = new TLegend(0.75,0.9-0.06,0.9,0.9);
				leg->AddEntry(contGraphSigma1, " 1 #sigma ", "f");
			} else if (args.contour_option == "sigma3") {
				TGraph* contGraphSigma3 = MinuitTools::Contour(fcn, minimum, "cont_sigma3", "x", "y", 5.9145790, 30, 1, "silence"); // 3 sigma
				TGraph* contGraphSigma1 = MinuitTools::Contour(fcn, minimum, "cont_sigma1", "x", "y", 1.1478745, 30, 1, "silence"); // 1 sigma

				contGraphSigma3->SetLineColor(44);
				contGraphSigma3->SetLineWidth(3);
				contGraphSigma3->SetFillStyle(3001);
				contGraphSigma3->SetFillColor(42);
				contGraphSigma1->SetLineColor(45);
				contGraphSigma1->SetLineWidth(3);
				contGraphSigma1->SetFillStyle(3244);
				contGraphSigma1->SetFillColor(46);

				plotCFA(contGraphSigma3);
				plotCFA(contGraphSigma1, 0);

				leg = new TLegend(0.75,0.9-0.06*2,0.9,0.9);
				leg->AddEntry(contGraphSigma1, " 1 #sigma ", "f");
				leg->AddEntry(contGraphSigma3, " 3 #sigma ", "f");
			} else if (args.contour_option == "sigma5") {
				TGraph* contGraphSigma5 = MinuitTools::Contour(fcn, minimum, "cont_sigma5", "x", "y", 14.371851, 30, 1, "silence"); // 5 sigma
				TGraph* contGraphSigma3 = MinuitTools::Contour(fcn, minimum, "cont_sigma3", "x", "y", 5.9145790, 30, 1, "silence"); // 3 sigma
				TGraph* contGraphSigma1 = MinuitTools::Contour(fcn, minimum, "cont_sigma1", "x", "y", 1.1478745, 30, 1, "silence"); // 1 sigma


				// the color setting for contours is to be finished
				contGraphSigma3->SetLineColor(43);
				contGraphSigma3->SetFillStyle(1001);
				contGraphSigma3->SetFillColor(41);
				contGraphSigma3->SetLineColor(44);
				contGraphSigma3->SetFillStyle(3001);
				contGraphSigma3->SetFillColor(42);
				contGraphSigma1->SetLineColor(45);
				contGraphSigma1->SetFillStyle(3244);
				contGraphSigma1->SetFillColor(46);

				plotCFA(contGraphSigma5);
				plotCFA(contGraphSigma3, 0);
				plotCFA(contGraphSigma1, 0);

				leg = new TLegend(0.75,0.9-0.06*3,0.9,0.9);
				leg->AddEntry(contGraphSigma1, " 1 #sigma ", "f");
				leg->AddEntry(contGraphSigma3, " 3 #sigma ", "f");
				leg->AddEntry(contGraphSigma5, " 5 #sigma ", "f");
			} else if (args.contour_option == "sigmafull") {
				TGraph* contGraphSigma5 = MinuitTools::Contour(fcn, minimum, "cont_sigma5", "x", "y", 14.371851, 30, 1, "silence"); // 5 sigma
				TGraph* contGraphSigma4 = MinuitTools::Contour(fcn, minimum, "cont_sigma4", "x", "y", 9.6669543, 30, 1, "silence"); // 4 sigma
				TGraph* contGraphSigma3 = MinuitTools::Contour(fcn, minimum, "cont_sigma3", "x", "y", 5.9145790, 30, 1, "silence"); // 3 sigma
				TGraph* contGraphSigma2 = MinuitTools::Contour(fcn, minimum, "cont_sigma2", "x", "y", 3.0900372, 30, 1, "silence"); // 2 sigma
				TGraph* contGraphSigma1 = MinuitTools::Contour(fcn, minimum, "cont_sigma1", "x", "y", 1.1478745, 30, 1, "silence"); // 1 sigma

				// the color setting for contours is to be finished

				plotCFA(contGraphSigma5);
				plotCFA(contGraphSigma4, 0);
				plotCFA(contGraphSigma3, 0);
				plotCFA(contGraphSigma2, 0);
				plotCFA(contGraphSigma1, 0);

				leg = new TLegend(0.75,0.9-0.06*5,0.9,0.9);
				leg->AddEntry(contGraphSigma1, " 1 #sigma ", "f");
				leg->AddEntry(contGraphSigma2, " 2 #sigma ", "f");
				leg->AddEntry(contGraphSigma3, " 3 #sigma ", "f");
				leg->AddEntry(contGraphSigma4, " 4 #sigma ", "f");
				leg->AddEntry(contGraphSigma5, " 5 #sigma ", "f");
			}



			leg->SetBorderSize(0);
			leg->SetFillStyle(0);
			leg->Draw();

			

			outfilename = args.outdir + outprefix + "-contour";
			Print::Canvas(c5,outfilename);
		}

		
		if (args.interactive) {
			std::cout << "Press Crtl+C to terminate" << std::endl;
			myapp->Run();
		}
		

		// plot decay-time distribution


	}


	return 0;
}

