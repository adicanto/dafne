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
#include <tools/ConfigFile.h>
#include <tools/Arguments.h>
#include <physics/Amplitudes.h>
#include <tools/ArbitraryBinningHistogram2D.h>
using namespace dafne;

// Configuration parameters
const unsigned nevents(100000);
std::string outprefix("fit-KSpipi-time-dependent");

// Define the arguments of the amplitude
declarg(MSqZero  , double)
declarg(MSqPlus  , double)
declarg(MSqMinus , double)
declarg(DecayTime, double)
using namespace hydra::arguments;

// Main
int main(int argc, char** argv)
{
	//---------------------------------------------------------------------------------------
	// Read command line arguments
	//---------------------------------------------------------------------------------------
	Arguments args("Fit to D0 -> K0_S pi+ pi- dataset with BABAR 2010 amplitude model","0.1","");
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
	auto phsp = D0ToKsPiPi_FVECTOR_BABAR::PhaseSpaceWithTime();

	// build baseline amplitudes model
	auto Adir = D0ToKsPiPi_FVECTOR_BABAR::Amplitude<MSqPlus,MSqMinus>(phsp);
	
	// config the model according to configuration file
	ConfigFile config(args.config_file.c_str(), (args.prlevel>3) );
	config.ConfigureModel(Adir);
	
	// Time-dependent parameters
	auto tau = hydra::Parameter::Create("tau").Value(Tau::D0).Fixed();
	auto x   = hydra::Parameter::Create("x").Value(0.003).Error(0.0001).Limits(-1.,1.);
	auto y   = hydra::Parameter::Create("y").Value(0.006).Error(0.0001).Limits(-1.,1.);
	auto qop = hydra::Parameter::Create("qop").Value(1.0).Error(0.0001).Limits(0.,10.);
	auto phi = hydra::Parameter::Create("phi").Value(0.0).Error(0.0001).Limits(-1.,1.);
	config.ConfigureMixingParameters(tau, x, y, qop, phi);
	
	// Total amplitude for the decay that undergoes mixing
	auto Abar = hydra::wrap_lambda( [&Adir] __hydra_dual__ (MSqPlus a, MSqMinus b) {
		MSqPlus switched_a = b;
		MSqMinus switched_b = a;
		return Adir(switched_a,switched_b);
	});

	// efficiency plane described by irregular binning 2D histogram
	ArbitraryBinningHistogram2D efficiency_hist = config.ConfigureEfficiencyHistogram();

	// move to the plotting part after the plotting part is ready
	TCanvas cefficiency("cefficiency", "cefficiency", 800, 600);
	gStyle->SetOptStat(0);
	gPad->SetRightMargin(0.15);
	efficiency_hist.Draw((outprefix + "_efficiency_hist").c_str(), 100);
	Print::Canvas(cefficiency,  args.outdir + outprefix + "efficiency_hist");
	gStyle->SetOptStat(1);

	// time dependent efficiency is ignored for the moment
	auto efficiency = hydra::wrap_lambda(
		[phsp, efficiency_hist] __hydra_dual__ (DecayTime tau, MSqPlus m2p, MSqMinus m2m) {

		// judge whether in phase space or not
		if (!phsp.Contains<2,3>(m2p, m2m)) return 0.0;

		double m2z = phsp.MSqJK(m2p, m2m);
		double helicity_z = phsp.HelicityJK<2,3>(m2p, m2m);

		double x = m2z;
		double y = helicity_z;
		return efficiency_hist.GetValue(x, y);

	}); 

	// D0 rate
	auto model_dz = time_dependent_rate<Flavor::Positive,DecayTime>(tau,x,y,qop,phi,Adir,Abar)*efficiency;
	
	// D0bar rate
	auto model_db = time_dependent_rate<Flavor::Negative,DecayTime>(tau,x,y,qop,phi,Adir,Abar)*efficiency;

	//---------------------------------------------------------------------------------------
	// get input data from ROOT file
	//---------------------------------------------------------------------------------------
	std::cout << "***** Input data" << std::endl;
	std::cout << "Creating data containers ... ...  " << std::endl;
	hydra::multivector<hydra::tuple<MSqPlus,MSqMinus,MSqZero,DecayTime>, hydra::device::sys_t> data_dz; // create data container for D0 events
	hydra::multivector<hydra::tuple<MSqPlus,MSqMinus,MSqZero,DecayTime>, hydra::device::sys_t> data_db; // create data container for D0-bar events
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

	double m2p, m2m, m2z, t;
	int flavor;
	ntp->SetBranchAddress("mSqP",&m2p);
	ntp->SetBranchAddress("mSqM",&m2m);
	ntp->SetBranchAddress("mSqZ",&m2z);
	ntp->SetBranchAddress("t",&t);
	ntp->SetBranchAddress("flavor",&flavor);

	for (auto i=0; i<nentries; ++i) {
		ntp->GetEntry(i);

		if (i % 100000 == 0) std::cout << "Reading " << i << "events" << std::endl;

		if (flavor > 0) {
			data_dz.push_back(hydra::make_tuple(MSqPlus(m2p),MSqMinus(m2m),MSqZero(m2z),DecayTime(t)));
		} else if (flavor < 0) {
			data_db.push_back(hydra::make_tuple(MSqPlus(m2p),MSqMinus(m2m),MSqZero(m2z),DecayTime(t)));
		} else {
			std::cout << "entry " << i << " : wrong flavor value!" << std::endl; 
		}

	}
	
	file->Close();

	auto ncands_dz = data_dz.size();
	auto ncands_db = data_db.size();
	if (ncands_dz + ncands_db < 1) {
		std::cout << "The events number in hydra container is less than 1, some thing must be wrong!" << std::endl;
		return -1;
	}
	std::cout << "Read " << ncands_dz << " D0 data candidates and " << ncands_db << " D0-bar events " << std::endl;

	std::cout << "D0 dataset test: " << std::endl;
	std::cout << data_dz[0] << std::endl;
	std::cout << data_dz[1] << std::endl;
	std::cout << data_dz[2] << std::endl;
	std::cout << "D0-bar data test: " << std::endl;
	std::cout << data_db[0] << std::endl;
	std::cout << data_db[1] << std::endl;
	std::cout << data_db[2] << std::endl;


	//---------------------------------------------------------------------------------------
	// Build pdf and log-likelihood function from model
	//---------------------------------------------------------------------------------------
	std::cout << "***** Build likelihood" << std::endl;
	auto pdf_dz = hydra::make_pdf( model_dz, phsp.IntegratorWithTime(10*ncands_dz) );
	std::cout << "Initial normalization for D0 PDF: "<< pdf_dz.GetNorm() << " +/- " << pdf_dz.GetNormError() << std::endl;
	
	auto pdf_db = hydra::make_pdf( model_db, phsp.IntegratorWithTime(10*ncands_db) );
	std::cout << "Initial normalization for D0bar PDF: "<< pdf_db.GetNorm() << " +/- " << pdf_db.GetNormError() << std::endl;
	
	auto fcn_dz = hydra::make_loglikehood_fcn( pdf_dz, data_dz.begin(), data_dz.end() );
	auto fcn_db = hydra::make_loglikehood_fcn( pdf_db, data_db.begin(), data_db.end() );
	
	auto fcn = hydra::make_simultaneous_fcn(fcn_dz, fcn_db);
	Print::CheckFCN(fcn);
	
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
	Print::MinimizerStatus(minimum);

	if ( !minimum.IsValid() || !minimum.HasAccurateCovar() ) {
		std::cout << "Fit did not converge or covariance matrix is not accurate." << std::endl;
		return -2;
	}

	// auto parameters = minimum.UserParameters();
	// auto covariance = minimum.UserCovariance();
	// std::cout << "***** Fit results:\n" << parameters << covariance << std::endl;
	std::cout << "***** Fit results:\n" << minimum.UserState() << std::endl;
	Print::CovarianceMatrixStatus(minimum);

	//---------------------------------------------------------------------------------------
	// Plot fit result
	//---------------------------------------------------------------------------------------
	if (args.plot) {
		std::cout << "***** Plot data and fit result" << std::endl;

		// get the fitted model_dz
		fcn_dz.GetParameters().UpdateParameters(minimum);
		auto model_dz_fitted = fcn_dz.GetPDF().GetFunctor().GetFunctor(hydra::placeholders::_0);
		model_dz_fitted.PrintRegisteredParameters();

		// data_dz + data_db are plotted with model_dz ignoring the CPV
		auto data = data_dz;
		data.insert(data.end(), data_db.begin(), data_db.end());

		// plotting procedure
		TApplication myapp("myapp",0,0);

		// plot dalitz distribution 
		auto plotterWithTime = DalitzPlotterWithTime<MSqPlus, MSqMinus, MSqZero, DecayTime>(phsp,"#it{K}^{0}_{S}","#it{#pi}^{+}","#it{#pi}^{#minus}",(args.prlevel>3));

		std::string outfilename = args.outdir + outprefix + "-HIST.root";
		plotterWithTime.FillHistograms(data, model_dz_fitted, outfilename, args.plotnbins); 
		plotterWithTime.SetCustomAxesTitles("#it{m}^{2}_{+} [GeV^{2}/#it{c}^{4}]","#it{m}^{2}_{#minus} [GeV^{2}/#it{c}^{4}]","#it{m}^{2}_{#it{#pi#pi}} [GeV^{2}/#it{c}^{4}]");

		// 1D Projection for dalitz distribution
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
		plotterWithTime.Plot1DProjectionData(0, "e1");
		plotterWithTime.Plot1DProjectionModel(0, "histo same");
		pad2->cd();
		plotterWithTime.Plot1DPull(0);

		pad3->cd();
		plotterWithTime.Plot1DProjectionData(1, "e1");
		plotterWithTime.Plot1DProjectionModel(1, "histo same");
		pad4->cd();
		plotterWithTime.Plot1DPull(1);

		pad5->cd();
		TH1D* h1_data = plotterWithTime.Plot1DProjectionData(2, "e1");
		TH1D* h1_model = plotterWithTime.Plot1DProjectionModel(2, "histo same");
		TLegend* leg = new TLegend(0.6,0.7,0.8,0.85);
		leg->SetBorderSize(0);
		leg->SetFillStyle(0);
		leg->AddEntry(h1_data, "Data", "pe");
		leg->AddEntry(h1_model, "Model Sum", "l");
		leg->Draw();

		pad6->cd();
		plotterWithTime.Plot1DPull(2);

		outfilename = args.outdir + outprefix + "-1d-projection";
		Print::Canvas(c1,outfilename);

		// 2D Projection for dalitz distribution
		TCanvas c2("c2","c2",1500,500);
		c2.Divide(3,1);

		c2.cd(1);
		plotterWithTime.Plot2DProjectionData(0,1); 
		// TH2D* h2_data = plotterWithTime.Plot2DProjectionData(0,1); // some effort to solve the -2d-projection.pdf plotting issue
		// h2_data->Draw("colz");

		c2.cd(2);
		plotterWithTime.Plot2DProjectionData(1,2);
		// h2_data = plotterWithTime.Plot2DProjectionData(1,2);
		// h2_data->Draw("colz");

		c2.cd(3);
		plotterWithTime.Plot2DPull(0,1);

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
		h1_data = plotterWithTime.Plot1DProjectionData(3, "e1");
		h1_model = plotterWithTime.Plot1DProjectionModel(3, "histo same");

		pad8->cd();
		TH1D* h1_pull = plotterWithTime.Plot1DPull(3);

		pad9->cd();
		h1_data->Draw("e1");
		h1_model->Draw("histo same");
		pad9->SetLogy();

		leg = new TLegend(0.75,0.75,0.9,0.9);
		leg->SetBorderSize(0);
		leg->SetFillStyle(0);
		leg->AddEntry(h1_data, h1_data->GetTitle(), "pe");
		leg->AddEntry(h1_model, h1_model->GetTitle(), "l");

		pad10->cd();
		h1_pull->Draw();

		outfilename = args.outdir + outprefix + "-decay-time";
		Print::Canvas(c4,outfilename);

		// plot contour
		ROOT::Minuit2::MnContours contours(fcn, minimum);
		double up = ROOT::Math::chisquared_quantile(0.6827,2)/2.;//1sigma
		fcn.SetErrorDef(up);
		std::vector<std::pair<double,double> > cont1 = contours(minimum.UserState().Index("x"), minimum.UserState().Index("y"), 30);

		TCanvas c5("c5","c5",600,500);
		c3.SetRightMargin(.14);

		auto ContourToTGraph = [](std::vector<std::pair<double,double> > contourPoints, const char * name="contour")->TGraph* {
			std::vector<double> xs(contourPoints.size()+1);
			std::vector<double> ys(contourPoints.size()+1);

			for (size_t i = 0; i < contourPoints.size(); ++i) {
				xs[i] = contourPoints[i].first;
				ys[i] = contourPoints[i].second;
			}
			xs[contourPoints.size()] = xs[0];
			ys[contourPoints.size()] = ys[0];

			TGraph* tg = new TGraph(contourPoints.size()+1, xs.data(), ys.data());
			tg->SetName(name);
			return tg;
		};

		auto cont1graph = ContourToTGraph(cont1, "cont1");
		// auto cont2graph = ContourToTGraph(cont2, "cont2");
		// auto cont3graph = ContourToTGraph(cont3, "cont3");
		// auto cont4graph = ContourToTGraph(cont4, "cont4");
		// auto cont5graph = ContourToTGraph(cont5, "cont5");

		cont1graph->Draw("al");
		// cont2graph->Draw("al");
		// cont3graph->Draw("al");
		// cont4graph->Draw("al");
		// cont5graph->Draw("al");

		outfilename = args.outdir + outprefix + "-contour";
		Print::Canvas(c5,outfilename);

		
		if (args.interactive) {
			std::cout << "Press Crtl+C to terminate" << std::endl;
			myapp.Run();
		}
		

		// plot decay-time distribution


	}


	return 0;
}

