#include <hydra/host/System.h>
#include <hydra/device/System.h>
#include <hydra/Placeholders.h>
#include <hydra/functions/JohnsonSUShape.h>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TStyle.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TLegend.h>

#include <iostream>
#include <iomanip>
#include <chrono>

#include <physics/PDG.h>
#include <physics/Rate.h>
#include <physics/ThreeBodyPhaseSpace.h>
#include <models/D0ToKsPiPi_FVECTOR_BABAR.h>
#include <tools/Plotter.h>
#include <tools/Printer.h>
#include <tools/ConfigFile.h>
#include <tools/Arguments.h>
#include <tools/ArbitraryBinningHistogram2D.h>


using namespace dafne;


// Define the arguments of the amplitude
declarg(MSqPlus , double)
declarg(MSqMinus , double)
declarg(MSqZero , double)
declarg(DecayTime, double)
declarg(DecayTimeError, double)
using namespace hydra::arguments;


// output prefix
std::string outprefix = "generate-KSpipi-time-dependent-time-resolution";


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
	auto phsp = D0ToKsPiPi_FVECTOR_BABAR::PhaseSpaceWithTimeAndTimeError();

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

 	// build and check the efficiency plane
	TCanvas cefficiency("cefficiency", "cefficiency", 800, 600);
	gStyle->SetOptStat(0);
	gPad->SetRightMargin(0.15);
	efficiency_hist.GetTH2D((outprefix + "_efficiency_hist").c_str(), 
	                    (outprefix + "_efficiency_hist").c_str(),
	                    "m^{2}_{#it{#pi#pi}} [GeV^{2}/#it{c}^{4}]", "cos(#theta_{#it{#pi#pi}})")->Draw("COLZ");
	Print::Canvas(cefficiency,  args.outdir + outprefix + "_efficiency_hist");
	gStyle->SetOptStat(1);

	// time dependent efficiency is ignored for the moment
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

	// time resolution test
	auto b   = hydra::Parameter::Create("b").Value(0.0).Error(0.0001).Limits(-1.,1.);
	auto s   = hydra::Parameter::Create("s").Value(1.0).Error(0.0001).Limits(0.9,1.1);

	auto johnson_delta  = hydra::Parameter::Create().Name("johnson_delta" ).Value(1.65335e+00).Error(0.01);
	auto johnson_lambda = hydra::Parameter::Create().Name("johnson_lambda").Value(1.87922e-02).Error(0.001);
	auto johnson_gamma  = hydra::Parameter::Create().Name("johnson_gamma" ).Value(-2.57429e+00).Error(0.01);
	auto johnson_xi     = hydra::Parameter::Create().Name("johnson_xi").Value(4.27580e-02).Error(0.001);

	config.ConfigureTimeResolutionParameters({&b, &s, &johnson_delta, &johnson_lambda, &johnson_gamma, &johnson_xi});

	auto johnson_su = hydra::JohnsonSU<DecayTimeError>(johnson_gamma, johnson_delta, johnson_xi, johnson_lambda);


	// auto model_dz = time_dependent_rate_with_time_resolution<Flavor::Positive,MSqPlus,MSqMinus,DecayTime,DecayTimeError>(tau,x,y,qop,phi,b,s,Adir,Abar,johnson_su)*efficiency; 
	// auto model_db = time_dependent_rate_with_time_resolution<Flavor::Negative,MSqPlus,MSqMinus,DecayTime,DecayTimeError>(tau,x,y,qop,phi,b,s,Adir,Abar,johnson_su)*efficiency;

	// "turth level" model without sigma_t, but ... ... includes efficiency plane
	auto model_truth_dz = time_dependent_rate<Flavor::Positive,DecayTime>(tau,x,y,qop,phi,Adir,Abar)*efficiency; 
	auto model_truth_db = time_dependent_rate<Flavor::Negative,DecayTime>(tau,x,y,qop,phi,Adir,Abar)*efficiency;

	// for checking the parameters order when debugging
	// typename decltype(model_dz)::argument_type  test{};
	// std::cout << test.dummy << '\n';

	//---------------------------------------------------------------------------------------
	// Generate data
	//---------------------------------------------------------------------------------------
	std::cout << "***** Generation" << std::endl;

	auto start = std::chrono::high_resolution_clock::now();	

	// auto data_dz = phsp.GenerateDataWithTimeAndTimeError<MSqPlus,MSqMinus,MSqZero,DecayTime,DecayTimeError>(model_dz,args.nevents,args.seed,0.,0.5,0.06,0.09,(args.prlevel>3));
	auto data_dz = phsp.GenerateDataWithTimeAndTimeError<MSqPlus,MSqMinus,MSqZero,DecayTime,DecayTimeError>(model_truth_dz,tau(),y(),b(),s(),johnson_su,args.nevents,args.seed,(args.prlevel>3));
	std::cout << "Generated " << data_dz.size() << " D0 candidates." << std::endl;
	// auto data_db = phsp.GenerateDataWithTimeAndTimeError<MSqPlus,MSqMinus,MSqZero,DecayTime,DecayTimeError>(model_db,args.nevents,args.seed+1,0.,0.5,0.06,0.09,(args.prlevel>3));
	auto data_db = phsp.GenerateDataWithTimeAndTimeError<MSqPlus,MSqMinus,MSqZero,DecayTime,DecayTimeError>(model_truth_db,tau(),y(),b(),s(),johnson_su,args.nevents,args.seed+1,(args.prlevel>3));
	std::cout << "Generated " << data_db.size() << " D0bar candidates." << std::endl;

	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> elapsed = end - start;
	std::cout << "Time elapsed (ms):"<< elapsed.count() << std::endl;
	
	//---------------------------------------------------------------------------------------
	// Save to ROOT file
	//---------------------------------------------------------------------------------------
	{
		std::string outfilename = args.outdir + outprefix + "-data.root";
		std::cout << "Saving output to " << outfilename << std::endl;
		TFile *outfile = new TFile(outfilename.c_str(),"recreate");

		double m2p, m2m, m2z, t, sigmat;
		int flavor(+1);
		TTree *ntp = new TTree("ntp","ntp");
		ntp->Branch("mSqP",&m2p);
		ntp->Branch("mSqM",&m2m);
		ntp->Branch("mSqZ",&m2z);
		ntp->Branch("t",&t);
		ntp->Branch("sigmat",&sigmat);
		ntp->Branch("flavor",&flavor);
		
		for( auto event : data_dz )
		{
			MSqPlus  	   a = hydra::get<0>(event);
			MSqMinus 	   b = hydra::get<1>(event);
			MSqZero  	   c = hydra::get<2>(event);
			DecayTime 	   d = hydra::get<3>(event);
			DecayTimeError f = hydra::get<4>(event);
			
			m2p 		= a.Value();
			m2m 		= b.Value();
			m2z 		= c.Value();
			t     		= d.Value();
			sigmat      = f.Value();
			
			ntp->Fill();
		}

		flavor = -1;
		for( auto event : data_db )
		{
			MSqPlus  	   a = hydra::get<0>(event);
			MSqMinus 	   b = hydra::get<1>(event);
			MSqZero 	   c = hydra::get<2>(event);
			DecayTime 	   d = hydra::get<3>(event);
			DecayTimeError f = hydra::get<4>(event);
			
			m2p 		= a.Value();
			m2m		    = b.Value();
			m2z		    = c.Value();
			t   	    = d.Value();
			sigmat      = f.Value();
			
			ntp->Fill();
		}
		
		outfile->Write();
		outfile->Close();
	}

	//---------------------------------------------------------------------------------------
	// Plot the model and toy MC
	//---------------------------------------------------------------------------------------
	if (args.plot) {
		std::cout << "***** Plot data and model" << std::endl;

		// data_dz + data_db are plotted with model_dz ignoring the CPV
		auto data = data_dz;
		data.insert(data.end(), data_db.begin(), data_db.end());

		// plot dalitz distribution 
		auto plotter = DalitzPlotterWithTimeAndTimeError<MSqPlus, MSqMinus, MSqZero, DecayTime, DecayTimeError>(phsp,"#it{K}^{0}_{S}","#it{#pi}^{+}","#it{#pi}^{#minus}",(args.prlevel>3));

		std::string outfilename = args.outdir + outprefix + "-HIST.root";
		// plotter.FillHistograms(data, model_dz, tau, y, b, s, outfilename, args.plotnbins); 
		plotter.FillDataHistogram(data, args.plotnbins);
		plotter.FillModelHistogramFast1(model_truth_dz, tau(), y(), b(), s(), johnson_su, args.plotnbins);
		if (outfilename != "") plotter.SaveHistograms(outfilename);
		plotter.SetCustomAxesTitles("#it{m}^{2}_{+} [GeV^{2}/#it{c}^{4}]","#it{m}^{2}_{#minus} [GeV^{2}/#it{c}^{4}]","#it{m}^{2}_{#it{#pi#pi}} [GeV^{2}/#it{c}^{4}]");

		// plotting procedure
		TApplication myapp("myapp",0,0);

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
		plotter.Plot1DProjectionData(0, "e1");
		plotter.Plot1DProjectionModel(0, "histo same");
		pad2->cd();
		TH1D* h1_pull = plotter.Plot1DPull(0);
		plotter.PlotPullLines(h1_pull->GetXaxis()->GetXmin(), h1_pull->GetXaxis()->GetXmax());

		pad3->cd();
		plotter.Plot1DProjectionData(1, "e1");
		plotter.Plot1DProjectionModel(1, "histo same");
		pad4->cd();
		h1_pull = plotter.Plot1DPull(1);
		plotter.PlotPullLines(h1_pull->GetXaxis()->GetXmin(), h1_pull->GetXaxis()->GetXmax());


		pad5->cd();
		TH1D* h1_data = plotter.Plot1DProjectionData(2, "e1");
		TH1D* h1_model = plotter.Plot1DProjectionModel(2, "histo same");
		TLegend* leg = new TLegend(0.6,0.7,0.8,0.85);
		leg->SetBorderSize(0);
		leg->SetFillStyle(0);
		leg->AddEntry(h1_data, "Data", "pe");
		leg->AddEntry(h1_model, "Model Sum", "l");
		leg->Draw();

		pad6->cd();
		h1_pull = plotter.Plot1DPull(2);
		plotter.PlotPullLines(h1_pull->GetXaxis()->GetXmin(), h1_pull->GetXaxis()->GetXmax());

		outfilename = args.outdir + outprefix + "-1d-projection";
		Print::Canvas(c1,outfilename);

		// 2D Projection for dalitz distribution
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
		h1_data = plotter.Plot1DProjectionData(3, "e1");
		h1_model = plotter.Plot1DProjectionModel(3, "histo same");

		pad8->cd();
		h1_pull = plotter.Plot1DPull(3);
		plotter.PlotPullLines(h1_pull->GetXaxis()->GetXmin(), h1_pull->GetXaxis()->GetXmax());

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
		h1_data = plotter.Plot1DProjectionData(4, "e1");
		h1_model = plotter.Plot1DProjectionModel(4, "histo same");

		pad8->cd();
		h1_pull = plotter.Plot1DPull(4);
		plotter.PlotPullLines(h1_pull->GetXaxis()->GetXmin(), h1_pull->GetXaxis()->GetXmax());

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
		plotter.PlotPullLines(h1_pull->GetXaxis()->GetXmin(), h1_pull->GetXaxis()->GetXmax());

		outfilename = args.outdir + outprefix + "-sigmat";
		Print::Canvas(c5,outfilename);
		
		if (args.interactive) {
			std::cout << "Press Crtl+C to terminate" << std::endl;
			myapp.Run();
		}
		
	}
	
	return 0;
}

