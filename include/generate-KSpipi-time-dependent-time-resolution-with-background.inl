#include <hydra/host/System.h>
#include <hydra/device/System.h>
#include <hydra/Placeholders.h>
#include <hydra/functions/JohnsonSUShape.h>
#include <hydra/Function.h>
#include <hydra/FunctorArithmetic.h>

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
#include <physics/Resolution.h>
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
std::string outprefix = "generate-KSpipi-time-dependent-time-resolution-with-background";


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
	auto model_truth_dz = time_dependent_rate<Flavor::Positive,DecayTime>(tau,x,y,qop,phi,Adir,Abar); 
	auto model_truth_db = time_dependent_rate<Flavor::Negative,DecayTime>(tau,x,y,qop,phi,Adir,Abar);

	// for checking the parameters order when debugging
	// typename decltype(model_dz)::argument_type  test{};
	// std::cout << test.dummy << '\n';

	// add background

	auto f_rnd = hydra::Parameter::Create("f_rnd").Value(0.25).Error(0.0001).Limits(0.,1.);
	auto tau_rnd = hydra::Parameter::Create("tau_rnd").Value(Tau::D0).Fixed();
	auto b_rnd   = hydra::Parameter::Create("b_rnd").Value(0.0).Fixed();
	auto s_rnd   = hydra::Parameter::Create("s_rnd").Value(1.0).Fixed();

	auto f_cmb = hydra::Parameter::Create("f_cmb").Value(0.25).Error(0.0001).Limits(0.,1.);
	auto tau_cmb = hydra::Parameter::Create("tau_cmb").Value(Tau::D0*0.9).Fixed();
	auto b_cmb   = hydra::Parameter::Create("b_cmb").Value(0.0).Fixed();
	auto s_cmb   = hydra::Parameter::Create("s_cmb").Value(1.0).Fixed();


	auto random_background_without_time = hydra::wrap_lambda( 
		[phsp] __hydra_dual__ (MSqPlus m2p, MSqMinus m2m) {
			// judge whether in phase space or not
			if (!phsp.Contains<2,3>(m2p, m2m)) return 0.0;
			return 1.;
		}
	);


	auto psi0_rnd = Psi0<DecayTime, DecayTimeError>(tau_rnd, s_rnd, b_rnd);
	auto random_background = random_background_without_time * psi0_rnd * johnson_su ;

	auto combinatorial_background_without_time = hydra::wrap_lambda( 
		[phsp] __hydra_dual__ (MSqPlus m2p, MSqMinus m2m) {
			// judge whether in phase space or not
			if (!phsp.Contains<2,3>(m2p, m2m)) return 0.0;
			return 1.;
		}
	);

	auto psi0_cmb = Psi0<DecayTime, DecayTimeError>(tau_cmb, s_cmb, b_cmb);
	auto combinatorial_background = combinatorial_background_without_time * psi0_cmb * johnson_su ;


	//---------------------------------------------------------------------------------------
	// Generate data
	//---------------------------------------------------------------------------------------
	std::cout << "***** Generation" << std::endl;

	std::cout << "Generating signal ... ... " << std::endl;
	auto start = std::chrono::high_resolution_clock::now();	
	// auto data_dz = phsp.GenerateDataWithTimeAndTimeError<MSqPlus,MSqMinus,MSqZero,DecayTime,DecayTimeError>(model_dz,args.nevents,args.seed,0.,0.5,0.06,0.09,(args.prlevel>3));
	auto data_dz = phsp.GenerateDataWithTimeAndTimeError<MSqPlus,MSqMinus,MSqZero,DecayTime,DecayTimeError>(
						model_truth_dz,efficiency,tau(),y(),b(),s(),
						johnson_su,(1-f_cmb.GetValue()-f_rnd.GetValue()) * args.nevents,
						args.seed,(args.prlevel>3));
	std::cout << "Generated " << data_dz.size() << " D0 signal candidates." << std::endl;
	// auto data_db = phsp.GenerateDataWithTimeAndTimeError<MSqPlus,MSqMinus,MSqZero,DecayTime,DecayTimeError>(model_db,args.nevents,args.seed+1,0.,0.5,0.06,0.09,(args.prlevel>3));
	auto data_db = phsp.GenerateDataWithTimeAndTimeError<MSqPlus,MSqMinus,MSqZero,DecayTime,DecayTimeError>(
						model_truth_db,efficiency,tau(),y(),b(),s(),
						johnson_su,(1-f_cmb.GetValue()-f_rnd.GetValue()) * args.nevents,
						args.seed+1,(args.prlevel>3));
	std::cout << "Generated " << data_db.size() << " D0bar signal candidates." << std::endl;

	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> elapsed = end - start;
	std::cout << "Time elapsed (ms):"<< elapsed.count() << std::endl;


	std::cout << "Generating combinatorial background ... ... " << std::endl;
	start = std::chrono::high_resolution_clock::now();

	auto dummy_efficiency = ConstantFunctor(1);

	auto data_cmb = phsp.GenerateDataWithTimeAndTimeError<MSqPlus,MSqMinus,MSqZero,DecayTime,DecayTimeError>(combinatorial_background_without_time, 
																					dummy_efficiency, tau_cmb(), 0, b_cmb(), s_cmb(), johnson_su, 
																					2* f_cmb.GetValue()*args.nevents,  // generate for _dz and _ab
																					args.seed); 

	end = std::chrono::high_resolution_clock::now();
	elapsed = end - start;
	std::cout << "Generated " << data_cmb.size() << " combinatorial background events." << std::endl;
	std::cout << "Time elapsed (ms):"<< elapsed.count() << std::endl;

	std::cout << "Generating random background ... ... " << std::endl;
	start = std::chrono::high_resolution_clock::now();

	auto data_rnd = phsp.GenerateDataWithTimeAndTimeError<MSqPlus,MSqMinus,MSqZero,DecayTime,DecayTimeError>(random_background_without_time, 
																					dummy_efficiency, tau_rnd(), 0, b_rnd(), s_rnd(), johnson_su,
																					2* f_rnd.GetValue()*args.nevents,  // generate for _dz and _ab
																					args.seed);
	
	end = std::chrono::high_resolution_clock::now();
	elapsed = end - start;
	std::cout << "Generated " << data_rnd.size() << " random background events." << std::endl;
	std::cout << "Time elapsed (ms):"<< elapsed.count() << std::endl;

	// merge and shuffle
	std::cout << "Distributing the backgrounds to data_dz and data_db ... ... " << std::endl;
	data_dz.insert( data_dz.end(), data_rnd.begin(), data_rnd.begin() + data_rnd.size()/2);
	data_dz.insert( data_dz.end(), data_cmb.begin(), data_cmb.begin() + data_cmb.size()/2);
	std::random_shuffle(data_dz.begin(), data_dz.end());

	data_db.insert( data_db.end(), data_rnd.begin() + data_rnd.size()/2 , data_rnd.end());
	data_db.insert( data_db.end(), data_cmb.begin() + data_cmb.size()/2 , data_cmb.end());
	std::random_shuffle(data_db.begin(), data_db.end());

	std::cout << "Generated " << data_dz.size() << " data_dz events in total." << std::endl;
	std::cout << "Generated " << data_db.size() << " data_db events in total." << std::endl;


	
	//---------------------------------------------------------------------------------------
	// Save to ROOT file
	//---------------------------------------------------------------------------------------

	// prepare the resolution for smearing the dalitz variables, then smear before storing
	Resolution2D resolution2D;
	config.ConfigureResolution2D(resolution2D);

	auto smear_dalitz = [&resolution2D, &phsp](double m2z, double helicity_z)->std::tuple<double,double,double,double>  {

		double m2p_s, m2m_s, m2z_s, helicity_z_s;

		do {

			auto smear_result = resolution2D.Smear(m2z, helicity_z);
			m2z_s = std::get<0>(smear_result);
			helicity_z_s = std::get<1>(smear_result);

			auto inv_helicity_z = phsp.invHelicityJK<2,3>(helicity_z_s, m2z_s);
			m2p_s = std::get<0>(inv_helicity_z); 
			m2m_s = std::get<1>(inv_helicity_z); 

		} while (!phsp.Contains<2,3>(m2p_s, m2m_s));

		return std::make_tuple(m2p_s, m2m_s, m2z_s, helicity_z_s);
	};


	{
		std::string outfilename = args.outdir + outprefix + "-data.root";
		std::cout << "Saving output to " << outfilename << std::endl;
		TFile *outfile = new TFile(outfilename.c_str(),"recreate");

		double m2p, m2m, m2z, helicity_z;
		double m2p_s, m2m_s, m2z_s, helicity_z_s; // smeared dalitz variables
		double t, sigmat;
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
			helicity_z  = phsp.HelicityJK<2,3>(m2p, m2m);

			auto smear_result 	= smear_dalitz(m2z, helicity_z);
			m2p_s 			= std::get<0>(smear_result);
			m2m_s 			= std::get<1>(smear_result);
			m2z_s 			= std::get<2>(smear_result);
			helicity_z_s 	= std::get<3>(smear_result);

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
			helicity_z  = phsp.HelicityJK<2,3>(m2p, m2m);

			auto smear_result 	= smear_dalitz(m2z, helicity_z);
			m2p_s 			= std::get<0>(smear_result);
			m2m_s 			= std::get<1>(smear_result);
			m2z_s 			= std::get<2>(smear_result);
			helicity_z_s 	= std::get<3>(smear_result);

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
	// if (args.plot) {
	// 	std::cout << "***** Plot data and model" << std::endl;

	// 	// data_dz + data_db are plotted with model_dz ignoring the CPV
	// 	auto data = data_dz;
	// 	data.insert(data.end(), data_db.begin(), data_db.end());

	// 	// plot dalitz distribution 
	// 	auto plotter = DalitzPlotterWithTimeAndTimeError<MSqPlus, MSqMinus, MSqZero, DecayTime, DecayTimeError>(phsp,"#it{K}^{0}_{S}","#it{#pi}^{+}","#it{#pi}^{#minus}",(args.prlevel>3));

	// 	std::string outfilename = args.outdir + outprefix + "-HIST.root";
	// 	// plotter.FillHistograms(data, model_dz, tau, y, b, s, outfilename, args.plotnbins); 
	// 	plotter.FillDataHistogram(data, args.plotnbins);
	// 	plotter.FillModelHistogram(model_truth_dz, efficiency, tau(), y(), b(), s(), johnson_su, args.plotnbins);
	// 	if (outfilename != "") plotter.SaveHistograms(outfilename);
	// 	plotter.SetCustomAxesTitles("#it{m}^{2}_{+} [GeV^{2}/#it{c}^{4}]","#it{m}^{2}_{#minus} [GeV^{2}/#it{c}^{4}]","#it{m}^{2}_{#it{#pi#pi}} [GeV^{2}/#it{c}^{4}]");

	// 	// plotting procedure
	// 	TApplication myapp("myapp",0,0);

	// 	// 1D Projection for dalitz distribution
	// 	TCanvas c1("c1","c1",1800,700);
	// 	TPad *pad1 = new TPad("pad1","pad1",0.01,0.25,0.33,0.99);
	// 	TPad *pad2 = new TPad("pad2","pad2",0.01,0.01,0.33,0.25);
	// 	TPad *pad3 = new TPad("pad3","pad3",0.33,0.25,0.65,0.99);
	// 	TPad *pad4 = new TPad("pad4","pad4",0.33,0.01,0.65,0.25);
	// 	TPad *pad5 = new TPad("pad5","pad5",0.65,0.25,0.97,0.99);
	// 	TPad *pad6 = new TPad("pad6","pad6",0.65,0.01,0.97,0.25);
	// 	pad1->Draw();
	// 	pad2->Draw();
	// 	pad3->Draw();
	// 	pad4->Draw();
	// 	pad5->Draw();
	// 	pad6->Draw();
	// 	pad1->SetLeftMargin(0.15);
	// 	pad2->SetLeftMargin(0.15);
	// 	pad3->SetLeftMargin(0.15);
	// 	pad4->SetLeftMargin(0.15);
	// 	pad5->SetLeftMargin(0.15);
	// 	pad6->SetLeftMargin(0.15);

	// 	pad1->cd();
	// 	plotter.Plot1DProjectionData(0, "e1");
	// 	plotter.Plot1DProjectionModel(0, "histo same");
	// 	pad2->cd();
	// 	TH1D* h1_pull = plotter.Plot1DPull(0);
	// 	plotter.PlotPullLines(h1_pull->GetXaxis()->GetXmin(), h1_pull->GetXaxis()->GetXmax());

	// 	pad3->cd();
	// 	plotter.Plot1DProjectionData(1, "e1");
	// 	plotter.Plot1DProjectionModel(1, "histo same");
	// 	pad4->cd();
	// 	h1_pull = plotter.Plot1DPull(1);
	// 	plotter.PlotPullLines(h1_pull->GetXaxis()->GetXmin(), h1_pull->GetXaxis()->GetXmax());


	// 	pad5->cd();
	// 	TH1D* h1_data = plotter.Plot1DProjectionData(2, "e1");
	// 	TH1D* h1_model = plotter.Plot1DProjectionModel(2, "histo same");
	// 	TLegend* leg = new TLegend(0.6,0.7,0.8,0.85);
	// 	leg->SetBorderSize(0);
	// 	leg->SetFillStyle(0);
	// 	leg->AddEntry(h1_data, "Data", "pe");
	// 	leg->AddEntry(h1_model, "Model Sum", "l");
	// 	leg->Draw();

	// 	pad6->cd();
	// 	h1_pull = plotter.Plot1DPull(2);
	// 	plotter.PlotPullLines(h1_pull->GetXaxis()->GetXmin(), h1_pull->GetXaxis()->GetXmax());

	// 	outfilename = args.outdir + outprefix + "-1d-projection";
	// 	Print::Canvas(c1,outfilename);

	// 	// 2D Projection for dalitz distribution
	// 	TCanvas c2("c2","c2",1500,500);
	// 	c2.Divide(3,1);

	// 	c2.cd(1);
	// 	plotter.Plot2DProjectionData(0,1); 

	// 	c2.cd(2);
	// 	plotter.Plot2DProjectionData(1,2);

	// 	c2.cd(3);
	// 	plotter.Plot2DPull(0,1);

	// 	c2.Update();
	// 	outfilename = args.outdir + outprefix + "-2d-projection";
	// 	Print::Canvas(c2,outfilename);
		
	// 	// Plot phase difference
	// 	TCanvas c3("c3","c3",600,500);
	// 	c3.SetRightMargin(.14);

	// 	auto phspWithoutTime = D0ToKsPiPi_FVECTOR_BABAR::PhaseSpace();
	// 	auto plotterWithoutTime = DalitzPlotter<MSqPlus, MSqMinus, MSqZero>(phsp,"#it{K}^{0}_{S}","#it{#pi}^{+}","#it{#pi}^{#minus}",(args.prlevel>3));
	// 	plotterWithoutTime.PlotPhaseDifference(Adir,Abar);
		
	// 	outfilename = args.outdir + outprefix + "-phase-difference";
	// 	Print::Canvas(c3,outfilename);

	// 	// time distribution
	// 	TCanvas c4("c4","c4",1200,500);
	// 	TPad *pad7 = new TPad("pad7","pad7",0.01,0.25,0.49,0.99);
	// 	TPad *pad8 = new TPad("pad8","pad8",0.01,0.01,0.49,0.25);
	// 	TPad *pad9 = new TPad("pad9","pad9",0.5,0.25,0.99,0.99);
	// 	TPad *pad10 = new TPad("pad10","pad10",0.5,0.01,0.99,0.25);	
	// 	pad7->Draw();
	// 	pad8->Draw();
	// 	pad9->Draw();
	// 	pad10->Draw();
	// 	pad7->SetLeftMargin(0.15);
	// 	pad8->SetLeftMargin(0.15);
	// 	pad9->SetLeftMargin(0.15);
	// 	pad10->SetLeftMargin(0.15);

	// 	pad7->cd();
	// 	h1_data = plotter.Plot1DProjectionData(3, "e1");
	// 	h1_model = plotter.Plot1DProjectionModel(3, "histo same");

	// 	pad8->cd();
	// 	h1_pull = plotter.Plot1DPull(3);
	// 	plotter.PlotPullLines(h1_pull->GetXaxis()->GetXmin(), h1_pull->GetXaxis()->GetXmax());

	// 	pad9->cd();
	// 	h1_data->Draw("e1");
	// 	h1_model->Draw("histo same");
	// 	pad9->SetLogy();

	// 	leg = new TLegend(0.75,0.75,0.9,0.9);
	// 	leg->SetBorderSize(0);
	// 	leg->SetFillStyle(0);
	// 	leg->AddEntry(h1_data, h1_data->GetTitle(), "pe");
	// 	leg->AddEntry(h1_model, h1_model->GetTitle(), "l");

	// 	pad10->cd();
	// 	h1_pull->Draw();
	// 	plotter.PlotPullLines(h1_pull->GetXaxis()->GetXmin(), h1_pull->GetXaxis()->GetXmax());

	// 	outfilename = args.outdir + outprefix + "-decay-time";
	// 	Print::Canvas(c4,outfilename);


	// 	// sigmat distribution
	// 	TCanvas c5("c5","c5",1200,500);
	// 	pad7 = new TPad("pad7","pad7",0.01,0.25,0.49,0.99);
	// 	pad8 = new TPad("pad8","pad8",0.01,0.01,0.49,0.25);
	// 	pad9 = new TPad("pad9","pad9",0.5,0.25,0.99,0.99);
	// 	pad10 = new TPad("pad10","pad10",0.5,0.01,0.99,0.25);	
	// 	pad7->Draw();
	// 	pad8->Draw();
	// 	pad9->Draw();
	// 	pad10->Draw();
	// 	pad7->SetLeftMargin(0.15);
	// 	pad8->SetLeftMargin(0.15);
	// 	pad9->SetLeftMargin(0.15);
	// 	pad10->SetLeftMargin(0.15);

	// 	pad7->cd();
	// 	h1_data = plotter.Plot1DProjectionData(4, "e1");
	// 	h1_model = plotter.Plot1DProjectionModel(4, "histo same");

	// 	pad8->cd();
	// 	h1_pull = plotter.Plot1DPull(4);
	// 	plotter.PlotPullLines(h1_pull->GetXaxis()->GetXmin(), h1_pull->GetXaxis()->GetXmax());

	// 	pad9->cd();
	// 	h1_data->Draw("e1");
	// 	h1_model->Draw("histo same");
	// 	pad9->SetLogy();

	// 	leg = new TLegend(0.75,0.75,0.9,0.9);
	// 	leg->SetBorderSize(0);
	// 	leg->SetFillStyle(0);
	// 	leg->AddEntry(h1_data, h1_data->GetTitle(), "pe");
	// 	leg->AddEntry(h1_model, h1_model->GetTitle(), "l");

	// 	pad10->cd();
	// 	h1_pull->Draw();
	// 	plotter.PlotPullLines(h1_pull->GetXaxis()->GetXmin(), h1_pull->GetXaxis()->GetXmax());

	// 	outfilename = args.outdir + outprefix + "-sigmat";
	// 	Print::Canvas(c5,outfilename);
		
	// 	if (args.interactive) {
	// 		std::cout << "Press Crtl+C to terminate" << std::endl;
	// 		myapp.Run();
	// 	}
		
	// }
	
	return 0;
}

