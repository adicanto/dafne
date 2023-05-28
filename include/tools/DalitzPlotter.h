
template<typename MSq12, typename MSq13, typename MSq23>
class DalitzPlotter: public PlotterBase<3>
{
private:
	// Phase-space object
	ThreeBodyPhaseSpace _phsp;

	// Fill amplitude histogram
	template<typename Model>
	THnSparseD* fill_histogram(Model const &model, size_t n,  
		size_t nbins12=200, size_t nbins13=200, size_t nbins23=200)
	{
		double msq12_min=_phsp.MSqMin<1,2>(); 
		double msq12_max=_phsp.MSqMax<1,2>();
		double msq13_min=_phsp.MSqMin<1,3>(); 
		double msq13_max=_phsp.MSqMax<1,3>();
		double msq23_min=_phsp.MSqMin<2,3>(); 
		double msq23_max=_phsp.MSqMax<2,3>(); 

		return fill_histogram(model, n,  
							  nbins12, nbins13, nbins23,
							  msq12_min, msq12_max,
							  msq13_min, msq13_max,
							  msq23_min, msq23_max);
	}

	// Fill amplitude histogram
	template<typename Model>
	THnSparseD* fill_histogram(Model const &model, size_t n,  
		size_t nbins12, size_t nbins13, size_t nbins23,
		double msq12_min, double msq12_max,
		double msq13_min, double msq13_max,
		double msq23_min, double msq23_max, size_t rndseed=0)
	{
		auto histo = _phsp.GenerateSparseHistogram<MSq12, MSq13, MSq23>(model, n, nbins12, nbins13, nbins23, 
																		 msq12_min, msq12_max,
																		 msq13_min, msq13_max,
																		 msq23_min, msq23_max, rndseed);
		auto h = _phsp.RootHistogram(_labels[0].c_str(), _labels[1].c_str(), _labels[2].c_str(), 
									 nbins12, nbins13, nbins23, 
									 msq12_min, msq12_max,
									 msq13_min, msq13_max,
									 msq23_min, msq23_max);
		fill_root_histogram(h,histo);
		return h;
	}

public:
	DalitzPlotter() = delete;

	DalitzPlotter(ThreeBodyPhaseSpace &phsp, std::vector<std::string> particles_labels, bool debug=false) :
		PlotterBase(particles_labels, debug), _phsp(phsp)
	{}

	DalitzPlotter(ThreeBodyPhaseSpace &phsp, const char *labelA, const char *labelB, const char *labelC, bool debug=false) :
		DalitzPlotter(phsp,{labelA,labelB,labelC}, debug)
	{}
	

	template<typename T>
	THnSparseD* FillDataHistogram(T &data, size_t nbins12=200, size_t nbins13=200, size_t nbins23=200)
	{
		double msq12_min=_phsp.MSqMin<1,2>(); 
		double msq12_max=_phsp.MSqMax<1,2>();
		double msq13_min=_phsp.MSqMin<1,3>(); 
		double msq13_max=_phsp.MSqMax<1,3>();
		double msq23_min=_phsp.MSqMin<2,3>(); 
		double msq23_max=_phsp.MSqMax<2,3>();

		return FillDataHistogram(data, nbins12, nbins13, nbins23,
								 msq12_min, msq12_max,
								 msq13_min, msq13_max,
								 msq23_min, msq23_max);
	}

	template<typename T>
	THnSparseD* FillDataHistogram(T &data, 
		size_t nbins12, size_t nbins13, size_t nbins23,
		double msq12_min, double msq12_max,
		double msq13_min, double msq13_max,
		double msq23_min, double msq23_max,
		Int_t lineColor=kBlack, Int_t lineStyle=kSolid, Int_t markerStyle=20)
	{
		std::cout << "Filling histogram for data..." << std::flush;

		_h_data = _phsp.RootHistogram(_labels[0].c_str(), _labels[1].c_str(), _labels[2].c_str(), 
									 nbins12, nbins13, nbins23, 
									 msq12_min, msq12_max,
									 msq13_min, msq13_max,
									 msq23_min, msq23_max);

		_h_data->SetName("h_data");
		_h_data->SetTitle("Data");
		_att_keeper_data.SetLineColor(lineColor);
		_att_keeper_data.SetLineStyle(lineStyle);
		_att_keeper_data.SetMarkerStyle(markerStyle);

		auto histo = _phsp.SparseHistogram(nbins12, nbins13, nbins23, 
											msq12_min, msq12_max,
											msq13_min, msq13_max,
											msq23_min, msq23_max);
		histo.Fill(data);
		fill_root_histogram(_h_data,histo);
		std::cout << " done" << std::endl ;

		return _h_data;
	}

	template<typename T>
	THnSparseD* FillModelHistogram(T & model, size_t nbins12=200, size_t nbins13=200, size_t nbins23=200)
	{
		double msq12_min=_phsp.MSqMin<1,2>(); 
		double msq12_max=_phsp.MSqMax<1,2>();
		double msq13_min=_phsp.MSqMin<1,3>(); 
		double msq13_max=_phsp.MSqMax<1,3>();
		double msq23_min=_phsp.MSqMin<2,3>(); 
		double msq23_max=_phsp.MSqMax<2,3>();

		return FillModelHistogram(model, nbins12, nbins13, nbins23,
								  msq12_min, msq12_max,
								  msq13_min, msq13_max,
								  msq23_min, msq23_max);
	}

	template<typename T>
	THnSparseD* FillModelHistogram(T & model,
						size_t nbins12, size_t nbins13, size_t nbins23,
						double msq12_min, double msq12_max,
						double msq13_min, double msq13_max,
						double msq23_min, double msq23_max,
						Int_t lineColor=kRed, Int_t lineStyle=kSolid, size_t rndseed=0)
	{
		std::cout << "Filling histogram for model..." << std::flush;

		if (!_h_data) {
			std::cout << "WARNING: data histograms has not been filled, the model histogram cannot be correctly normalized." << std::endl;
			_h_model = fill_histogram(model,5000000,
									  nbins12, nbins13, nbins23, 
									  msq12_min, msq12_max,
									  msq13_min, msq13_max,
									  msq23_min, msq23_max, rndseed);
		} else {
			double ndata = get_integral(_h_data);
			double nevents = (ndata<5e5) ? 5e6 : 10.*ndata;
			_h_model = fill_histogram(model, nevents, 
									  nbins12, nbins13, nbins23, 
									  msq12_min, msq12_max,
									  msq13_min, msq13_max,
									  msq23_min, msq23_max, rndseed);
			_h_model->Scale( ndata/get_integral(_h_model) );
		}
		_h_model->SetName("h_model");
		_h_model->SetTitle("Total");
		_att_keeper_model.SetLineColor(lineColor);
		_att_keeper_model.SetLineStyle(lineStyle);

		std::cout << " done" << std::endl ;

		return _h_model;
	}

	template<unsigned int I, typename T, typename EFFICIENCY>
	THnSparseD* FillComponentHistogram(T & amp, EFFICIENCY efficiency=ConstantFunctor(1), const double signal_fraction=1.0, 
		size_t nbins12=200, size_t nbins13=200, size_t nbins23=200) 
	{
		double msq12_min=_phsp.MSqMin<1,2>(); 
		double msq12_max=_phsp.MSqMax<1,2>();
		double msq13_min=_phsp.MSqMin<1,3>(); 
		double msq13_max=_phsp.MSqMax<1,3>();
		double msq23_min=_phsp.MSqMin<2,3>(); 
		double msq23_max=_phsp.MSqMax<2,3>();

		return FillComponentHistogram<I, T, EFFICIENCY>(amp, efficiency, signal_fraction, 
										nbins12, nbins13, nbins23,
										msq12_min, msq12_max,
										msq13_min, msq13_max,
										msq23_min, msq23_max);

	}

	template<unsigned int I, typename T, typename EFFICIENCY>
	THnSparseD* FillComponentHistogram(T & amp, EFFICIENCY efficiency, const double signal_fraction, 
		size_t nbins12, size_t nbins13, size_t nbins23,
		double msq12_min, double msq12_max,
		double msq13_min, double msq13_max,
		double msq23_min, double msq23_max,
		size_t rndseed=0)
	{
		if (!_h_model) {
			std::cout << "The model sum (_h_model) is not filled, please fill it by calling ";
			std::cout << "FillModelHistogram() at first! And be sure about the efficiency is "; 
			std::cout << "consistent between model sum and components. " << std::endl;
			std::cout << "Otherwise, I can not normlize the components correctly." << std::endl;
		}

		const hydra::placeholders::placeholder<I> X;

		// remove the rate functor and get the component
		auto component_amplitude = amp.GetFunctor(X);
		const std::string amplitudeName = component_amplitude.Name();

		if (component_amplitude.IsRemoved()) {
			if (_debug) std::cout << "Component " << amplitudeName << " is removed from the model, skipping it." << std::endl;
			return 0;
		}
		const std::string amplitudeLabel = component_amplitude.Label();
		const char *label = (amplitudeLabel.empty()) ? amplitudeName.c_str() : amplitudeLabel.c_str();

		std::cout << "Filling histogram for component " << amplitudeName << " ( " << label << " )..." << std::flush;


		auto model = rate(amp) * efficiency;
		auto component = rate(component_amplitude) * efficiency;


		auto amp_FF = _phsp.FitFraction(model, component);
		double nmodel = amp_FF.first * get_integral(_h_model);
		double nevents = (nmodel<1e6) ? 1e6 : nmodel;
		auto h = fill_histogram(component, nevents,
								  nbins12, nbins13, nbins23, 
								  msq12_min, msq12_max,
								  msq13_min, msq13_max,
								  msq23_min, msq23_max,  rndseed);
		h->Scale( nmodel/get_integral(h) * signal_fraction);

		h->SetName( Form("h_%s",amplitudeName.c_str()) );
		h->SetTitle(label);
		
		TH1D att;
		att.SetName( Form("att_%s",amplitudeName.c_str()) );
		att.SetLineColor(component_amplitude.Color());
		att.SetLineStyle(component_amplitude.Style());

		std::cout << " done (fit fraction = " << amp_FF.first << ")" << std::endl;

		_h_components.push_back(h);
		_att_keeper_components.push_back(att);
		return h;
	}



	template<typename T, typename EFFICIENCY, int ...Is>
	void FillComponentHistogramsHelper(T & amp, EFFICIENCY & efficiency, const double signal_fraction, 
										size_t nbins12, size_t nbins13, size_t nbins23,
										std::integer_sequence<int, Is...>)
	{
		((FillComponentHistogram<Is>(amp, efficiency, signal_fraction, nbins12, nbins13, nbins23)), ...);
	}


	template<typename ...Fs, typename EFFICIENCY>
	void FillComponentHistograms(hydra::Sum<Fs...> & amp, EFFICIENCY & efficiency, const double signal_fraction,
								 size_t nbins12=200, size_t nbins13=200, size_t nbins23=200)
	{
		FillComponentHistogramsHelper(amp, efficiency, signal_fraction, nbins12, nbins13, nbins23, std::make_integer_sequence<int, sizeof...(Fs)>{});
	}




	template<typename T, typename EFFICIENCY, int ...Is>
	void FillComponentHistogramsHelper(T & amp, EFFICIENCY & efficiency, const double signal_fraction, 
										size_t nbins12, size_t nbins13, size_t nbins23,
										double msq12_min, double msq12_max,
										double msq13_min, double msq13_max,
										double msq23_min, double msq23_max,
										size_t rndseed, 
										std::integer_sequence<int, Is...>)
	{
		((FillComponentHistogram<Is>(amp, efficiency, signal_fraction, nbins12, nbins13, nbins23,
									 msq12_min, msq12_max,
									 msq13_min, msq13_max,
									 msq23_min, msq23_max,
									 rndseed)), ...);
	}


	template<typename ...Fs, typename EFFICIENCY>
	void FillComponentHistograms(hydra::Sum<Fs...> & amp, EFFICIENCY & efficiency, const double signal_fraction,
								 size_t nbins12, size_t nbins13, size_t nbins23,
								 double msq12_min, double msq12_max,
								 double msq13_min, double msq13_max,
								 double msq23_min, double msq23_max,
								 size_t rndseed=0)
	{
		FillComponentHistogramsHelper(amp, efficiency, signal_fraction, nbins12, nbins13, nbins23,
									 msq12_min, msq12_max,
									 msq13_min, msq13_max,
									 msq23_min, msq23_max,
									 rndseed, std::make_integer_sequence<int, sizeof...(Fs)>{});
	}


	// template<typename T, int ...Is>
	// void FillComponentsHistograms(T & amp, std::integer_sequence<int, Is...>)
	// {
	// 	((FillComponentHistogram<Is>(amp)), ...);
	// }

	// template<typename ...Fs>
	// void FillModelAndComponentsHistograms(hydra::Sum<Fs...> & amp, EFFICIENCY efficiency=ConstantFunctor(1), const double signal_fraction=1.0)
	// {
	// 	FillModelHistogram(amp);
	// 	FillComponentsHistogramsWithEfficiency(amp, efficiency, signal_fraction);
	// }

	// template<typename T1, typename T2>
	// void FillHistograms(T1 & data, T2 & model, const std::string outfilename="", const bool plotComponents=1)
	// {
	// 	//first reduce data to the phase-space variables only
	// 	auto reduce = hydra::wrap_lambda( [] __hydra_dual__ (MSq12 a, MSq13 b, MSq23 c)
	// 	{
	// 		return hydra::make_tuple(a,b,c);
	// 	} );

	// 	hydra::multivector<hydra::tuple<MSq12,MSq13,MSq23>, hydra::device::sys_t> reduced_data = data | reduce;

	// 	FillDataHistogram(reduced_data);

	// 	if (plotComponents) {
	// 		FillModelAndComponentsHistograms(model);
	// 	} else {
	// 		FillModelHistogram(model);
	// 	}

	// 	if (outfilename != "") SaveHistograms(outfilename);
	// }


	template<typename T>
	THnSparseD* FillOtherHistogram(const std::string name, const std::string title, T & functor, const double fraction, Int_t lineColor, Int_t lineStyle, Int_t fillColor, size_t nbins12=200, size_t nbins13=200, size_t nbins23=200) 
	{

		double msq12_min=_phsp.MSqMin<1,2>(); 
		double msq12_max=_phsp.MSqMax<1,2>();
		double msq13_min=_phsp.MSqMin<1,3>(); 
		double msq13_max=_phsp.MSqMax<1,3>();
		double msq23_min=_phsp.MSqMin<2,3>(); 
		double msq23_max=_phsp.MSqMax<2,3>();

		return FillOtherHistogram(name, title, functor, fraction, lineColor, lineStyle, fillColor, 
									nbins12, nbins13, nbins23,
									msq12_min, msq12_max,
									msq13_min, msq13_max,
									msq23_min, msq23_max);
	}


	template<typename T>
	THnSparseD* FillOtherHistogram(const std::string name, const std::string title, T & functor, const double fraction, Int_t lineColor, Int_t lineStyle, Int_t fillColor, size_t nbins12, size_t nbins13, size_t nbins23,
		double msq12_min, double msq12_max,
		double msq13_min, double msq13_max,
		double msq23_min, double msq23_max,
		size_t rndseed=0)
	{
		std::cout << "Filling histogram for other histogram: " << name << "..." << std::flush;

		_h_others[name] = nullptr;
		// TH1D att_keeper_temp();
		// _att_keeper_others[name] = att_keeper_temp;


		if (!_h_data) {
			std::cout << "WARNING: data histograms has not been filled, other histogram: " << name << " cannot be correctly normalized." << std::endl;
			_h_others[name] = fill_histogram(functor, 5000000, 
											  nbins12, nbins13, nbins23, 
											  msq12_min, msq12_max,
											  msq13_min, msq13_max,
											  msq23_min, msq23_max, rndseed);
		} else {
			double ndata = get_integral(_h_data);
			double nevents = (ndata<5e5) ? 5e6 : 10.*ndata;
			_h_others[name] = fill_histogram(functor, nevents, 
											  nbins12, nbins13, nbins23, 
											  msq12_min, msq12_max,
											  msq13_min, msq13_max,
											  msq23_min, msq23_max, rndseed);
			_h_others[name]->Scale( fraction*ndata/get_integral(_h_others[name]) );
		}
		_h_others[name]->SetName(Form("h_other_%s", name.c_str()));
		_h_others[name]->SetTitle(title.c_str());
		_att_keeper_others[name].SetLineColor(lineColor);
		_att_keeper_others[name].SetLineStyle(lineStyle);
		_att_keeper_others[name].SetFillColor(fillColor);

		std::cout << " done" << std::endl ;

		return _h_others[name];
	}


	void Plot1DProjections(const int xdim, bool legendOn=false)
	{
		// add legend
		TLegend* leg(0);
		if (legendOn) {
			leg = new TLegend(0.45,0.65,0.9,0.9);
			leg->SetEntrySeparation(0.25);
			leg->SetNColumns(3);
			leg->SetBorderSize(0);
			leg->SetFillStyle(0);
		}

		// plot data
		auto h1d_data = Plot1DProjectionData(xdim);
		if (leg) leg->AddEntry(h1d_data,_h_data->GetTitle(),"pe");

		for (auto &hmap : _h_others) {
			if (_debug) std::cout << hmap.second->GetName() << "..." << std::endl;
			auto h1d = hmap.second->Projection(xdim);
			apply_attributes(h1d, &_att_keeper_others[hmap.first]); 
			h1d->Draw("histo same");
			if (leg) leg->AddEntry(h1d,hmap.second->GetTitle(),"l");
		}


		// plot total model sum
		auto h1d_model = Plot1DProjectionModel(xdim, "histo same");
		if (leg) leg->AddEntry(h1d_model,_h_model->GetTitle(),"l");


		// plot components
		unsigned i(0);
		for (auto h : _h_components) {
			if (_debug) std::cout << h->GetName() << "..." << std::endl;
			auto h1d = h->Projection(xdim);
			apply_attributes(h1d, &(_att_keeper_components[i]));
			h1d->Draw("histo same");
			if (leg) leg->AddEntry(h1d,h->GetTitle(),"l");
			++i;
		}


		if (legendOn) {
			if (_debug) std::cout << "legend..." << std::endl;
			leg->Draw();
		}

		if (_debug) std::cout << "done" << std::endl;
	}

	template<typename ADIR, typename ABAR>
	TH2D* PlotPhaseDifference(ADIR const &adir, ABAR const &abar, const char *goption="colz", const char* name="phase", const int nbins12=300, const int nbins13=300, const int nbins23=300)
	{
		double msq12_min=_phsp.MSqMin<1,2>(); 
		double msq12_max=_phsp.MSqMax<1,2>();
		double msq13_min=_phsp.MSqMin<1,3>(); 
		double msq13_max=_phsp.MSqMax<1,3>();
		double msq23_min=_phsp.MSqMin<2,3>(); 
		double msq23_max=_phsp.MSqMax<2,3>();

		auto h = (TH2D*) _phsp.RootHistogram(_labels[0].c_str(), _labels[1].c_str(), _labels[2].c_str(), 
								nbins12, nbins13, nbins23,
								msq12_min, msq12_max,
								msq13_min, msq13_max,
								msq23_min, msq23_max)->Projection(0,1); // the h have (x,y) = (m2m,m2p)

		h->SetName(name);
		h->SetZTitle("Phase difference [rad]");

		for (int i=1; i<=nbins12; ++i) {
			for (int j=1; j<=nbins13; ++j) {
				double x = h->GetXaxis()->GetBinCenter(i);
				double y = h->GetYaxis()->GetBinCenter(j);
				hydra::complex<double> A  = adir(y,x);
				hydra::complex<double> Ab = abar(y,x);
				double phase = hydra::arg( A * hydra::conj( Ab ) );
				h->SetBinContent(i, j, phase);
			}
		}
		h->Draw(goption);

		return h;
	}

	template<typename Amplitude>
	TH2D* PlotMagnitude(Amplitude const &amp, const char *goption="colz", const char* name="magnitude", const int nbins12=300, const int nbins13=300, const int nbins23=300)
	{
		double msq12_min=_phsp.MSqMin<1,2>(); 
		double msq12_max=_phsp.MSqMax<1,2>();
		double msq13_min=_phsp.MSqMin<1,3>(); 
		double msq13_max=_phsp.MSqMax<1,3>();
		double msq23_min=_phsp.MSqMin<2,3>(); 
		double msq23_max=_phsp.MSqMax<2,3>();

		auto h = (TH2D*) _phsp.RootHistogram(_labels[0].c_str(), _labels[1].c_str(), _labels[2].c_str(), 
								nbins12, nbins13, nbins23,
								msq12_min, msq12_max,
								msq13_min, msq13_max,
								msq23_min, msq23_max)->Projection(0,1); // the h have (x,y) = (m2m,m2p)
		h->SetName(name);
		h->SetZTitle("Magnitude");

		for (int i=1; i<=nbins12; ++i) {
			for (int j=1; j<=nbins23; ++j) {
				double x = h->GetXaxis()->GetBinCenter(i);
				double y = h->GetYaxis()->GetBinCenter(j);
				hydra::complex<double> A = amp(y,x);
				h->SetBinContent(i, j, hydra::norm(A));
			}
		}
		h->Draw(goption);

		return h;
	}


	template<typename ADIR, typename ABAR>
	void GetBinnedPhaseInformation(ADIR const &adir, ABAR const &abar, const std::string binningFileName, 
								   std::vector<double> & _Fi, std::vector<double> & _Fmi,
								   std::vector<hydra::complex<double>> & _Xi, const int nsteps=20000)
	{

		// get Dalitz binning
	    std::string binningHistogramName = binningFileName;
	    binningHistogramName = binningHistogramName.substr(binningHistogramName.find_last_of("/")+1);
	    binningHistogramName = "h_" + binningHistogramName.substr(0, binningHistogramName.size()-5);

	    TFile *fbin = TFile::Open(binningFileName.c_str());
	    if (!fbin) {
	    	std::cout << "Failed to open the binning file : " << binningFileName << std::endl;
	    	std::cout << "Exit!" << std::endl;
	    	exit(-1);
	    }
	    TH2F *hbinning = (TH2F*) fbin->Get(binningHistogramName.c_str());
	    if (!hbinning) {
	    	std::cout << "Failed to reach the binning histogram : " << binningHistogramName << std::endl;
	    	std::cout << "Exit!" << std::endl;
	    	exit(-1);
	    }
	    // loop and compute Fi, Xi
	    int n_phaseBins = 8; 
	    std::vector<double> Fi(n_phaseBins);
	    std::vector<double> Fmi(n_phaseBins);
			std::vector<hydra::complex<double>> Xi(n_phaseBins); 

		// grid
		double msq12_min=_phsp.MSqMin<1,2>(); 
		double msq12_max=_phsp.MSqMax<1,2>();
		double msq13_min=_phsp.MSqMin<1,3>(); 
		double msq13_max=_phsp.MSqMax<1,3>();
		double msq23_min=_phsp.MSqMin<2,3>(); 
		double msq23_max=_phsp.MSqMax<2,3>();

		auto h = (TH2D*) _phsp.RootHistogram(_labels[0].c_str(), _labels[1].c_str(), _labels[2].c_str(), 
											 nsteps, nsteps, nsteps,
											 msq12_min, msq12_max,
											 msq13_min, msq13_max,
											 msq23_min, msq23_max)->Projection(0,1);

		double dm2m = h->GetXaxis()->GetBinCenter(1);
		double dm2p = h->GetXaxis()->GetBinCenter(1);

		for (int i=1; i<=nsteps; ++i) {
			for (int j=1; j<=nsteps; ++j) {
				double m2m = h->GetXaxis()->GetBinCenter(i); // the h have (x,y) = (m2m,m2p)
				double m2p = h->GetYaxis()->GetBinCenter(j);

				int d_bin = hbinning->GetBinContent(hbinning->FindBin(m2p,m2m));
		    if (d_bin<=0 || d_bin>n_phaseBins) continue;

		    if (m2p>m2m) {
		      hydra::complex<double> A  = adir(m2p,m2m);
					hydra::complex<double> Ab = abar(m2p,m2m);
					hydra::complex<double> AsAb = hydra::conj( A ) * Ab;
					double A2 = hydra::norm(A);
					double Ab2 = hydra::norm(Ab);

					Fi[d_bin - 1] += A2 * dm2m * dm2p; 
					Fmi[d_bin - 1] += Ab2 * dm2m * dm2p;
					Xi[d_bin - 1] = Xi[d_bin - 1] + AsAb * dm2m * dm2p; 
		    }
			}
		}

		for (int i = 0; i < n_phaseBins; ++i) {
			Xi[i] = Xi[i] / sqrt(Fi[i]*Fmi[i]);	
		}

		double F_sum = 0;
		for (int i = 0; i < n_phaseBins; ++i) {
			F_sum += Fi[i] + Fmi[i];	
		}		

		for (int i = 0; i < n_phaseBins; ++i) {
			Fi[i] = Fi[i]/F_sum;
			Fmi[i] = Fmi[i]/F_sum;
		}		

		_Fi = Fi;
		_Fmi = Fmi;
		_Xi = Xi;

	}

};
