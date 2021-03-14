
template<typename MSq12, typename MSq13, typename MSq23>
class DalitzPlotter: public PlotterBase<3>
{
private:
	// Phase-space object
	ThreeBodyPhaseSpace _phsp;

	// Fill amplitude histogram
	template<typename Model>
	THnSparseD* fill_histogram(Model const &model, size_t n, size_t nbins=200)
	{
		auto histo = _phsp.GenerateSparseHistogram<MSq12, MSq13, MSq23>(model, n, nbins);
		auto h = _phsp.RootHistogram(_labels[0].c_str(), _labels[1].c_str(), _labels[2].c_str(), nbins);
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
	THnSparseD* FillDataHistogram(T &data, size_t nbins=200, Int_t lineColor=kBlack, Int_t lineStyle=kSolid, Int_t markerStyle=20)
	{
		std::cout << "Filling histogram for data..." << std::flush;

		_h_data = _phsp.RootHistogram(_labels[0].c_str(), _labels[1].c_str(), _labels[2].c_str(), nbins);

		_h_data->SetName("h_data");
		_h_data->SetTitle("Data");
		_att_keeper_data.SetLineColor(lineColor);
		_att_keeper_data.SetLineStyle(lineStyle);
		_att_keeper_data.SetMarkerStyle(markerStyle);

		auto histo = _phsp.SparseHistogram(nbins);
		histo.Fill(data);
		fill_root_histogram(_h_data,histo);
		std::cout << " done" << std::endl ;

		return _h_data;
	}

	template<typename T>
	THnSparseD* FillModelHistogram(T & model, size_t nbins=200, Int_t lineColor=kRed, Int_t lineStyle=kSolid)
	{
		std::cout << "Filling histogram for model..." << std::flush;

		if (!_h_data) {
			std::cout << "WARNING: data histograms has not been filled, the model histogram cannot be correctly normalized." << std::endl;
			_h_model = fill_histogram(model,5000000,nbins);
		} else {
			double ndata = get_integral(_h_data);
			double nevents = (ndata<5e5) ? 5e6 : 10.*ndata;
			_h_model = fill_histogram(model, nevents);
			_h_model->Scale( ndata/get_integral(_h_model) );
		}
		_h_model->SetName("h_model");
		_h_model->SetTitle("Total");
		_att_keeper_model.SetLineColor(lineColor);
		_att_keeper_model.SetLineStyle(lineStyle);

		std::cout << " done" << std::endl ;

		return _h_model;
	}

	template<unsigned int I, typename T>
	THnSparseD* FillComponentHistogram(T & model)
	{
		if (!_h_model) FillModelHistogram(model);

		const hydra::placeholders::placeholder<I> X;

		auto & amplitude = model.GetFunctor(X);
		const std::string amplitudeName = amplitude.Name();
		if (amplitude.IsRemoved()) {
			if (_debug) std::cout << "Component " << amplitudeName << " is removed from the model, skipping it." << std::endl;
			return 0;
		}
		const std::string amplitudeLabel = amplitude.Label();
		const char *label = (amplitudeLabel.empty()) ? amplitudeName.c_str() : amplitudeLabel.c_str();

		std::cout << "Filling histogram for component " << amplitudeName << " ( " << label << " )..." << std::flush;

		auto amp_FF = _phsp.template FitFraction<MSq12, MSq13, T>(model, amplitude);
		double nmodel = amp_FF.first * get_integral(_h_model);
		double nevents = (nmodel<1e6) ? 1e6 : nmodel;
		auto h = fill_histogram(amplitude, nevents);
		h->Scale( nmodel/get_integral(h) );

		h->SetName( Form("h_%s",amplitudeName.c_str()) );
		h->SetTitle(label);
		
		TH1D att;
		att.SetName( Form("att_%s",amplitudeName.c_str()) );
		att.SetLineColor(amplitude.Color());
		att.SetLineStyle(amplitude.Style());

		std::cout << " done (fit fraction = " << amp_FF.first << ")" << std::endl;

		_h_components.push_back(h);
		_att_keeper_components.push_back(att);
		return h;
	}

	template<typename T, int ...Is>
	void FillComponentsHistograms(T & model, std::integer_sequence<int, Is...>)
	{
		((FillComponentHistogram<Is>(model)), ...);
	}

	template<typename ...Fs>
	void FillModelAndComponentsHistograms(hydra::Sum<Fs...> & model)
	{
		FillModelHistogram(model);
		FillComponentsHistograms(model, std::make_integer_sequence<int, sizeof...(Fs)>{});
	}

	template<typename T1, typename T2>
	void FillHistograms(T1 & data, T2 & model, const std::string outfilename="", const bool plotComponents=1)
	{
		//first reduce data to the phase-space variables only
		auto reduce = hydra::wrap_lambda( [] __hydra_dual__ (MSq12 a, MSq13 b, MSq23 c)
		{
			return hydra::make_tuple(a,b,c);
		} );

		hydra::multivector<hydra::tuple<MSq12,MSq13,MSq23>, hydra::device::sys_t> reduced_data = data | reduce;

		FillDataHistogram(reduced_data);

		if (plotComponents) {
			FillModelAndComponentsHistograms(model);
		} else {
			FillModelHistogram(model);
		}

		if (outfilename != "") SaveHistograms(outfilename);
	}



	template<typename T>
	THnSparseD* FillOtherHistogram(const std::string name, const std::string title, T & functor, const double fraction, Int_t lineColor, Int_t lineStyle, Int_t fillColor, size_t nbins=200)
	{
		std::cout << "Filling histogram for other histogram: " << name << "..." << std::flush;

		_h_others[name] = nullptr;
		// TH1D att_keeper_temp();
		// _att_keeper_others[name] = att_keeper_temp;


		if (!_h_data) {
			std::cout << "WARNING: data histograms has not been filled, other histogram: " << name << " cannot be correctly normalized." << std::endl;
			_h_others[name] = fill_histogram(functor,5000000,nbins);
		} else {
			double ndata = get_integral(_h_data);
			double nevents = (ndata<5e5) ? 5e6 : 10.*ndata;
			_h_others[name] = fill_histogram(functor, nevents);
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

		
		for (auto &hmap : _h_others) {
			if (_debug) std::cout << hmap.second->GetName() << "..." << std::endl;
			auto h1d = hmap.second->Projection(xdim);
			apply_attributes(h1d, &_att_keeper_others[hmap.first]); 
			h1d->Draw("histo same");
			if (leg) leg->AddEntry(h1d,hmap.second->GetTitle(),"l");
		}

		if (legendOn) {
			if (_debug) std::cout << "legend..." << std::endl;
			leg->Draw();
		}

		if (_debug) std::cout << "done" << std::endl;
	}

	template<typename ADIR, typename ABAR>
	TH2D* PlotPhaseDifference(ADIR const &adir, ABAR const &abar, const char *goption="colz", const char* name="phase", const int nbins=300)
	{
		auto h = (TH2D*) _phsp.RootHistogram(_labels[0].c_str(), _labels[1].c_str(), _labels[2].c_str(), nbins)->Projection(0,1);
		h->SetName(name);
		h->SetZTitle("Phase difference [rad]");

		for (int i=1; i<=nbins; ++i) {
			for (int j=1; j<=nbins; ++j) {
				double x = h->GetXaxis()->GetBinCenter(i);
				double y = h->GetYaxis()->GetBinCenter(j);
				hydra::complex<double> A  = adir(x,y);
				hydra::complex<double> Ab = abar(x,y);
				double phase = hydra::arg( A * hydra::conj( Ab ) );
				h->SetBinContent(i, j, phase);
			}
		}
		h->Draw(goption);

		return h;
	}

	template<typename Amplitude>
	TH2D* PlotMagnitude(Amplitude const &amp, const char *goption="colz", const char* name="magnitude", const int nbins=300)
	{
		auto h = (TH2D*) _phsp.RootHistogram(_labels[0].c_str(), _labels[1].c_str(), _labels[2].c_str(), nbins)->Projection(0,1);
		h->SetName(name);
		h->SetZTitle("Magnitude");

		for (int i=1; i<=nbins; ++i) {
			for (int j=1; j<=nbins; ++j) {
				double x = h->GetXaxis()->GetBinCenter(i);
				double y = h->GetYaxis()->GetBinCenter(j);
				hydra::complex<double> A = amp(x,y);
				h->SetBinContent(i, j, hydra::norm(A));
			}
		}
		h->Draw(goption);

		return h;
	}

};
