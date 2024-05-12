
template<typename MSq12, typename MSq13, typename MSq23, typename Time>
class DalitzPlotterWithTime: public PlotterBase<4>
{
private:
	// Phase-space object
	ThreeBodyPhaseSpaceWithTime _phsp;

	// Fill amplitude histogram
	template<typename Model>
	THnSparseD *fill_histogram(Model const &model, double tau, double y, size_t n, 
		size_t nbins12=200, size_t nbins13=200, size_t nbins23=200, size_t nbinst=200) 
	{
		double msq12_min=_phsp.MSqMin<1,2>(); 
		double msq12_max=_phsp.MSqMax<1,2>();
		double msq13_min=_phsp.MSqMin<1,3>(); 
		double msq13_max=_phsp.MSqMax<1,3>();
		double msq23_min=_phsp.MSqMin<2,3>(); 
		double msq23_max=_phsp.MSqMax<2,3>();

		return fill_histogram(model, tau, y, n, 
								nbins12, nbins13, nbins23, nbinst, 
								msq12_min, msq12_max,
								msq13_min, msq13_max,
								msq23_min, msq23_max);
	}

	template<typename Model>
	THnSparseD *fill_histogram(Model const &model, double tau, double y, size_t n, 
		size_t nbins12, size_t nbins13, size_t nbins23, size_t nbinst, 
		double msq12_min, double msq12_max,
		double msq13_min, double msq13_max,
		double msq23_min, double msq23_max, size_t rndseed=0)
	{
		auto histo = _phsp.GenerateSparseHistogramWithTime<MSq12, MSq13, MSq23, Time>(model, tau, y, n, nbins12, nbins13, nbins23, nbinst, 
								 msq12_min, msq12_max,
								 msq13_min, msq13_max,
								 msq23_min, msq23_max, rndseed);

		auto h = _phsp.RootHistogramWithTime(_labels[0].c_str(), _labels[1].c_str(), _labels[2].c_str(), nbins12, nbins13, nbins23, nbinst, 
				 msq12_min, msq12_max,
				 msq13_min, msq13_max,
				 msq23_min, msq23_max);
		fill_root_histogram(h,histo);
		return h;
	}

public:
	DalitzPlotterWithTime() = delete;

	DalitzPlotterWithTime(ThreeBodyPhaseSpaceWithTime &phsp, std::vector<std::string> particles_labels, bool debug=false) :
		PlotterBase(particles_labels, debug), _phsp(phsp)
	{}

	DalitzPlotterWithTime(ThreeBodyPhaseSpaceWithTime &phsp, const char *labelA, const char *labelB, const char *labelC, bool debug=false) :
		DalitzPlotterWithTime(phsp,{labelA,labelB,labelC}, debug)
	{}


	template<typename T>
	THnSparseD* FillDataHistogram(T &data, 
		size_t nbins12=200, size_t nbins13=200, size_t nbins23=200, size_t nbinst=200)
	{

		double msq12_min=_phsp.MSqMin<1,2>(); 
		double msq12_max=_phsp.MSqMax<1,2>();
		double msq13_min=_phsp.MSqMin<1,3>(); 
		double msq13_max=_phsp.MSqMax<1,3>();
		double msq23_min=_phsp.MSqMin<2,3>(); 
		double msq23_max=_phsp.MSqMax<2,3>(); 

		return FillDataHistogram(data, nbins12, nbins13, nbins23, nbinst, 
								 msq12_min, msq12_max,
								 msq13_min, msq13_max,
								 msq23_min, msq23_max);
	}


	template<typename T>
	THnSparseD* FillDataHistogram(T &data, 
		size_t nbins12, size_t nbins13, size_t nbins23, size_t nbinst, 
		double msq12_min, double msq12_max,
		double msq13_min, double msq13_max,
		double msq23_min, double msq23_max, 
		Int_t lineColor=kBlack, Int_t lineStyle=kSolid, Int_t markerStyle=20)
	{
		std::cout << "Filling histogram for data..." << std::flush;

		_h_data = _phsp.RootHistogramWithTime(_labels[0].c_str(), _labels[1].c_str(), _labels[2].c_str(), nbins12, nbins13, nbins23, nbinst,
				 msq12_min, msq12_max,
				 msq13_min, msq13_max,
				 msq23_min, msq23_max);

		_h_data->SetName("h_data");
		_h_data->SetTitle("Data");
		_att_keeper_data.SetLineColor(lineColor);
		_att_keeper_data.SetLineStyle(lineStyle);
		_att_keeper_data.SetMarkerStyle(markerStyle);

		auto histo = _phsp.SparseHistogramWithTime(nbins12, nbins13, nbins23, nbinst,
													 msq12_min, msq12_max,
													 msq13_min, msq13_max,
													 msq23_min, msq23_max);
		histo.Fill(data);
		fill_root_histogram(_h_data,histo);
		std::cout << " done" << std::endl ;

		return _h_data;
	}


	template<typename T>
	THnSparseD* FillModelHistogram(T & model, double tau, double y, 
				size_t nbins12=200, size_t nbins13=200, size_t nbins23=200, size_t nbinst=200,
				Int_t lineColor=kRed, Int_t lineStyle=kSolid, size_t rndseed=0)
	{
		double msq12_min=_phsp.MSqMin<1,2>(); 
		double msq12_max=_phsp.MSqMax<1,2>();
		double msq13_min=_phsp.MSqMin<1,3>(); 
		double msq13_max=_phsp.MSqMax<1,3>();
		double msq23_min=_phsp.MSqMin<2,3>(); 
		double msq23_max=_phsp.MSqMax<2,3>();

		return FillModelHistogram(model, tau, y, 
								nbins12, nbins13, nbins23, nbinst,
								msq12_min, msq12_max,
								msq13_min, msq13_max,
								msq23_min, msq23_max); 
	}


	template<typename T>
	THnSparseD* FillModelHistogram(T & model, double tau, double y, 
				size_t nbins12, size_t nbins13, size_t nbins23, size_t nbinst, 
				double msq12_min, double msq12_max,
				double msq13_min, double msq13_max,
				double msq23_min, double msq23_max, 
				Int_t lineColor=kRed, Int_t lineStyle=kSolid, size_t rndseed=0)
	{
		std::cout << "Filling histogram for model..." << std::flush;

		if (!_h_data) {
			std::cout << "WARNING: data histograms has not been filled, the model histogram cannot be correctly normalized." << std::endl;
			_h_model = fill_histogram(model, tau, y, 5000000,
									  nbins12, nbins13, nbins23, nbinst,
									  msq12_min, msq12_max,
									  msq13_min, msq13_max,
									  msq23_min, msq23_max, rndseed);
		} else {
			double ndata = get_integral(_h_data);
			double nevents = (ndata<5e5) ? 10*5e5 : 10*ndata;
			_h_model = fill_histogram(model, tau, y, nevents, 
									  nbins12, nbins13, nbins23, nbinst,
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

	template<typename T1, typename T2>
	void FillHistograms(T1 & data, T2 & model, double tau, double y, const std::string outfilename="", 
		size_t nbins12=200, size_t nbins13=200, size_t nbins23=200, size_t nbinst=200, 
		const bool plotComponents=false)
	{

		double msq12_min=_phsp.MSqMin<1,2>(); 
		double msq12_max=_phsp.MSqMax<1,2>();
		double msq13_min=_phsp.MSqMin<1,3>(); 
		double msq13_max=_phsp.MSqMax<1,3>();
		double msq23_min=_phsp.MSqMin<2,3>(); 
		double msq23_max=_phsp.MSqMax<2,3>(); 

		FillDataHistogram(data, nbins12, nbins13, nbins23, nbinst,
							  msq12_min, msq12_max,
							  msq13_min, msq13_max,
							  msq23_min, msq23_max);

		FillModelHistogram(model, tau, y, nbins12, nbins13, nbins23, nbinst,
							  msq12_min, msq12_max,
							  msq13_min, msq13_max,
							  msq23_min, msq23_max);

		if (plotComponents) {
			std::cout << "WARNING: components plotting for time-dependent model is not supported yet." << std::endl;
		}

		if (outfilename != "") SaveHistograms(outfilename);
	}
};
