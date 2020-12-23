
template<typename MSq12, typename MSq13, typename MSq23, typename Time>
class DalitzPlotterWithTime: public PlotterBase<4>
{
private:
	// Phase-space object
	ThreeBodyPhaseSpaceWithTime _phsp;

	// Fill amplitude histogram
	template<typename Model>
	THnSparseD *fill_histogram(Model const &model, size_t n, size_t nbins=200)
	{
		auto histo = _phsp.GenerateSparseHistogramWithTime<MSq12, MSq13, MSq23, Time>(model, n, nbins);

		auto h = _phsp.RootHistogramWithTime(_labels[0].c_str(), _labels[1].c_str(), _labels[2].c_str(), nbins);
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
	THnSparseD* FillDataHistogram(T &data, size_t nbins=200, Int_t lineColor=kBlack, Int_t lineStyle=kSolid, Int_t markerStyle=20)
	{
		std::cout << "Filling histogram for data..." << std::flush;

		_h_data = _phsp.RootHistogramWithTime(_labels[0].c_str(), _labels[1].c_str(), _labels[2].c_str(), nbins);

		_h_data->SetName("h_data");
		_h_data->SetTitle("Data");
		_att_keeper_data.SetLineColor(lineColor);
		_att_keeper_data.SetLineStyle(lineStyle);
		_att_keeper_data.SetMarkerStyle(markerStyle);

		auto histo = _phsp.SparseHistogramWithTime(nbins);
		histo.Fill(data);
		fill_root_histogram(_h_data,histo);
		std::cout << " done" << std::endl ;

		return _h_data;
	}

	template<typename T>
	THnSparseD* FillModelHistogram(T & model, const size_t nbins=200, Int_t lineColor=kRed, Int_t lineStyle=kSolid)
	{
		std::cout << "Filling histogram for model..." << std::flush;

		if (!_h_data) {
			std::cout << "WARNING: data histograms has not been filled, the model histogram cannot be correctly normalized." << std::endl;
			_h_model = fill_histogram(model,5000000, nbins);
		} else {
			double ndata = get_integral(_h_data);
			double nevents = (ndata<5e5) ? 5e6 : 10.*ndata;
			_h_model = fill_histogram(model, nevents, nbins);
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
	void FillHistograms(T1 & data, T2 & model, const std::string outfilename="", const size_t nbins=200, const bool plotComponents=false)
	{
		FillDataHistogram(data, nbins);

		FillModelHistogram(model, nbins);

		if (plotComponents) {
			std::cout << "WARNING: components plotting for time-dependent model is not supported yet." << std::endl;
		}

		if (outfilename != "") SaveHistograms(outfilename);
	}
};
