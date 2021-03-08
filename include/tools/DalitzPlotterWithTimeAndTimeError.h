
template<typename MSq12, typename MSq13, typename MSq23, typename Time, typename TimeError>
class DalitzPlotterWithTimeAndTimeError: public PlotterBase<5>
{
private:
	// Phase-space object
	ThreeBodyPhaseSpaceWithTimeAndTimeError _phsp;


	// Fill amplitude histogram
	template<typename ModelTruth, typename PDFSIGMAT>
	THnSparseD *fill_histogram(ModelTruth const &modelTruth, double tau, double y, double b, double s, PDFSIGMAT const& pdf_sigma_t, size_t n, size_t nbins=200, size_t rndseed=0)
	{
		auto histo = _phsp.GenerateSparseHistogramWithTimeAndTimeError<MSq12, MSq13, MSq23, Time, TimeError>(modelTruth, tau, y, b, s, pdf_sigma_t, n, nbins, rndseed);

		auto h = _phsp.RootHistogramWithTimeAndTimeError(_labels[0].c_str(), _labels[1].c_str(), _labels[2].c_str(), nbins);
		fill_root_histogram(h,histo);
		return h;
	}


public:
	DalitzPlotterWithTimeAndTimeError() = delete;

	DalitzPlotterWithTimeAndTimeError(ThreeBodyPhaseSpaceWithTimeAndTimeError &phsp, std::vector<std::string> particles_labels, bool debug=false) :
		PlotterBase(particles_labels, debug), _phsp(phsp)
	{}

	DalitzPlotterWithTimeAndTimeError(ThreeBodyPhaseSpaceWithTimeAndTimeError &phsp, const char *labelA, const char *labelB, const char *labelC, bool debug=false) :
		DalitzPlotterWithTimeAndTimeError(phsp,{labelA,labelB,labelC}, debug)
	{}

	template<typename T>
	THnSparseD* FillDataHistogram(T &data, size_t nbins=200, Int_t lineColor=kBlack, Int_t lineStyle=kSolid, Int_t markerStyle=20)
	{
		std::cout << "Filling histogram for data..." << std::flush;

		_h_data = _phsp.RootHistogramWithTimeAndTimeError(_labels[0].c_str(), _labels[1].c_str(), _labels[2].c_str(), nbins);

		_h_data->SetName("h_data");
		_h_data->SetTitle("Data");
		_att_keeper_data.SetLineColor(lineColor);
		_att_keeper_data.SetLineStyle(lineStyle);
		_att_keeper_data.SetMarkerStyle(markerStyle);

		auto histo = _phsp.SparseHistogramWithTimeAndTimeError(nbins);
		histo.Fill(data);
		fill_root_histogram(_h_data,histo);
		std::cout << " done" << std::endl ;

		return _h_data;
	}

	template<typename ModelTruth, typename PDFSIGMAT>
	THnSparseD* FillModelHistogram(ModelTruth & model, double tau, double y, double b, double s, PDFSIGMAT const& pdf_sigma_t, const size_t nbins=200, size_t rndseed=0, Int_t lineColor=kRed, Int_t lineStyle=kSolid)
	{
		std::cout << "Filling histogram for model..." << std::flush;

		if (!_h_data) {
			std::cout << "WARNING: data histograms has not been filled, the model histogram cannot be correctly normalized." << std::endl;
			_h_model = fill_histogram(model, tau, y, b, s, pdf_sigma_t, 500000, nbins, rndseed);
		} else {
			double ndata = get_integral(_h_data);
			double nevents = (ndata<5e5) ? 1e7 : 20.*ndata;
			_h_model = fill_histogram(model, tau, y, b, s, pdf_sigma_t, nevents, nbins, rndseed);
			_h_model->Scale( ndata/get_integral(_h_model) );
		}

		_h_model->SetName("h_model");
		_h_model->SetTitle("Total");
		_att_keeper_model.SetLineColor(lineColor);
		_att_keeper_model.SetLineStyle(lineStyle);

		std::cout << " done" << std::endl ;

		return _h_model;
	}


	template<typename DATA, typename ModelTruth, typename PDFSIGMAT>
	void FillHistograms(DATA & data, ModelTruth & model,  double tau, double y, double b, double s, PDFSIGMAT const& pdf_sigma_t, const std::string outfilename="", const size_t nbins=200, const bool plotComponents=false)
	{
		FillDataHistogram(data, nbins);

		FillModelHistogram(model, tau, y, b, s, pdf_sigma_t, nbins);

		if (plotComponents) {
			std::cout << "WARNING: components plotting for time-dependent model is not supported yet." << std::endl;
		}

		if (outfilename != "") SaveHistograms(outfilename);
	}
};
