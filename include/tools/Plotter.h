#include <hydra/Function.h>
#include <hydra/FunctorArithmetic.h>
#include <hydra/functions/Utils.h>
#include <hydra/Placeholders.h>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TStyle.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TH2D.h>
#include <TLine.h>

#include <physics/ThreeBodyPhaseSpace.h>
#include <tools/FunctorTools.h>

#include <map>

namespace dafne {

template<unsigned NDIM>
class PlotterBase
{
protected:
	// LaTex labels for final-state particles
	std::vector<std::string> _labels;

	THnSparseD* _h_data;
	THnSparseD* _h_model;
	std::vector<THnSparseD*> _h_components;
	std::map<std::string, THnSparseD*> _h_others;

	// Since THnSparseD does not inherit from TAttLine TAttFill TAttMarker,
	// we use the following members to keep the plotting attributes
	TH1D _att_keeper_data;
	TH1D _att_keeper_model;
	std::vector<TH1D> _att_keeper_components;
	std::map<std::string, TH1D> _att_keeper_others;
	
	bool _debug;

	double get_integral(THnSparseD *h) {
		auto proj = h->Projection(0);
		double integral = proj->Integral();
		proj->Delete();
		return integral;
	}
	
	// Convert sparse histogram into ROOT histogram
	template<typename HydraHisto>
	void fill_root_histogram(THnSparseD *h, HydraHisto &histo)
	{
		for (auto entry : histo)
		{
			size_t bin     = hydra::get<0>(entry);
			double content = hydra::get<1>(entry);
			int bins[NDIM];
			histo.GetIndexes(bin, bins);

			for (size_t i=0; i<NDIM; ++i) bins[i]++;
			h->SetBinContent(bins, content);
		}
	}
	
	template<typename HistType, typename KeeperType>
	void apply_attributes(HistType* h, const KeeperType* att_keeper)
	{
		((TAttLine*)att_keeper)->Copy(*h);
		((TAttFill*)att_keeper)->Copy(*h);
		((TAttMarker*)att_keeper)->Copy(*h);
	}
	
	TH2D* fill_pull_histogram(const TH2D *hdat, const TH2D *hfit, const char *name)
	{
		TH2D *hpull = (TH2D*) hdat->Clone(name);
		hpull->Reset();

		for (int i = 1; i <= hpull->GetXaxis()->GetNbins(); ++i) {
			for (int j = 1; j <= hpull->GetYaxis()->GetNbins(); ++j) {
				double ndat = hdat->GetBinContent(i,j);
				double edat = hdat->GetBinError(i,j);
				double nfit = hfit->GetBinContent(i,j);
				double efit = hfit->GetBinError(i,j);

				double sigma = sqrt(edat*edat+efit*efit);
				if ( sigma>0. ) hpull->SetBinContent(i, j, (ndat-nfit)/sigma);
			}
		}

		hpull->GetZaxis()->SetRangeUser(-5.,5.);
		hpull->GetZaxis()->SetTitle("Pull");

		return hpull;
	}

	TH1D* fill_pull_histogram(const TH1D *hdat, const TH1D *hfit, const char *name)
	{
		TH1D *hpull = (TH1D*) hdat->Clone(name);
		hpull->Reset();
		
		for (int i = 1; i <= hpull->GetNbinsX(); ++i) {
			double ndat = hdat->GetBinContent(i);
			double edat = hdat->GetBinError(i);
			double nfit = hfit->GetBinContent(i);
			double efit = hfit->GetBinError(i);

			double sigma = sqrt(edat*edat+efit*efit);
			if ( sigma>0. ) hpull->SetBinContent(i, (ndat-nfit)/sigma);
		}

		hpull->GetYaxis()->SetRangeUser(-5.,5.);
		hpull->GetYaxis()->SetTitle("Pull");
		hpull->GetYaxis()->SetTitleSize(0.13);
		hpull->GetYaxis()->SetLabelSize(0.13);
		hpull->GetYaxis()->SetTitleOffset(0.3);
		hpull->GetXaxis()->SetTitle("");
		hpull->GetXaxis()->SetLabelSize(0.15);
		hpull->GetXaxis()->SetTitleSize(0);
		hpull->GetXaxis()->SetLabelSize(0);

		return hpull;
	}
	
public:
	PlotterBase() = delete;

	PlotterBase(std::vector<std::string> particles_labels, bool debug=false) :  _labels(particles_labels), _h_data(nullptr), _h_model(nullptr), _h_components(), _h_others(), _att_keeper_model(), _att_keeper_data(), _att_keeper_components(), _att_keeper_others(), _debug(debug)
	{
		// set style
		gROOT->SetStyle("Modern");
		gStyle->SetOptTitle(0);
		gStyle->SetOptStat(0);
		gStyle->SetOptFit(0);
		gStyle->SetHistLineWidth(2.0);
		gStyle->SetHistMinimumZero();
	}

	void SetCustomAxesTitles(const char *xtitle, const char *ytitle, const char *ztitle, const char *ttitle=0)
	{
		_h_data->GetAxis(0)->SetTitle(xtitle);
		_h_data->GetAxis(1)->SetTitle(ytitle);
		_h_data->GetAxis(2)->SetTitle(ztitle);
		if(ttitle) _h_data->GetAxis(3)->SetTitle(ttitle);

		_h_model->GetAxis(0)->SetTitle(xtitle);
		_h_model->GetAxis(1)->SetTitle(ytitle);
		_h_model->GetAxis(2)->SetTitle(ztitle);
		if(ttitle) _h_model->GetAxis(3)->SetTitle(ttitle);
	}
	
	TH2D* Plot2DProjectionData(const int xdim, const int ydim, const char *goption="colz")
	{
		// plot data
		if (_debug) {
			std::cout << "--- Plotting 2D projection of data on " << xdim <<", "<< ydim << " axes..." << std::endl;
			std::cout << _h_data->GetName() << "..." << std::endl;
		}

		if (!_h_data) {
			std::cout << "Cannot plot 2D projection, fill data histogram first" << std::endl;
			return nullptr;
		}

		auto h2d_data = (TH2D*) gDirectory->FindObject(Form("data_proj%d%d",xdim,ydim));
		if (!h2d_data) {
			h2d_data = (TH2D*)_h_data->Projection(xdim, ydim);
			h2d_data->SetName(Form("data_proj%d%d",xdim,ydim));
		}
		h2d_data->Draw(goption);

		if (_debug) std::cout << "done" << std::endl;

		return h2d_data;
	}

	TH2D* Plot2DProjectionModel(const int xdim, const int ydim, const char *goption="colz")
	{
		// plot model
		if (_debug) {
			std::cout << "--- Plotting 2D projection of model on " << xdim <<", "<< ydim << " axes..." << std::endl;
			std::cout << _h_model->GetName() << "..." << std::endl;
		}

		if (!_h_model) {
			std::cout << "Cannot plot 2D projection, fill model histogram first" << std::endl;
			return nullptr;
		}

		auto h2d_model = (TH2D*) gDirectory->FindObject(Form("model_proj%d%d",xdim,ydim));
		if (!h2d_model) {
			h2d_model = (TH2D*)_h_model->Projection(xdim, ydim);
			h2d_model->SetName(Form("model_proj%d%d",xdim,ydim));
		}
		h2d_model->Draw(goption);

		if (_debug) std::cout << "done" << std::endl;

		return h2d_model;
	}

	TH2D* Plot2DProjectionOther(const std::string name, const int xdim, const int ydim, const char *goption="colz")
	{
		if (_debug) {
			std::cout << "--- Plotting 2D projection of model on " << xdim <<", "<< ydim << " axes..." << std::endl;
			std::cout << "other histogram: " << name << "..." << std::endl;
		}

		auto it = _h_others.find(name);
		if (it == _h_others.end()) {
			std::cout << "Cannot find other histogram: " << name << std::endl;
			return nullptr;
		}

		if (!it->second) {
			std::cout << "Only found a null pointer to other histogram: " << name << std::endl;
			return nullptr;
		}

		auto h2d_other = (TH2D*) gDirectory->FindObject(Form("other_%s_proj%d%d",name.c_str(),xdim,ydim));
		if (!h2d_other) {
			h2d_other = (TH2D*)(it->second->Projection(xdim, ydim));
			h2d_other->SetName(Form("other_%s_proj%d%d",xdim,ydim));
		}
		h2d_other->Draw(goption);

		if (_debug) std::cout << "done" << std::endl;

		return h2d_other;
	}


	TH1D* Plot1DProjectionData(const int xdim, const char *goption="pe")
	{
		// plot data
		if (_debug) {
			std::cout << "--- Plotting 1D projection of data on " << xdim << " axis..." << std::endl;
			std::cout << _h_data->GetName() << "..." << std::endl;
		}

		if (!_h_data) {
			std::cout << "Cannot plot 1D projection, fill data histogram first" << std::endl;
			return nullptr;
		}

		auto h1d_data = (TH1D*) gDirectory->FindObject(Form("data_proj%d",xdim));
		if (!h1d_data) {
			h1d_data = (TH1D*)_h_data->Projection(xdim);
			h1d_data->SetName(Form("data_proj%d",xdim));
			h1d_data->SetYTitle(Form("Candidates per %.3g GeV^{2}/#it{c}^{4}",h1d_data->GetBinWidth(1)));
			apply_attributes(h1d_data, &_att_keeper_data);
		}
		h1d_data->Draw(goption);

		if (_debug) std::cout << "done" << std::endl;

		return h1d_data;
	}

	TH1D* Plot1DProjectionModel(const int xdim, const char *goption="histo")
	{
		// plot model
		if (_debug) {
			std::cout << "--- Plotting 1D projection of model on " << xdim << " axis..." << std::endl;
			std::cout << _h_model->GetName() << "..." << std::endl;
		}

		if (!_h_model) {
			std::cout << "Cannot plot 1D projection, fill model histogram first" << std::endl;
			return nullptr;
		}

		auto h1d_model = (TH1D*) gDirectory->FindObject(Form("model_proj%d",xdim));
		if (!h1d_model) {
			h1d_model = (TH1D*)_h_model->Projection(xdim);
			h1d_model->SetName(Form("model_proj%d",xdim));
			h1d_model->SetYTitle(Form("Candidates per %.3g GeV^{2}/#it{c}^{4}",h1d_model->GetBinWidth(1)));
			apply_attributes(h1d_model, &_att_keeper_model);
		}
		h1d_model->Draw(goption);

		if (_debug) std::cout << "done" << std::endl;

		return h1d_model;
	}

	TH1D* Plot1DProjectionOther(const std::string name, const int xdim, const char *goption="histo")
	{
		if (_debug) {
			std::cout << "--- Plotting 1D projection of other histogram: " << name << " on " << xdim << " axes..." << std::endl;
			std::cout << "other histogram: " << name << "..." << std::endl;
		}

		auto it = _h_others.find(name);
		if (it == _h_others.end()) {
			std::cout << "Cannot find other histogram: " << name << std::endl;
			return nullptr;
		}

		if (!it->second) {
			std::cout << "Only found a null pointer to other histogram: " << name << std::endl;
			return nullptr;
		}

		auto h1d_other = (TH1D*) gDirectory->FindObject(Form("other_%s_proj%d",name.c_str(),xdim));
		if (!h1d_other) {
			h1d_other = (TH1D*)(it->second->Projection(xdim));
			h1d_other->SetName(Form("other_%s_proj%d",name.c_str(),xdim));
			h1d_other->SetYTitle(Form("Candidates per %.3g GeV^{2}/#it{c}^{4}",h1d_other->GetBinWidth(1)));
			apply_attributes(h1d_other, &_att_keeper_others);
		}
		h1d_other->Draw(goption);

		if (_debug) std::cout << "done" << std::endl;

		return h1d_other;
	}


	TH1D* Plot1DPull(const int xdim, const char *goption="histo")
	{
		if (_debug) {
			std::cout << "--- Plotting 1D pull of model and data on " << xdim << " axis..." << std::endl;
		}

		if (!_h_model) {
			std::cout << "Empty histogram, please fill model histogram first" << std::endl;
			return nullptr;
		}

		if (!_h_data) {
			std::cout << "Empty histogram, please fill data histogram first" << std::endl;
			return nullptr;
		}

		auto h_model = (TH1D*) gDirectory->FindObject(Form("model_proj%d",xdim));
		if (!h_model) h_model = (TH1D*) _h_model->Projection(xdim);
		auto h_data  = (TH1D*) gDirectory->FindObject(Form("data_proj%d",xdim));
		if (!h_data) h_data = (TH1D*) _h_data->Projection(xdim);

		auto h = fill_pull_histogram(h_data, h_model, Form("h_pull_%d", xdim));
		h->SetLineWidth(0);
		h->SetFillColor(1);
		h->Draw(goption);
		
		return h;
	}

	TH2D* Plot2DPull(const int xdim, const int ydim, const char *goption="colz")
	{
		if (_debug) {
			std::cout << "--- Plotting 1D pull of model and data on " << xdim << ", " << ydim << " axis..." << std::endl;
		}

		if (!_h_model) {
			std::cout << "Empty histogram, please fill model histogram first" << std::endl;
			return nullptr;
		}

		if (!_h_data) {
			std::cout << "Empty histogram, please fill data histogram first" << std::endl;
			return nullptr;
		}

		auto h_model = (TH2D*) gDirectory->FindObject(Form("model_proj%d%d",xdim,ydim));
		if (!h_model) {
			h_model = (TH2D*) _h_model->Projection(xdim, ydim);
			h_model->SetName("tmp_model");
		}
		auto h_data  = (TH2D*) gDirectory->FindObject(Form("data_proj%d%d",xdim,ydim));
		if (!h_data) {
			h_data = (TH2D*) _h_data->Projection(xdim, ydim);
			h_data->SetName("tmp_data");
		}

		auto h = fill_pull_histogram(h_data, h_model, Form("h_pull_%d_%d", xdim, ydim));
		h->Draw(goption);

		if (std::string(h_model->GetName())=="tmp_model") h_model->Delete();
		if (std::string(h_data->GetName()) =="tmp_data" ) h_data->Delete();
		
		return h;
	}

	void PlotPullLines(double mFitMin, double mFitMax)
	{
		TLine *lineUP = new TLine(mFitMin, 3, mFitMax, 3);
		lineUP->SetLineColor(38);
		lineUP->SetLineStyle(2);
		lineUP->SetLineWidth(1);
		lineUP->Draw();
		TLine *lineDOWN = new TLine(mFitMin, -3, mFitMax, -3);
		lineDOWN->SetLineColor(38);
		lineDOWN->SetLineStyle(2);
		lineDOWN->SetLineWidth(1);
		lineDOWN->Draw();
		TLine *lineZERO = new TLine(mFitMin, 0., mFitMax, 0.);
		lineZERO->SetLineColor(1);
		lineZERO->SetLineStyle(0);
		lineZERO->SetLineWidth(1);
		lineZERO->Draw();
	}

	const THnSparseD* DataHistogram() const { return _h_data; }

	const THnSparseD* ModelHistogram() const { return _h_model; }

	const std::vector<THnSparseD*> & ComponentsHistograms() const { return _h_components; }

	void SaveHistograms(const std::string filename)
	{
		TFile *f = TFile::Open(filename.c_str(), "recreate");
		if (!f) {
			std::cout << "Cannot create file " << filename << std::endl;
			return;
		}
		if (_h_data) _h_data->Write();
		if (_h_model) _h_model->Write();
		if (_h_components.size() > 0) for (auto h : _h_components) h->Write();
		if (_h_others.size() > 0) for (auto &hmap : _h_others) hmap.second->Write();
		
		f->Close();
		delete f;
	}

};

#include <tools/DalitzPlotter.h>
#include <tools/DalitzPlotterWithTime.h>
#include <tools/DalitzPlotterWithTimeAndTimeError.h>

}
