#pragma once

#include<vector>
#include<algorithm>
#include<iostream>
#include<fstream>
#include<cmath>
#include<TColor.h>
#include<TCanvas.h>
#include<TBox.h>
#include<TPave.h>
#include<TPaveStats.h>
#include<TPad.h>
#include<TStyle.h>
#include<TH2D.h>
#include<TGaxis.h>
#include<TTree.h>
#include<TFile.h>


#include<vector>
#include<algorithm>
#include<iostream>
#include<fstream>
#include<cmath>
#include<TColor.h>
#include<TCanvas.h>
#include<TBox.h>
#include<TPave.h>
#include<TPaveStats.h>
#include<TPad.h>
#include<TStyle.h>
#include<TH2D.h>
#include<TGaxis.h>
#include<TTree.h>
#include<TFile.h>



// a class to pass the ROOT TH2D information to hydra device
class ArbitraryBinningHistogram2D {

public:
	ArbitraryBinningHistogram2D() = delete;

	ArbitraryBinningHistogram2D(const TH2D& th2):fXTicks(th2.GetXaxis()->GetXbins()->GetSize(), 0), fYTicks(th2.GetYaxis()->GetXbins()->GetSize(), 0) {

		const double * XTicks = th2.GetXaxis()->GetXbins()->GetArray();
		for (int i = 0; i < GetNXTicks(); ++i) {
			fXTicks[i] = XTicks[i];
		}

		const double * YTicks = th2.GetYaxis()->GetXbins()->GetArray();
		for (int i = 0; i < GetNYTicks(); ++i) {
			fYTicks[i] = YTicks[i];
		}

		// initialize fZs and fZErrors
		fZs.clear();
		fZErrors.clear();

		std::vector<double> init_line(GetNYTicks(),0);
		for (int i = 0; i < GetNXTicks(); ++i) {
			fZs.push_back(init_line); 
			fZErrors.push_back(init_line); 
		}

		for (int i = 0; i < GetNXBins(); ++i)
		for (int j = 0; j < GetNYBins(); ++j) {
			fZs[i][j] = th2.GetBinContent(i+1, j+1);
			fZErrors[i][j] = th2.GetBinError(i+1, j+1);
		}

	}

	ArbitraryBinningHistogram2D(const std::vector<double> XTicks, const std::vector<double> YTicks):fXTicks(XTicks), fYTicks(YTicks){

		// initialize fZs and fZErrors
		fZs.clear();
		fZErrors.clear();

		std::vector<double> init_line(GetNYTicks(),0);
		for (int i = 0; i < GetNXTicks(); ++i) {
			fZs.push_back(init_line); 
			fZErrors.push_back(init_line); 
		}

	}

	ArbitraryBinningHistogram2D(const std::vector<double> XTicks, const std::vector<double> YTicks, std::vector< std::vector<double> > const & Zs, std::vector< std::vector<double> > const & ZErrors):fXTicks(XTicks), fYTicks(YTicks), fZs(Zs), fZErrors(ZErrors) {

		if (GetNXBins() != Zs.size() || GetNXBins() != ZErrors.size()) {
			std::cout << "NXBins seems different from the Zs/ZErrors X direction size, exit!" << std::endl;
			std::cout << "NXBins: " << GetNXBins() << " Zs X size: " << Zs.size() << " ZErrors X size: " << ZErrors.size() << std::endl;
			exit(-1);
		}

		for (unsigned ix = 0; ix < GetNXBins(); ++ix) {
			if (GetNYBins() != Zs[ix].size() || GetNYBins() != ZErrors[ix].size()) {
				std::cout << "NYBins seems different from the Zs/ZErrors Y direction size, exit!" << std::endl;

				std::cout << "NYBins: " << GetNYBins() << " Zs Y size: " << Zs.size() << " ZErrors Y size: " << ZErrors.size() << "(ix=" << ix << ")" << std::endl;
				exit(-1);
			}
		}

	}

	__hydra_dual__ inline
	ArbitraryBinningHistogram2D(ArbitraryBinningHistogram2D const & other): fXTicks(other.fXTicks), fYTicks(other.fYTicks), fZs(other.fZs), fZErrors(other.fZErrors) {}

	__hydra_dual__ inline
	ArbitraryBinningHistogram2D& operator=(ArbitraryBinningHistogram2D const & other)
	{
		if(this == &other) return *this;

		fXTicks = other.fXTicks;
		fYTicks = other.fYTicks; 
		fZs = other.fZs; 
		fZErrors = other.fZErrors;

	}

	void SetValue(const int ix, const int iy, const double value) 
	{
		fZs[ix][iy] = value;
	}

	void SetError(const int ix, const int iy, const double error) 
	{
		fZErrors[ix][iy] = error;
	}

	__hydra_dual__ inline
	bool InRange(int ix, int iy) const 
	{
		if (ix < 0 || iy < 0) {
			std::cout << "NXBins: " << GetNXBins() << ", NYBins: " << GetNYBins() << std::endl;
			std::cout << "ix: " << ix << ", iy: " << iy << std::endl;	
			return 0;
		}

		if (ix >= GetNXBins() || iy >= GetNYBins()) {
			std::cout << "NXBins: " << GetNXBins() << ", NYBins: " << GetNYBins() << std::endl;
			std::cout << "ix: " << ix << ", iy: " << iy << std::endl;
			return 0;
		}

		return 1;
	}

	__hydra_dual__ inline
	double GetValue(const double x, const double y) const
	{
		// judge if in range
  		if (!(XMin()<x&&x<XMax()&&YMin()<y&&y<YMax())) 
  			return 0;

		auto ix = unsigned(std::upper_bound(fXTicks.begin(), fXTicks.end(), x)-fXTicks.begin()) - 1;
		auto iy = unsigned(std::upper_bound(fYTicks.begin(), fYTicks.end(), y)-fYTicks.begin()) - 1;

		if (!InRange(ix, iy)) {
			std::cout << "ArbitraryBinningHistogram2D:GetValue(): ix or iy out of range. Something must be wrong!" << std::endl;
			std::cout << "x: " << x << ", y: " << y << std::endl;
			exit(-1);
		}

		return fZs[ix][iy];
	}

	__hydra_dual__ inline
	double GetError(const double x, const double y) const
	{
		// judge if in range
  		if (!(XMin()<x&&x<XMax()&&YMin()<y&&y<YMax())) 
  			return 0;

		auto ix = unsigned(std::upper_bound(fXTicks.begin(), fXTicks.end(), x)-fXTicks.begin()) - 1;
		auto iy = unsigned(std::upper_bound(fYTicks.begin(), fYTicks.end(), y)-fYTicks.begin()) - 1;

		if (!InRange(ix, iy)) {
			std::cout << "ArbitraryBinningHistogram2D:GetError(): ix or iy out of range. Something must be wrong!" << std::endl;
			std::cout << "x: " << x << ", y: " << y << std::endl;
			exit(-1);
		}

		return fZErrors[ix][iy];
	}

	__hydra_dual__ inline
	int GetBinValue(const int ix, const int iy) 
	{
		return fZs[ix][iy];
	}

	__hydra_dual__ inline
	int GetBinError(const int ix, const int iy) 
	{
		if (fZs[ix][iy] != 0 && fZErrors[ix][iy] == 0) // if error is not set, then return sqrt(z);
			return sqrt(fZs[ix][iy]);

		return fZErrors[ix][iy];
	}


	__hydra_dual__ inline
	int GetNXTicks() const {return fXTicks.size();}

	__hydra_dual__ inline
	int GetNYTicks() const {return fYTicks.size();}

	__hydra_dual__ inline
	int GetNXBins() const {return fXTicks.size()-1;}

	__hydra_dual__ inline
	int GetNYBins() const {return fYTicks.size()-1;}

	__hydra_dual__ inline
	double XMin() const {return fXTicks[0];}

	__hydra_dual__ inline
	double XMax() const {return fXTicks[fXTicks.size()-1];}

	__hydra_dual__ inline
	double YMin() const {return fYTicks[0];}

	__hydra_dual__ inline
	double YMax() const {return fYTicks[fYTicks.size()-1];}

	__hydra_dual__ inline
	double Sum() const
	{
		double sum = 0;
		for (int i = 0; i < GetNXBins(); ++i)
		for (int j = 0; j < GetNYBins(); ++j) {
			sum += fZs[i][j];
		}
	}

	void Print(std::ostream & outstream=std::cout)
	{
		outstream << "XTicks: ";
		for (int i = 0; i < fXTicks.size(); ++i) outstream << fXTicks[i] << ",";
		outstream << std::endl;

		outstream << "YTicks: ";
		for (int i = 0; i < fYTicks.size(); ++i) outstream << fYTicks[i] << ",";
		outstream << std::endl;
	
		outstream << std::endl;
		outstream << "Zs: " << std::endl;	
		for (int iy = 0; iy < GetNYBins(); ++iy) {
			for (int ix = 0; ix < GetNXBins(); ++ix) outstream << fZs[ix][iy] << ",";
			outstream << std::endl;
		}

		outstream << std::endl;
		outstream << "ZErrors: " << std::endl;	
		for (int iy = 0; iy < GetNYBins(); ++iy) {
			for (int ix = 0; ix < GetNXBins(); ++ix) outstream << fZErrors[ix][iy] << ",";
			outstream << std::endl;
		}

		outstream << std::endl;
		outstream << "# y |" << std::endl;
		outstream << "#   |" << std::endl;
		outstream << "#   |" << std::endl;
		outstream << "#   | ---------" << std::endl;
		outstream << "#             x" << std::endl;
		outstream << "# For Zs and ZErrors, the horizontal axis is x, while vertical axis is y." << std::endl;

	}

	TH2D* GetTH2D(const char* name="th2", const char* description="th2", const char* xtitle="", const char* ytitle="")
	{
		TH2D* th2 = new TH2D(name, description, GetNXBins(), fXTicks.data(), GetNYBins(), fYTicks.data());

		th2->GetXaxis()->CenterTitle();
		th2->SetXTitle(xtitle);
		th2->GetYaxis()->CenterTitle();
		th2->SetYTitle(ytitle);

		for (int i = 0; i < GetNXBins(); ++i)
		for (int j = 0; j < GetNYBins(); ++j) {
			th2->SetBinContent(i+1, j+1, fZs[i][j]);
			th2->SetBinError(i+1, j+1, fZErrors[i][j]);
		}

		return th2;
	}


private:

	std::vector<double> fXTicks;
	std::vector<double> fYTicks;

	std::vector< std::vector<double> > fZs;
	std::vector< std::vector<double> > fZErrors;

};
