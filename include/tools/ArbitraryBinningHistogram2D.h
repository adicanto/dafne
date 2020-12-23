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


class PaletteAxis {
public:
	PaletteAxis() = delete;

	PaletteAxis(const double ZMin, const double ZMax, unsigned Nbins, long FI): fZMin(ZMin), fZMax(ZMax), fNbins(Nbins), fFI(FI), fColors(Nbins)

	{
		for (unsigned i = 0; i < Nbins; ++i) fColors[i] = FI+i;
	}

	PaletteAxis(PaletteAxis const & other): fZMin(other.fZMin), fZMax(other.fZMax), fNbins(other.fNbins), fFI(other.fFI), fColors(other.fColors) {}

	std::vector<double> GetColors() {return fColors;}

	double GetColor(const double Z) {
		if (Z >= fZMax) return fColors[fNbins-1];
		if (Z <= fZMin) return fColors[0];

		double ZStep = (fZMax-fZMin) / double(fNbins);
		for (unsigned i = 0; i < fNbins; ++i) {
			double lowBound = fZMin + i*ZStep;
			double highBound = fZMin + (i+1)*ZStep;
			if (lowBound<Z && Z<=highBound) return fColors[i];
		}

		return 0;
	}

	void Draw(const double xmin, const double ymin, const double xmax, const double ymax, const double xbias)
	{
		double ystep = (ymax-ymin)/double(fNbins); 
		for (unsigned i = 0; i < fNbins; i++) {
			TPave * b = new TPave(xmin-xbias, ymin+i*ystep, xmin, ymin+(i+1)*ystep, 0);
			b->SetFillColor(fColors[i]);
			b->Draw();
		}

		TGaxis * axis = new TGaxis(xmin, ymin, xmax, ymax, fZMin, fZMax, 505, "+L");
		axis->SetName("ArbitraryBinningHistogram2D_ZAxis");
		axis->Draw();
		
	}


private:
	double fZMin;
	double fZMax;
	unsigned fNbins;
	long fFI;
	std::vector<double> fColors;

};

class ArbitraryBinningHistogram2D {

public:
	ArbitraryBinningHistogram2D() = delete;

	ArbitraryBinningHistogram2D(const unsigned NBinsX, const unsigned NBinsY): fNBinsX(NBinsX), fNBinsY(NBinsY), fXLows(NBinsX,0), fXHighs(NBinsX,0), fYLows(NBinsY,0), fYHighs(NBinsY,0), fZs(), fZErrors(){

		// initialize to zero for insurance
		// for (int i = 0; i < fXLows.size(); ++i) fXLows[i] = 0;
		// for (int i = 0; i < fXHighs.size(); ++i) fXHighs[i] = 0;
		// for (int i = 0; i < fYLows.size(); ++i) fYLows[i] = 0;
		// for (int i = 0; i < fYHighs.size(); ++i) fYHighs[i] = 0;
		// for (int i = 0; i < fZs.size(); ++i) fZs[i] = 0;
		// for (int i = 0; i < fZErrors.size(); ++i) fZErrors[i] = 0;

		// initialize fZs and fZErrors
		fZs.clear();
		fZErrors.clear();

		std::vector<double> init_line(NBinsY,0);
		for (unsigned i = 0; i < NBinsX; ++i) {
			fZs.push_back(init_line); 
			fZErrors.push_back(init_line); 
		}


	}

	ArbitraryBinningHistogram2D(const unsigned NBinsX, const unsigned NBinsY, std::vector<double> const & XLows, std::vector<double> const & XHighs, std::vector<double> const & YLows, std::vector<double> const & YHighs, std::vector< std::vector<double> > const & Zs, std::vector< std::vector<double> > const & ZErrors): fNBinsX(NBinsX), fNBinsY(NBinsY), fXLows(XLows), fXHighs(XHighs), fYLows(YLows), fYHighs(YHighs), fZs(Zs), fZErrors(ZErrors) {

		if (fNBinsX != XLows.size() || fNBinsX != XHighs.size()) {
			std::cout << "fNBinsX seems different from the vector size, exit!" << std::endl;
			std::cout << "NBinsX: " << fNBinsX << " XLows size: " << XLows.size() << " XLows size: " << XHighs.size() << std::endl;
			exit(-1);
		}

		if (fNBinsY != YLows.size() || fNBinsY != YHighs.size()) {
			std::cout << "fNBinsY seems different from the vector size, exit!" << std::endl;
			std::cout << "NBinsY: " << fNBinsY << " YLows size: " << YLows.size() << " YLows size: " << YHighs.size() << std::endl;
			exit(-1);
		}

		if (fNBinsX != Zs.size() || fNBinsX != ZErrors.size()) {
			std::cout << "fNBinsX seems different from the Zs/ZErrors X direction size, exit!" << std::endl;
			std::cout << "NBinsX: " << fNBinsX << " Zs X size: " << Zs.size() << " ZErrors X size: " << ZErrors.size() << std::endl;
			exit(-1);
		}

		for (unsigned ix = 0; ix < fNBinsX; ++ix) {
			if (fNBinsY != Zs[ix].size() || fNBinsY != ZErrors[ix].size()) {
				std::cout << "fNBinsY seems different from the Zs/ZErrors Y direction size, exit!" << std::endl;

				std::cout << "NBinsY: " << fNBinsY << " Zs Y size: " << Zs.size() << " ZErrors Y size: " << ZErrors.size() << "(ix=" << ix << ")" << std::endl;
				exit(-1);
			}
		}


	}

	__hydra_dual__ inline
	ArbitraryBinningHistogram2D(ArbitraryBinningHistogram2D const & other): fNBinsX(other.fNBinsX), fNBinsY(other.fNBinsY), fXLows(other.fXLows), fXHighs(other.fXHighs), fYLows(other.fYLows), fYHighs(other.fYHighs), fZs(other.fZs), fZErrors(other.fZErrors) {}

	__hydra_dual__ inline
	ArbitraryBinningHistogram2D& operator=(ArbitraryBinningHistogram2D const & other)
	{
		if(this == &other) return *this;

		if(!(fNBinsX==other.fNBinsX && fNBinsY==other.fNBinsY)) {
			std::cout << " \"=\" operator is not allow for ArbitraryBinningHistogram2Ds with different binning numbers" << std::endl;
			exit(-1);
		}

		fXLows = other.fXLows;
		fXHighs = other.fXHighs; 
		fYLows = other.fYLows; 
		fYHighs = other.fYHighs; 
		fZs = other.fZs; 
		fZErrors = other.fZErrors;
		
		return *this;
	}

	// currently the filling function only support the basic method, without and optimization
	// When projecting the entries in a TTree, a trick to accelerate is to firstly fill a TH2D with smart binning, and Fill ArbitraryBinningHistogram2D with the information in TH2D
	int Fill(const double x, const double y, const double weight=1) 
	{
		// judge if in range
  		if (!(XMin()<x&&x<XMax()&&YMin()<y&&y<YMax())) 
  			return -1;
		
		auto ix = unsigned(std::upper_bound(fXHighs.begin(), fXHighs.end(), x)-fXHighs.begin());
		auto iy = unsigned(std::upper_bound(fYHighs.begin(), fYHighs.end(), y)-fYHighs.begin());

		if (ix >= fNBinsX || iy >= fNBinsY) {
			std::cout << "ArbitraryBinningHistogram2D:Fill(): ix, iy out of range!" << std::endl;
			std::cout << "fNBinsX: " << fNBinsX << ", fNBinsY: " << fNBinsY << std::endl;
			std::cout << "ix: " << ix << ", iy: " << iy << std::endl;
			std::cout << "x: " << x << ", y: " << y << std::endl;
			Print();
			exit(-1);
		}


	    fZs[ix][iy] += weight;
	    return 0; // more completed return value could be developed in the future

	}

	bool IsIncremental(std::vector<double> const & v) {
		for (unsigned i = 0; i < v.size()-1; ++i) {
			if ( !(v[i] < v[i+1]) ) return false;
		}
		return true;
	}

	void SetXHighs(std::vector<double> XHighs){
		if (XHighs.size() != fNBinsX) {
			std::cout << "SetXHighs() only allow vector input with size() == fNBinsX" << std::endl;
			exit(-1);
		}

		fXHighs = XHighs;
	}

	void SetXLows(std::vector<double> XLows){
		if (XLows.size() != fNBinsX) {
			std::cout << "SetXLows() only allow vector input with size() == fNBinsX" << std::endl;
			exit(-1);
		}

		fXLows = XLows;
	}

	void SetYHighs(std::vector<double> YHighs){
		if (YHighs.size() != fNBinsY) {
			std::cout << "SetYHighs() only allow vector input with size() == fNBinsY" << std::endl;
			exit(-1);
		}

		fYHighs = YHighs;
	}

	void SetYLows(std::vector<double> YLows){
		if (YLows.size() != fNBinsY) {
			std::cout << "SetYLows() only allow vector input with size() == fNBinsY" << std::endl;
			exit(-1);
		}

		fYLows = YLows;
	}

	void SetBinValue(const int ix, const int iy, const double value) 
	{
		fZs[ix][iy] = value;
	}

	void SetBinError(const int ix, const int iy, const double error) 
	{
		fZErrors[ix][iy] = error;
	}

	__hydra_dual__ inline
	double GetValue(const double x, const double y) const
	{
		// judge if in range
  		if (!(XMin()<x&&x<XMax()&&YMin()<y&&y<YMax())) 
  			return 0;

		auto ix = unsigned(std::upper_bound(fXHighs.begin(), fXHighs.end(), x)-fXHighs.begin());
		auto iy = unsigned(std::upper_bound(fYHighs.begin(), fYHighs.end(), y)-fYHighs.begin());

		if (ix >= fNBinsX || iy >= fNBinsY) {
			std::cout << "ArbitraryBinningHistogram2D:GetValue(): ix, iy out of range!" << std::endl;
			std::cout << "fNBinsX: " << fNBinsX << ", fNBinsY: " << fNBinsY << std::endl;
			std::cout << "ix: " << ix << ", iy: " << iy << std::endl;
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

		auto ix = unsigned(std::upper_bound(fXHighs.begin(), fXHighs.end(), x)-fXHighs.begin());
		auto iy = unsigned(std::upper_bound(fYHighs.begin(), fYHighs.end(), y)-fYHighs.begin());

		if (ix >= fNBinsX || iy >= fNBinsY) {
			std::cout << "ArbitraryBinningHistogram2D:GetValue(): ix, iy out of range!" << std::endl;
			std::cout << "fNBinsX: " << fNBinsX << ", fNBinsY: " << fNBinsY << std::endl;
			std::cout << "ix: " << ix << ", iy: " << iy << std::endl;
			std::cout << "x: " << x << ", y: " << y << std::endl;
			exit(-1);
		}

		return fZErrors[ix][iy];
	}

	int GetBinValue(const int ix, const int iy) 
	{
		return fZs[ix][iy];
	}

	int GetBinError(const int ix, const int iy) 
	{
		if (fZs[ix][iy] != 0 && fZErrors[ix][iy] == 0) // if error is not set, then return sqrt(z);
			return sqrt(fZs[ix][iy]);

		return fZErrors[ix][iy];
	}

	double GetNBinsX() {return fNBinsX;}
	double GetNBinsY() {return fNBinsY;}

	__hydra_dual__ inline
	double XMin() const {return fXLows[0];} 

	__hydra_dual__ inline
	double XMax() const {return fXHighs[fNBinsX-1];} 

	__hydra_dual__ inline
	double YMin() const {return fYLows[0];} 

	__hydra_dual__ inline
	double YMax() const {return fYHighs[fNBinsY-1];} 

	double Sum() 
	{
		double sum(0.);
		for (unsigned i = 0; i < fNBinsX; ++i) {
			for (unsigned j = 0; j < fNBinsY; ++j) {
				sum += fZs[i][j];
			}
		}
		return sum;
	}

	void Print(std::ostream & outstream=std::cout)
	{
		outstream << "XLows: ";
		for (unsigned i = 0; i < fXLows.size(); ++i) outstream << fXLows[i] << ",";
		outstream << std::endl;

		outstream << "XHighs: ";
		for (unsigned i = 0; i < fXHighs.size(); ++i) outstream << fXHighs[i] << ",";
		outstream << std::endl;

		outstream << "YLows: ";
		for (unsigned i = 0; i < fYLows.size(); ++i) outstream << fYLows[i] << ",";
		outstream << std::endl;
	
		outstream << "YHighs: ";
		for (unsigned i = 0; i < fYHighs.size(); ++i) outstream << fYHighs[i] << ",";
		outstream << std::endl;

		outstream << std::endl;
		outstream << "Zs: " << std::endl;	
		for (unsigned iy = 0; iy < fNBinsY; ++iy) {
			for (unsigned ix = 0; ix < fNBinsX; ++ix) outstream << fZs[ix][iy] << ",";
			outstream << std::endl;
		}

		outstream << std::endl;
		outstream << "ZErrors: " << std::endl;	
		for (unsigned iy = 0; iy < fNBinsY; ++iy) {
			for (unsigned ix = 0; ix < fNBinsX; ++ix) outstream << fZErrors[ix][iy] << ",";
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

	TH2D* Draw(const  char*  histName="ArbitraryBinningHistogram2D_THBase", const int NBinsZ=100, long FI=-9999, const double forceZMin=-9999, const double forceZMax=-9999)
	{
		// use kBird as default
		// double stops[9] = { 0.0000, 0.1250, 0.2500, 0.3750, 0.5000, 0.6250, 0.7500, 0.8750, 1.0000 };
		// double red[9]   = { 0.2082, 0.0592, 0.0780, 0.0232, 0.1802, 0.5301, 0.8186, 0.9956, 0.9764 };
		// double green[9] = { 0.1664, 0.3599, 0.5041, 0.6419, 0.7178, 0.7492, 0.7328, 0.7862, 0.9832 };
		// double blue[9]  = { 0.5293, 0.8684, 0.8385, 0.7914, 0.6425, 0.4662, 0.3499, 0.1968, 0.0539 };
		// then FI = TColor::CreateGradientColorTable(9, stops, red, green, blue, nSteps, 1) = 1179
		// more palette setting could be found in TColor::SetPalette()@TColor.cxx 
		// if you want to use other palette, please call CreateGradientColorTable at first outside the ArbitraryBinningHistogram2D.Draw()
		double stops[9] = { 0.0000, 0.1250, 0.2500, 0.3750, 0.5000, 0.6250, 0.7500, 0.8750, 1.0000 };
		double red[9]   = { 0.2082, 0.0592, 0.0780, 0.0232, 0.1802, 0.5301, 0.8186, 0.9956, 0.9764 };
		double green[9] = { 0.1664, 0.3599, 0.5041, 0.6419, 0.7178, 0.7492, 0.7328, 0.7862, 0.9832 };
		double blue[9]  = { 0.5293, 0.8684, 0.8385, 0.7914, 0.6425, 0.4662, 0.3499, 0.1968, 0.0539 };
		if (FI == -9999) FI = TColor::CreateGradientColorTable(9, stops, red, green, blue, NBinsZ, 1);

		// search for ZMin and ZMax
		double ZMin = fZs[0][0];
		double ZMax = fZs[0][0];

		for(unsigned ix = 0; ix < fNBinsX; ++ix) {
			double lineMin = *std::min_element(fZs[ix].begin(), fZs[ix].end());
			double lineMax = *std::max_element(fZs[ix].begin(), fZs[ix].end());
			if (lineMin < ZMin) ZMin = lineMin;
			if (lineMax > ZMax) ZMax = lineMax;
		}

		double ZSpan = ZMax - ZMin;
		double graphZMin = ZMin - 1.1*ZSpan;
		double graphZMax = ZMax + 1.1*ZSpan;
		if (forceZMin!=-9999 && forceZMax!=-9999) {
			graphZMax = forceZMax;
			graphZMin = forceZMin;
		}

		PaletteAxis palette(graphZMin, graphZMax, NBinsZ, FI);

		// draw the axes through TH2D
		double XSpan = XMax() - XMin();
		double YSpan = YMax() - YMin();
		double graphXMin = XMin()-0.05*XSpan; 
		double graphXMax = XMax()+0.05*XSpan;
		double graphYMin = YMin()-0.05*YSpan;
		double graphYMax = YMax()+0.05*YSpan;
		TH2D * th2 = new TH2D(histName, histName, 100, graphXMin, graphXMax, 100, graphYMin, graphYMax);
		th2->Draw();

		for(unsigned ix = 0; ix < fNBinsX; ++ix)
		for(unsigned iy = 0; iy < fNBinsY; ++iy) {
			double x1 = fXLows[ix]; 
			double x2 = fXHighs[ix]; 
			double y1 = fYLows[iy]; 
			double y2 = fYHighs[iy]; 
			TPave * box = new TPave(x1, y1, x2, y2);
			box->SetBorderSize(0);
			box->SetFillColor(palette.GetColor(fZs[ix][iy]));
			std::cout << "color: " << palette.GetColor(fZs[ix][iy]) << std::endl;
			box->Draw();

		}

		// draw the palette axis
		double XBias = XSpan/15.;
		palette.Draw(graphXMax+XBias, graphYMin, graphXMax+XBias, graphYMax, XBias);


		Print();
		std::cout << "graphZMin: " << graphZMin << ", graphZMax: " << graphZMax << std::endl;

		return th2;

	}

	// project TTree to ArbitraryBinningHistogram2D through a smart binned TH2D
	void Project(TTree* tree, const char* branchName, const char* cutStr, const double smartNBinsX, const double smartXMin, const double smartXMax, const double smartNBinsY, const double smartYMin, const double smartYMax, const char * th2Name="smartTH2D") 
	{
		TH2D th2(th2Name, th2Name, smartNBinsX, smartXMin, smartXMax, smartNBinsY, smartYMin, smartYMax);
		tree->Project(th2Name, branchName, cutStr);

		Project(&th2);
	}

	void Project(const TH2D* th2)
	{
		unsigned nbinsx = th2->GetXaxis()->GetNbins();
		unsigned nbinsy = th2->GetYaxis()->GetNbins();
		for (unsigned ix = 0; ix < nbinsx; ++ix) {
			for (unsigned iy = 0; iy < nbinsy; ++iy) {
				Fill(th2->GetXaxis()->GetBinCenter(ix), th2->GetYaxis()->GetBinCenter(iy), th2->GetBinContent(ix, iy));
			}
		}
	}

private:
	const unsigned fNBinsX;
	const unsigned fNBinsY;

	// the following four vectors could be optimized to fXTicks and fYTicks
	std::vector<double> fXLows;
	std::vector<double> fXHighs;
	std::vector<double> fYLows;
	std::vector<double> fYHighs;

	std::vector< std::vector<double> > fZs;
	std::vector< std::vector<double> > fZErrors;

};

// ArbitraryBinningHistogram2D make_arbitrary_binning_histogram2D(){

// }
