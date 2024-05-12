#pragma once

#include <random>

#include <string>
#include <fstream>
#include <cstdlib>
#include <streambuf>
#include <sstream>
#include <utility>
#include <TH3D.h>
#include <THn.h>
#include <TRandom3.h>
#include <physics/Utils.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TPad.h>

namespace dafne {

class Resolution2D
{
private:
	// xmin, xmax, ymin, ymax, rmsx, rmsy, rho
	std::vector<std::tuple<double, double, double, double, double, double, double>> items;
	TRandom3 rnd;

public:

	Resolution2D(double seed=0):rnd(seed) {};

	void Push(double const xmin, double const xmax, double const ymin, double const ymax, double const rmsx, double const rmsy, double const rho) 
	{
		items.push_back(std::make_tuple(xmin, xmax, ymin, ymax, rmsx, rmsy, rho));
	}

	auto Smear(double const x, double const y)
	{
		for (auto &item : items) {
			double xmin = std::get<0>(item);
			double xmax = std::get<1>(item);
			double ymin = std::get<2>(item);
			double ymax = std::get<3>(item);
			double rmsx = std::get<4>(item);
			double rmsy = std::get<5>(item);
			double rho = std::get<6>(item);

			if (xmin < x && x < xmax && ymin < y && y < ymax) {
				double dx0 = rnd.Gaus(0, rmsx);
				double dy0 = rnd.Gaus(0, rmsy);
				double dx = dx0;
				double dy = rmsy/rmsx*rho*dx0 + sqrt(1-rho*rho)*dy0;
				return std::make_pair(x + dx, y + dy);
			}
		}

		std::cout << "Resolution2D: failed to found a proper region for input (x, y) = (" << x << ", " << y << ") !" << std::endl;
		exit(-1);
		return std::make_pair(-9999., -9999.);
	}

	void Draw(const std::string xtitle="x", const std::string ytitle="y", const int n_sample=10000, const double max_sigmax=0.005, const double max_sigmay=0.005)
	{
		// find the x and y boundaries
		double xmin_global = std::get<0>(items[0]);
		double xmax_global = std::get<1>(items[0]);
		double ymin_global = std::get<2>(items[0]);
		double ymax_global = std::get<3>(items[0]);
		for (auto &item : items) {
			double xmin = std::get<0>(item);
			double xmax = std::get<1>(item);
			double ymin = std::get<2>(item);
			double ymax = std::get<3>(item);

			if (xmin < xmin_global) xmin_global = xmin;
			if (xmax > xmax_global) xmax_global = xmax;
			if (ymin < ymin_global) ymin_global = ymin;
			if (ymax > ymax_global) ymax_global = ymax;
		}

		gStyle->SetOptStat(0);
		gStyle->SetOptTitle(0);
		// create the global pad
		TPad * pad_global = new TPad("pad_global", "pad_global", 0.0, 0.0, 1.0, 1.0);
		pad_global->Draw();
		pad_global->cd();
		double lmargin = 0.15;
		double rmargin = 0.1;
		double tmargin = 0.1;
		double bmargin = 0.15;
		pad_global->SetLeftMargin(lmargin);
		pad_global->SetRightMargin(rmargin);
		pad_global->SetTopMargin(tmargin);
		pad_global->SetBottomMargin(bmargin);
		// create the global axes
		TH2D * h_global = new TH2D("h_global", "h_global", 100, xmin_global, xmax_global,
														   100, ymin_global, ymax_global);
		h_global->Draw();
		h_global->GetXaxis()->CenterTitle();
		h_global->GetXaxis()->SetTitle(TString::Format("#delta%s", xtitle.c_str()).Data());
		h_global->GetXaxis()->SetTitleSize(0.05);
		h_global->GetXaxis()->SetTitleOffset(1);
		h_global->GetXaxis()->SetLabelSize(0.05);
		h_global->GetYaxis()->CenterTitle();
		h_global->GetYaxis()->SetTitle(TString::Format("#delta%s", ytitle.c_str()).Data());
		h_global->GetYaxis()->SetTitleSize(0.05);
		h_global->GetYaxis()->SetTitleOffset(1);
		h_global->GetYaxis()->SetLabelSize(0.05);

		for (auto &item : items) {
			double xmin = std::get<0>(item);
			double xmax = std::get<1>(item);
			double ymin = std::get<2>(item);
			double ymax = std::get<3>(item);
			double rmsx = std::get<4>(item);
			double rmsy = std::get<5>(item);
			double rho = std::get<6>(item);


			double xmin_pad = (xmin-xmin_global) / (xmax_global-xmin_global);
			xmin_pad = lmargin+0.02 + xmin_pad*(1-rmargin-lmargin-0.02);
			double xmax_pad = (xmax-xmin_global) / (xmax_global-xmin_global);
			xmax_pad = lmargin+0.02 + xmax_pad*(1-rmargin-lmargin-0.02);
			double ymin_pad = (ymin-ymin_global) / (ymax_global-ymin_global);
			ymin_pad = bmargin+0.02 + ymin_pad*(1-tmargin-bmargin-0.02);
			double ymax_pad = (ymax-ymin_global) / (ymax_global-ymin_global);
			ymax_pad = bmargin+0.02 + ymax_pad*(1-tmargin-bmargin-0.02);

			// if one create a TPad with new TPad("pad","pad",0.0, 0.0, 0.01, 0.01);
			// , a TPad will be created on the bottom left corner in the TCanvas
			pad_global->cd();
			TString padname = TString::Format("pad__%.3f__%.3f__%.3f__%.3f",xmin,xmax,ymin,ymax);
			padname.ReplaceAll(".", "_"); 
			padname.ReplaceAll("-", "m"); 
			TPad * pad_region = new TPad(padname.Data(), padname.Data(), xmin_pad, ymin_pad, 
																		 xmax_pad, ymax_pad);
			pad_region->Draw();
			pad_region->cd();
			pad_region->SetBottomMargin(0.15);
			pad_region->SetLeftMargin(0.15);

			TString histname = TString::Format("h__%.3f__%.3f__%.3f__%.3f",xmin,xmax,ymin,ymax);
			histname.ReplaceAll(".", "_"); 
			histname.ReplaceAll("-", "m"); 
			TH2D * h_region = new TH2D(histname.Data(), histname.Data(), 
									   100, -3*max_sigmax, 3*max_sigmax, 
									   100, -3*max_sigmay, 3*max_sigmay);
			for (int i = 0; i < n_sample; ++i) {
				double dx0 = rnd.Gaus(0, rmsx);
				double dy0 = rnd.Gaus(0, rmsy);
				double dx = dx0;
				double dy = rmsy/rmsx*rho*dx0 + sqrt(1-rho*rho)*dy0;
				h_region->Fill(dx, dy);
			}

			h_region->GetXaxis()->CenterTitle();
			h_region->GetXaxis()->SetTitle(TString::Format("#delta%s", xtitle.c_str()).Data());
			h_region->GetXaxis()->SetTitleSize(0.05);
			h_region->GetXaxis()->SetTitleOffset(1);
			h_region->GetXaxis()->SetLabelSize(0.05);
			h_region->GetYaxis()->CenterTitle();
			h_region->GetYaxis()->SetTitle(TString::Format("#delta%s", ytitle.c_str()).Data());
			h_region->GetYaxis()->SetTitleSize(0.05);
			h_region->GetYaxis()->SetTitleOffset(1);
			h_region->GetYaxis()->SetLabelSize(0.05);
			h_region->Draw("COLZ");
		}

		
	}



};


}