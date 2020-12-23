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



};


}