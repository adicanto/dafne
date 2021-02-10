#pragma once

#include <algorithm>
#include <string>
#include <cmath>
#include <cstdlib>
#include <utility>

#include <TCanvas.h>
#include <TGraph.h>

namespace dafne {

namespace Print {

void Canvas(const TCanvas &c, std::string outfilename)
{
	// remove illegal char in ROOT macro naming convention
	std::size_t found = outfilename.find_last_of("/");
	std::replace(outfilename.begin()+found, outfilename.end(), '.', '_');
	std::replace(outfilename.begin()+found, outfilename.end(), '-', '_');
	c.Print((outfilename + ".pdf").c_str());
	c.Print((outfilename + ".C").c_str());
}

} // namespace Print

} // namespace dafne

