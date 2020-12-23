#pragma once

#include <hydra/detail/Config.h>
#include <cmath>

namespace dafne {

namespace Constants {

const double DegToRad = 0.017453293;
const double RadToDeg = 57.295780;

}

__hydra_dual__ inline
double kallen(double const &x, double const &y, double const &z)
{
	double w = x-y-z;
	return (w*w - 4.*y*z);
}

__hydra_dual__ inline
double p_mother_frame(double const &mSqMother, double const &mSq1, double const &mSq2)
{
	auto arg = kallen(mSqMother,mSq1,mSq2);
	return (arg>0) ? sqrt( arg / mSqMother )/2. : 0.;
}

}// namespace dafne
