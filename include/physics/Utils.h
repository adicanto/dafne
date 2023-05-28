#pragma once

#include <hydra/detail/Config.h>
#include <cmath>

#include <hydra/detail/BackendPolicy.h>
#include <hydra/Types.h>
#include <hydra/Function.h>
#include <hydra/Pdf.h>
#include <hydra/Integrator.h>
#include <hydra/Distribution.h>
#include <hydra/detail/utility/CheckValue.h>
#include <hydra/Parameter.h>
#include <hydra/Tuple.h>
#include <tuple>
#include <limits>
#include <stdexcept>
#include <assert.h>
#include <utility>

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




} // namespace dafne

namespace hydra {


template<typename ArgType, typename Signature=double(ArgType) >
class DoubleGaussian: public BaseFunctor<DoubleGaussian<ArgType>, Signature, 5>
{
	using BaseFunctor<DoubleGaussian<ArgType>, Signature, 5>::_par;

public:

	DoubleGaussian()=delete;

	DoubleGaussian(Parameter const& mean1, Parameter const& sigma1, Parameter const& fraction2, Parameter const& mean2, Parameter const& sigma2):
		BaseFunctor<DoubleGaussian<ArgType>, Signature, 5>({mean1, sigma1, fraction2, mean2, sigma2})
		{}

	__hydra_host__ __hydra_device__
	DoubleGaussian(DoubleGaussian<ArgType> const& other ):
		BaseFunctor<DoubleGaussian<ArgType>, Signature, 5>(other)
		{}

	__hydra_host__ __hydra_device__
	DoubleGaussian<ArgType>& operator=(DoubleGaussian<ArgType> const& other )
	{
		if(this==&other) return  *this;
		BaseFunctor<DoubleGaussian<ArgType>, Signature, 5>::operator=(other);
		return  *this;
	}

	__hydra_host__ __hydra_device__
	inline double Evaluate(ArgType x)  const
	{
		double m2 = ( x - _par[0])*(x - _par[0] );
		double s2 = _par[1]*_par[1];
		double mean1 = _par[0];
		double sigma1 = _par[1];
		double fraction2 = _par[2];
		double mean2 = _par[3];
		double sigma2 = _par[4];
		double gauss1 = exp(-((x -mean1)*(x -mean1)) / (2.0*sigma1*sigma1 ));
		double gauss2 = exp(-((x -mean2)*(x -mean2)) / (2.0*sigma2*sigma2));
		return (1-fraction2)*gauss1 + fraction2*gauss2;
	}

};

// template<typename ArgType>
// class IntegrationFormula< DoubleGaussian<ArgType>, 1>
// {

// protected:

// 	inline std::pair<double, double>
// 	EvalFormula(DoubleGaussian<ArgType>const& functor, double LowerLimit, double UpperLimit )const
// 	{
// 		double fraction = cumulative(functor[0], functor[1], UpperLimit)
// 							- cumulative(functor[0], functor[1], LowerLimit);

// 			return std::make_pair( CHECK_VALUE(fraction,
// 					" par[0] = %f par[1] = %f fLowerLimit = %f fUpperLimit = %f",
// 					functor[0], functor[1], LowerLimit, UpperLimit ) ,0.0);

// 	}
// private:

// 	inline double cumulative(const double mean, const double sigma, const double x) const
// 	{
// 		static const double sqrt_pi_over_two = 1.2533141373155002512079;
// 		static const double sqrt_two         = 1.4142135623730950488017;

// 		return sigma*sqrt_pi_over_two*(1.0 + erf( (x-mean)/( sigma*sqrt_two ) ) );
// 	}
// };


template<typename ArgType>
struct RngFormula< DoubleGaussian<ArgType> >
{

	typedef ArgType value_type;
	__hydra_host__ __hydra_device__
	inline unsigned NCalls( DoubleGaussian<ArgType>const&) const
	{
		return 1;
	}

	template< typename T>
	__hydra_host__ __hydra_device__
	inline unsigned NCalls( std::initializer_list<T>) const
	{
		return 1;
	}

	template<typename Engine>
	__hydra_host__ __hydra_device__
	inline value_type Generate(Engine& rng, DoubleGaussian<ArgType>const& functor) const
	{
		double mean1  = functor[0];
		double sigma1 = functor[1];
		double fraction2 = functor[2];
		double mean2  = functor[3];
		double sigma2 = functor[4];

		double x = 0;
		if (RngBase::uniform(rng) > fraction2) {
			x = mean1 + sigma1*RngBase::normal(rng);
		} else {
			x = mean2 + sigma2*RngBase::normal(rng);
		}

		return static_cast<value_type>(x);
	}

	template<typename Engine, typename T>
	__hydra_host__ __hydra_device__
	inline value_type Generate(Engine& rng, std::initializer_list<T> pars) const
	{
		double mean1  = pars.begin()[0];
		double sigma1 = pars.begin()[1];
		double fraction2 = pars.begin()[2];
		double mean2  = pars.begin()[3];
		double sigma2 = pars.begin()[4];

		double x = 0;
		if (RngBase::uniform(rng) > fraction2) {
			x = mean1 + sigma1*RngBase::normal(rng);
		} else {
			x = mean2 + sigma2*RngBase::normal(rng);
		}


		return static_cast<value_type>(x);
	}



};




} // namespace hydra
