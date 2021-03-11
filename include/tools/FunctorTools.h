#pragma once

#include <hydra/detail/Config.h>
#include <hydra/Types.h>
#include <hydra/detail/TypeTraits.h>
#include <hydra/detail/utility/Utility_Tuple.h>
#include <hydra/detail/base_functor.h>
#include <hydra/detail/Constant.h>
//#include <hydra/detail/CompositeBase.h>
#include <hydra/detail/FunctorTraits.h>
#include <hydra/detail/CompositeTraits.h>
#include <hydra/Parameter.h>
#include <hydra/Tuple.h>
#include <hydra/detail/BaseCompositeFunctor.h>
#include <hydra/detail/TupleUtility.h>
#include <hydra/detail/TupleTraits.h>
#include <hydra/detail/utility/Utility_Tuple.h>


namespace dafne {

// if hydra::Pdf<> could be used in functor arithmetic operation or be used to build Compose<>, then
// this class would be not needed anymore.

// What we want from BaseCompositeFunctor, is the ability to create a new type of functor
// from a existing functor, while inheriting all the parameters of the existing function. Besides,
// the predefined parameters hash key generation function is also useful.
// But, because BaseCompositeFunctor requires at least 2 input functor, and integrator is not 
// a functor, we need to put a dummy functor to fill the blank.
template<typename FUNCTOR, typename INTEGRATOR, typename DUMMYFUNCTOR>
class NumericalNormalizedFunctor: public hydra::BaseCompositeFunctor<
							NumericalNormalizedFunctor<FUNCTOR, INTEGRATOR, DUMMYFUNCTOR>, // typename Composite

							hydra_thrust::tuple<FUNCTOR, DUMMYFUNCTOR>, // typename FunctorList

							typename hydra::detail::merged_tuple<  // typename Signature, this part could be simplified  
								hydra_thrust::tuple<typename FUNCTOR::return_type>,  // return value

								typename hydra::detail::stripped_tuple<  // input value
								typename hydra::detail::merged_tuple<
									typename FUNCTOR::argument_type
								>::type
								>::type
							>::type
					     >

{

	typedef hydra::BaseCompositeFunctor<
				NumericalNormalizedFunctor<FUNCTOR, INTEGRATOR, DUMMYFUNCTOR>, // typename Composite

				hydra_thrust::tuple<FUNCTOR, DUMMYFUNCTOR>, // typename FunctorList

				typename hydra::detail::merged_tuple<  // typename Signature  
					hydra_thrust::tuple<typename FUNCTOR::return_type>,  // return value

					typename hydra::detail::stripped_tuple<  // input value
					typename hydra::detail::merged_tuple<
						typename FUNCTOR::argument_type
					>::type
					>::type
				>::type
		     > super_type;

public:

	NumericalNormalizedFunctor()=delete;

	NumericalNormalizedFunctor(FUNCTOR const& functor, INTEGRATOR const& integrator, DUMMYFUNCTOR const& dummy_functor): 
	super_type(functor, dummy_functor),
	fIntegrator(integrator),
	fNormCache(std::unordered_map<size_t, std::pair<double, double>>() )
	{Normalize();}

	__hydra_host__ __hydra_device__
	inline NumericalNormalizedFunctor(NumericalNormalizedFunctor<FUNCTOR,INTEGRATOR,DUMMYFUNCTOR> const& other):
    super_type(other),
    fIntegrator(other.fIntegrator),
    fNormCache(other.fNormCache)
	{Normalize();}

	__hydra_host__ __hydra_device__
	inline NumericalNormalizedFunctor<FUNCTOR,INTEGRATOR,DUMMYFUNCTOR>& operator=(NumericalNormalizedFunctor<FUNCTOR,INTEGRATOR,DUMMYFUNCTOR> const& other)
	{
		if(this==&other) return *this;
		super_type::operator=(other);
		fIntegrator=other.fIntegrator;
		fNormCache=other.fNormCache;
		return *this;
	}

	inline size_t  GetParametersKeyC() const { // a const version of GetParametersKeyC to be called in Normalize() const

		std::vector<hydra::Parameter*> _parameters;
		auto functors = this->GetFunctors();
		hydra::detail::add_parameters_in_tuple(_parameters, functors);

		std::vector<double> _temp(_parameters.size());

		for(size_t i=0; i< _parameters.size(); i++)
			_temp[i]= *(_parameters[i]);

		size_t key = hydra::detail::hash_range(_temp.begin(), _temp.end() );

		return key;
	}

	// I guess that the result of normalization could not be shared between threads, therefore 
	// maybe the functor needs to normalize inside each thread. If the Evaluate() of the functor
	// could be called from the host right before getting into device, in each loop, then, this
	// problem could be solved (the host would normalize it, then store some information in
	// fNormCache).
	inline	void Normalize( ) const
	{
		double Norm;
		double NormError;

		auto functor =  hydra::get<0>(this->GetFunctors());

		size_t key = GetParametersKeyC();

		auto search = fNormCache.find(key);
		if (search != fNormCache.end() && fNormCache.size()>0) {

			//std::cout << "found in cache "<< key << std::endl;
			std::tie(Norm, NormError) = search->second;

		}
		else {

			std::tie(Norm, NormError) =  fIntegrator(functor) ;
			fNormCache[key] = std::make_pair(Norm, NormError);

		}

		fNormThis = Norm;
		fNormErrorThis = NormError;

	}


	 template<typename ...T>
  	__hydra_host__ __hydra_device__
  	inline typename  super_type::return_type Evaluate(T... x ) const
  	{
  		Normalize();

  		auto functor =  hydra::get<0>(this->GetFunctors());

  		return functor(hydra_thrust::tie(x...))/fNormThis;
  	}


private:
	mutable INTEGRATOR fIntegrator;
	// double fNorm; // fNorm is a father class member
	mutable double fNormThis;  // define as mutable to be changed in the Evaluate
	mutable double fNormErrorThis;
	mutable std::unordered_map<size_t, std::pair<double, double>> fNormCache;

};

class ConstantFunctor: public hydra::BaseFunctor<ConstantFunctor, double(void), 0>
{
	typedef hydra::BaseFunctor<ConstantFunctor, double(void), 0> super_type;
public:
	ConstantFunctor()=delete;

	__hydra_dual__
	ConstantFunctor(const double Constant=0):fConstant(Constant){};

	__hydra_dual__ inline
	ConstantFunctor(ConstantFunctor const& other):fConstant(other.fConstant){};

	__hydra_dual__ inline
	ConstantFunctor& operator=(ConstantFunctor const& other) 
	{
		if(this==&other) return *this;
		return *this;
	}

	__hydra_dual__ inline
	double Evaluate(void) const
	{
		return fConstant;
	}

private:
	double fConstant;
	
};

template <typename FUNCTOR, typename INTEGRATOR>
inline NumericalNormalizedFunctor<FUNCTOR, INTEGRATOR, ConstantFunctor> make_numerical_normalized_functor(FUNCTOR functor, INTEGRATOR integrator) {
		ConstantFunctor dummy_functor(0);
		return NumericalNormalizedFunctor(functor, integrator, dummy_functor);
}


class SingleValue: public hydra::BaseFunctor<SingleValue, double(void), 1>
{
	typedef hydra::BaseFunctor<SingleValue, double(void), 1> super_type;
	using super_type::_par;

public:
	SingleValue() = delete;
     
	__hydra_dual__
	SingleValue(hydra::Parameter const& v) : super_type({v})
	{}
	
	__hydra_dual__
	SingleValue(SingleValue const& other) : super_type(other)
	{}
	
	__hydra_dual__ inline
	SingleValue& operator=(SingleValue const& other)
	{
		if(this==&other) return *this;
		super_type::operator=(other);
		return *this;
	}
	   
	__hydra_dual__ inline
	double Evaluate(void) const
	{
		return _par[0];
	}
};


template<typename Amplitude>
auto rate(Amplitude const& amplitude)
{
	auto norm = hydra::wrap_lambda(
			  [=] __hydra_dual__ (hydra::complex<double> amp){
						 return hydra::norm( amp );
			  }
	);
	
	return hydra::compose(norm, amplitude);
}

template<typename ...Amplitudes>
auto rate(Amplitudes const& ... amplitudes)
{
	auto tot = hydra::sum(amplitudes...);
	
	auto norm = hydra::wrap_lambda(
			  [=] __hydra_dual__ (hydra::complex<double> amp){
						 return hydra::norm( amp );
			  }
	);
	
	return hydra::compose(norm, tot);
}

template<typename TA, typename Amplitude, typename TD>
auto divideBy(Amplitude const& amplitude, const TD denominator)
{
	auto _divide = hydra::wrap_lambda(
			  [=] __hydra_dual__ (TA amp){
						 return amp / denominator;
			  }
	);
	
	return hydra::compose(_divide, amplitude);
}

template<typename TA, typename Amplitude, typename TF>
auto multiplyBy(Amplitude const& amplitude, const TF factor)
{
	auto _multiply = hydra::wrap_lambda(
			  [=] __hydra_dual__ (TA amp){
						 return amp * factor;
			  }
	);
	
	return hydra::compose(_multiply, amplitude);
}

template<typename Amplitude>
auto real_part(Amplitude const& amplitude)
{
	auto _real = hydra::wrap_lambda(
			  [=] __hydra_dual__ (hydra::complex<double> amp){
						 return amp.real();
			  }
	);
	
	return hydra::compose(_real, amplitude);
}

template<typename Amplitude>
auto imag_part(Amplitude const& amplitude)
{
	auto _imag = hydra::wrap_lambda(
			  [=] __hydra_dual__ (hydra::complex<double> amp){
						 return amp.imag();
			  }
	);
	
	return hydra::compose(_imag, amplitude);
}


template<typename Amplitude>
auto conjugate(Amplitude const& amplitude)
{
	auto _conj = hydra::wrap_lambda(
			  [=] __hydra_dual__ (hydra::complex<double> amp){
						 return hydra::conj(amp);
			  }
	);
	
	return hydra::compose(_conj, amplitude);
}

template<typename BACKEND>
struct ConstantIntegrator;

template<hydra::detail::Backend BACKEND>
class ConstantIntegrator<hydra::detail::BackendPolicy<BACKEND>>: 
public hydra::Integral<ConstantIntegrator<hydra::detail::BackendPolicy<BACKEND>>>
{
	typedef hydra::detail::BackendPolicy<BACKEND> system_t;
public:
	ConstantIntegrator()=delete;

	ConstantIntegrator(double norm):fNorm(norm){};

	// double is not allow to be template parameter, therefore we need to treat the norm as member variables
	ConstantIntegrator( ConstantIntegrator<hydra::detail::BackendPolicy<BACKEND>> const& other):fNorm(other.fNorm){}

	ConstantIntegrator<hydra::detail::BackendPolicy<BACKEND>>&
	operator=( ConstantIntegrator<hydra::detail::BackendPolicy<BACKEND>> const& other){fNorm = other.fNorm;}

	template<hydra::detail::Backend BACKEND2>
	ConstantIntegrator( ConstantIntegrator<hydra::detail::BackendPolicy<BACKEND2>> const& other):fNorm(other.fNorm){}

	template<hydra::detail::Backend BACKEND2>
	ConstantIntegrator<hydra::detail::BackendPolicy<BACKEND>>&
    operator=( ConstantIntegrator<hydra::detail::BackendPolicy<BACKEND2>> const& other){fNorm = other.fNorm;}

    template<typename FUNCTOR> 
    inline std::pair<double, double> Integrate(FUNCTOR const& fFunctor) 
    {
    	return std::make_pair(fNorm, 0.0); // constant has no error
    }
private:
	const double fNorm;
};


} // namespace dafne