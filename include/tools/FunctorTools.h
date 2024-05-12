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


#include <hydra/LogLikelihoodFCN.h>


namespace dafne {


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



template <typename FUNCTOR, typename INTEGRATOR>
inline NumericalNormalizedFunctor<FUNCTOR, INTEGRATOR, ConstantFunctor> make_numerical_normalized_functor(FUNCTOR functor, INTEGRATOR integrator) {
		ConstantFunctor dummy_functor(0);
		return NumericalNormalizedFunctor(functor, integrator, dummy_functor);
}



template<typename FUNCTOR, typename INTEGRATOR, typename DUMMYFUNCTOR>
class NumericalNormalizationFunctor: public hydra::BaseCompositeFunctor<
							NumericalNormalizationFunctor<FUNCTOR, INTEGRATOR, DUMMYFUNCTOR>, // typename Composite

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
				NumericalNormalizationFunctor<FUNCTOR, INTEGRATOR, DUMMYFUNCTOR>, // typename Composite

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

	NumericalNormalizationFunctor()=delete;

	NumericalNormalizationFunctor(FUNCTOR const& functor, INTEGRATOR const& integrator, DUMMYFUNCTOR const& dummy_functor): 
	super_type(functor, dummy_functor),
	fIntegrator(integrator),
	fNormCache(std::unordered_map<size_t, std::pair<double, double>>() )
	{Normalize();}

	__hydra_host__ __hydra_device__
	inline NumericalNormalizationFunctor(NumericalNormalizationFunctor<FUNCTOR,INTEGRATOR,DUMMYFUNCTOR> const& other):
    super_type(other),
    fIntegrator(other.fIntegrator),
    fNormCache(other.fNormCache)
	{Normalize();}

	__hydra_host__ __hydra_device__
	inline NumericalNormalizationFunctor<FUNCTOR,INTEGRATOR,DUMMYFUNCTOR>& operator=(NumericalNormalizationFunctor<FUNCTOR,INTEGRATOR,DUMMYFUNCTOR> const& other)
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

  		return fNormThis;
  	}


private:
	mutable INTEGRATOR fIntegrator;
	// double fNorm; // fNorm is a father class member
	mutable double fNormThis;  // define as mutable to be changed in the Evaluate
	mutable double fNormErrorThis;
	mutable std::unordered_map<size_t, std::pair<double, double>> fNormCache;

};


template <typename FUNCTOR, typename INTEGRATOR>
inline NumericalNormalizationFunctor<FUNCTOR, INTEGRATOR, ConstantFunctor> make_numerical_normalization_functor(FUNCTOR functor, INTEGRATOR integrator) {
		ConstantFunctor dummy_functor(0);
		return NumericalNormalizationFunctor(functor, integrator, dummy_functor);
}



class PassParameter: public hydra::BaseFunctor<PassParameter, double(void), 1>
{
	typedef hydra::BaseFunctor<PassParameter, double(void), 1> super_type;
	using super_type::_par;

public:
	PassParameter() = delete;
     
	__hydra_dual__
	PassParameter(hydra::Parameter const& v) : super_type({v})
	{}
	
	__hydra_dual__
	PassParameter(PassParameter const& other) : super_type(other)
	{}
	
	__hydra_dual__ inline
	PassParameter& operator=(PassParameter const& other)
	{
		if(this==&other) return *this;
		return *this;
	}
	   
	__hydra_dual__ inline
	double Evaluate(void) const
	{
		return _par[0];
	}
};

template<typename VARIABLES, typename Signature=double(VARIABLES)>
class PassVariable: public hydra::BaseFunctor<PassVariable<VARIABLES>, Signature, 0>
{
	typedef hydra::BaseFunctor<PassVariable<VARIABLES>, Signature, 0> super_type;
	using super_type::_par;

public:
	__hydra_dual__
	PassVariable() {}
	
	__hydra_dual__
	PassVariable(PassVariable<VARIABLES, Signature> const& other) {}
	
	__hydra_dual__ inline
	PassVariable<VARIABLES, Signature>& operator=(PassVariable<VARIABLES, Signature> const& other)
	{
		if(this==&other) return *this;
		super_type::operator=(other);
		return *this;
	}
	   
	__hydra_dual__ inline
	double Evaluate(VARIABLES v) const
	{
		return v;
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




template<typename FUNCTOR, typename INTEGRATOR, typename FBKG, typename PDFBKG>
class PdfWithBackgroundFunctor: public hydra::BaseCompositeFunctor<
							PdfWithBackgroundFunctor<FUNCTOR, INTEGRATOR, FBKG, PDFBKG>, // typename Composite
							hydra_thrust::tuple<FUNCTOR, FBKG, PDFBKG>, // typename FunctorList
							typename hydra::detail::merged_tuple<  // typename Signature
								hydra_thrust::tuple<typename FUNCTOR::return_type>,  // return value
								typename hydra::detail::stripped_tuple<  // input value
								typename hydra::detail::merged_tuple<
									typename FUNCTOR::argument_type,
									typename FBKG::argument_type
								>::type
								>::type
							>::type
					     >
{
	typedef hydra::BaseCompositeFunctor<
				PdfWithBackgroundFunctor<FUNCTOR, INTEGRATOR, FBKG, PDFBKG>, // typename Composite
				hydra_thrust::tuple<FUNCTOR, FBKG, PDFBKG>, // typename FunctorList
				typename hydra::detail::merged_tuple<  // typename Signature  
					hydra_thrust::tuple<typename FUNCTOR::return_type>,  // return value
					typename hydra::detail::stripped_tuple<  // input value
					typename hydra::detail::merged_tuple<
						typename FUNCTOR::argument_type,
						typename FBKG::argument_type
					>::type
					>::type
				>::type
		     > super_type;


public:

	PdfWithBackgroundFunctor()=delete;

	PdfWithBackgroundFunctor(FUNCTOR const& functor, INTEGRATOR const& integrator, FBKG const& f_bkg, PDFBKG const& pdf_bkg): 
	super_type(functor, f_bkg, pdf_bkg),
	fIntegrator(integrator),
	fFunctorNormCache(std::unordered_map<size_t, double>() )
	{Normalize();}

	__hydra_host__ __hydra_device__
	inline PdfWithBackgroundFunctor(PdfWithBackgroundFunctor<FUNCTOR, INTEGRATOR, FBKG, PDFBKG> const& other):
    super_type(other),
    fIntegrator(other.fIntegrator),
    fFunctorNormCache(other.fFunctorNormCache)
	{Normalize();}

	__hydra_host__ __hydra_device__
	inline PdfWithBackgroundFunctor<FUNCTOR, INTEGRATOR, FBKG, PDFBKG>& operator=(PdfWithBackgroundFunctor<FUNCTOR, INTEGRATOR, FBKG, PDFBKG> const& other)
	{
		if(this==&other) return *this;
		super_type::operator=(other);
		fIntegrator=other.fIntegrator;
	    fFunctorNormCache=other.fFunctorNormCache;
		return *this;
	}

	inline size_t  GetDalitzParametersKeyC() const 
	{ // a const version of GetDalitzParametersKeyC to be called in Normalize() const

		std::vector<hydra::Parameter*> _parameters;
		auto functor = hydra::get<0>(this->GetFunctors());
		functor.AddUserParameters(_parameters);

		std::vector<double> _temp(_parameters.size());

		for(size_t i=0; i< _parameters.size(); i++)
			_temp[i]= *(_parameters[i]);

		size_t key = hydra::detail::hash_range(_temp.begin(), _temp.end() );

		return key;
	}

	inline	void Normalize( ) const
	{

		size_t key = GetDalitzParametersKeyC();

		double FunctorNorm;

		auto search1 = fFunctorNormCache.find(key);
		if (search1 != fFunctorNormCache.end() && fFunctorNormCache.size()>0) {
			//std::cout << "found in cache "<< key << std::endl;
			FunctorNorm = search1->second;

		} else {
			auto functor = hydra::get<0>(this->GetFunctors());

			double FunctorNorm = fIntegrator(functor).first;

			fFunctorNormCache[key] = FunctorNorm;
		}


		fFunctorNormThis = FunctorNorm;

	}


	 template<typename ...T>
  	__hydra_host__ __hydra_device__
  	inline typename  super_type::return_type Evaluate(T... x) const
  	{
  		Normalize();

  		auto functor = hydra::get<0>(this->GetFunctors());
		auto f_bkg = hydra::get<1>(this->GetFunctors());
		auto pdf_bkg =  hydra::get<2>(this->GetFunctors());


		double _pdf_sig = functor(hydra_thrust::tie(x...));
		_pdf_sig = _pdf_sig/fFunctorNormThis;
		double _f_bkg = f_bkg(hydra_thrust::tie(x...));
		double _pdf_bkg = pdf_bkg(hydra_thrust::tie(x...));



		double result = (1-_f_bkg) * _pdf_sig;
		result += _f_bkg * _pdf_bkg;

  		return result;
  	}


private:
	mutable INTEGRATOR fIntegrator;
	mutable double fFunctorNormThis;  // define as mutable to be changed in the Evaluate
	mutable std::unordered_map<size_t, double> fFunctorNormCache;
};


template<typename FUNCTOR, typename INTEGRATOR, typename FBKG, typename PDFBKG>
auto make_pdf_with_background_functor(FUNCTOR const& functor, INTEGRATOR const& integrator, FBKG const& f_bkg, PDFBKG const& pdf_bkg)
{
	return PdfWithBackgroundFunctor<FUNCTOR, INTEGRATOR, FBKG, PDFBKG>(functor, integrator, f_bkg, pdf_bkg);
}







} // namespace dafne





namespace hydra {

// This class is similar to hydra::Pdf.
// The difference is, the functor to build this class has Normalized() member function.

// This PdfFromNormalizedFunctor<Functor> + LogLikelihoodFCN< PdfFromNormalizedFunctor<Functor> ...>
// scheme, could be replaced with Pdf<Functor, Integrator> + LogLikelihoodFCN< Pdf<Functor, Integrator> ...>, 
// where the Functor should implement the Evaluate(x, cache) member function. By doing so, the cache in 
// the Functor, like MixingPdfWithBackgroundFunctor, can be in principle moved to the Pdf<Functor, Integrator>.

template<typename FUNCTOR>
class PdfFromNormalizedFunctor
{

	HYDRA_STATIC_ASSERT( (detail::is_hydra_functor<FUNCTOR>::value ||
			              detail::is_hydra_lambda<FUNCTOR>::value  ||
			              detail::is_hydra_composite_functor<FUNCTOR>::value) ,
			"Instantiation of hydra::PdfFromNormalizedFunctor<FUNCTOR> class, as the naming suggests,\n "
			" requires a compliant functor or lambda as first template parameter\n"
			" and a valid hydra integrator algorithm as second  template parameter\n.")

public:

	typedef void hydra_pdf_type;

	typedef FUNCTOR functor_type;


	PdfFromNormalizedFunctor(FUNCTOR const& functor):
	fFunctor(functor),
	fNormCache(std::unordered_map<size_t, std::pair<GReal_t, GReal_t>>() )
	{Normalize();}



	PdfFromNormalizedFunctor(PdfFromNormalizedFunctor<FUNCTOR> const& other):
		fFunctor(other.GetFunctor()),
		fNorm(other.GetNorm() ),
		fNormError(other.GetNormError() ),
		fNormCache(other.GetNormCache())
	{Normalize();}

	~PdfFromNormalizedFunctor(){};


	inline PdfFromNormalizedFunctor<FUNCTOR>&
	operator=(PdfFromNormalizedFunctor<FUNCTOR> const & other )
	{
		if(this == &other) return *this;

		this->fNorm  = other.GetNorm() ;
		this->fNormError  = other.GetNormError() ;
		this->fFunctor    = other.GetFunctor();
		this->fNormCache  = other.GetNormCache();

		return *this;
	}


	inline	void AddUserParameters(std::vector<hydra::Parameter*>& user_parameters )
	{
		fFunctor.AddUserParameters(user_parameters );
	}



	inline	void PrintRegisteredParameters()
	{
		HYDRA_CALLER ;
		HYDRA_MSG << "Registered parameters begin:" << HYDRA_ENDL;
		fFunctor.PrintRegisteredParameters();
		HYDRA_MSG <<"Registered parameters end."<< HYDRA_ENDL;
		HYDRA_MSG << HYDRA_ENDL;
	}


	inline	void SetParameters(const std::vector<double>& parameters){

		this->fFunctor.SetParameters(parameters);

		this->Normalize();

		return;
	}


	inline	const FUNCTOR& GetFunctor() const {return fFunctor;}

	inline	FUNCTOR& GetFunctor() {return fFunctor;}


	inline	void Normalize( )
	{
		size_t key = fFunctor.GetParametersKey();

		auto search = fNormCache.find(key);
		if (search != fNormCache.end() && fNormCache.size()>0) {

			//std::cout << "found in cache "<< key << std::endl;
			std::tie(fNorm, fNormError) = search->second;

		}
		else {

			// Normalize the functor
			fFunctor.Normalize();

			fNorm = 1;
			fNormError = 0;
			fNormCache[key] = std::make_pair(fNorm, fNormError);
		}
		fFunctor.SetNorm(1.0/fNorm);
	}


 	template<typename T1>
 	inline  GReal_t operator()(T1&& t ) const
  	{

  		return fFunctor.GetNorm()*fFunctor(t);

  	}



  	template<typename T1, typename T2>
  	inline  GReal_t operator()( T1&& t, T2&& cache) const
  	{

  		return fFunctor.GetNorm()*fFunctor(t, cache);
  	}


   template<typename T>
   inline  GReal_t operator()( T* x, T*) const
  	{

  	  		return fFunctor.GetNorm()*fFunctor(x);
  	}


	inline GReal_t GetNorm()  const   {
		return fNorm;
	}


	inline GReal_t GetNormError()  const  {

		return fNormError;
	}


	const std::unordered_map<size_t,std::pair<GReal_t,GReal_t> >& GetNormCache()const 	{
		return fNormCache;
	}

private:

  	mutable FUNCTOR fFunctor;
	GReal_t fNorm;
	GReal_t fNormError;
	std::unordered_map<size_t, std::pair<GReal_t, GReal_t>> fNormCache;

};


/**
 * \ingroup fit
 * \brief Build a hydra::Pdf given a shape described by a functor and a integrator
 *  (algorithm or functor).
 * \param functor shape.
 * \param integrator algorithm or functor.
 * \return a hydra::Pdf instance.
 */
template<typename FUNCTOR>
PdfFromNormalizedFunctor<FUNCTOR> make_pdf_from_normalized_functor( FUNCTOR const& functor)
{

	return PdfFromNormalizedFunctor<FUNCTOR>(functor);
}




template<typename Functor, typename IteratorD, typename ...IteratorW>
class LogLikelihoodFCN< PdfFromNormalizedFunctor<Functor> , IteratorD, IteratorW...>: public FCN<LogLikelihoodFCN< PdfFromNormalizedFunctor<Functor>,IteratorD, IteratorW... >, true >{

public:

	typedef void likelihood_estimator_type;

	LogLikelihoodFCN()=delete;


	/**
	 * @brief LogLikelihoodFCN constructor for non-cached models.
	 *
	 * @param functor hydra::PDF instance.
	 * @param begin  IteratorD pointing to the begin of the dataset.
	 * @param end   IteratorD pointing to the end of the dataset.
	 */
	LogLikelihoodFCN(PdfFromNormalizedFunctor<Functor> const& functor, IteratorD begin, IteratorD end, IteratorW ...wbegin):
		FCN<LogLikelihoodFCN< PdfFromNormalizedFunctor<Functor>, IteratorD, IteratorW...>, true>(functor,begin, end, wbegin...)
		{}

	LogLikelihoodFCN(LogLikelihoodFCN< PdfFromNormalizedFunctor<Functor>, IteratorD, IteratorW...>const& other):
		FCN<LogLikelihoodFCN< PdfFromNormalizedFunctor<Functor>, IteratorD, IteratorW...>, true>(other)
		{}

	LogLikelihoodFCN< PdfFromNormalizedFunctor<Functor>, IteratorD>&
	operator=(LogLikelihoodFCN< PdfFromNormalizedFunctor<Functor>, IteratorD, IteratorW...>const& other)
	{
		if(this==&other) return  *this;
		FCN<LogLikelihoodFCN< PdfFromNormalizedFunctor<Functor>, IteratorD, IteratorW...>, true>::operator=(other);

		return  *this;
	}

	template<size_t M = sizeof...(IteratorW)>
	inline typename std::enable_if<(M==0), double >::type
	Eval( const std::vector<double>& parameters ) const{

		using   hydra::thrust::system::detail::generic::select_system;
		typedef typename hydra::thrust::iterator_system<IteratorD>::type System;
		typedef typename PdfFromNormalizedFunctor<Functor>::functor_type functor_type;
		System system;

		// create iterators
		hydra::thrust::counting_iterator<size_t> first(0);
		hydra::thrust::counting_iterator<size_t> last = first + this->GetDataSize();

		GReal_t final;
		GReal_t init=0;

		if (INFO >= Print::Level()  )
		{
			std::ostringstream stringStream;
			for(size_t i=0; i< parameters.size(); i++){
				stringStream << "Parameter["<< i<<"] :  " << parameters[i]  << "  ";
			}
			HYDRA_LOG(INFO, stringStream.str().c_str() )
		}

		const_cast< LogLikelihoodFCN< PdfFromNormalizedFunctor<Functor>, IteratorD, IteratorW...>* >(this)->GetPDF().SetParameters(parameters);

		auto NLL = detail::LogLikelihood1<functor_type>(this->GetPDF().GetFunctor());

		final = hydra::thrust::transform_reduce(select_system(system),
				this->begin(), this->end(), NLL, init, hydra::thrust::plus<GReal_t>());

		return (GReal_t)this->GetDataSize() -final ;
	}

	template<size_t M = sizeof...(IteratorW)>
	inline typename std::enable_if<(M>0), double >::type
	Eval( const std::vector<double>& parameters ) const{

		using   hydra::thrust::system::detail::generic::select_system;
		typedef typename hydra::thrust::iterator_system<typename FCN<LogLikelihoodFCN< PdfFromNormalizedFunctor<Functor>, IteratorD, IteratorW...>>::iterator>::type System;
		typedef typename PdfFromNormalizedFunctor<Functor>::functor_type functor_type;
		System system;

		// create iterators
		hydra::thrust::counting_iterator<size_t> first(0);
		hydra::thrust::counting_iterator<size_t> last = first + this->GetDataSize();

		GReal_t final;
		GReal_t init=0;

		if (INFO >= Print::Level()  )
		{
			std::ostringstream stringStream;
			for(size_t i=0; i< parameters.size(); i++){
				stringStream << "Parameter["<< i<<"] :  " << parameters[i]  << "  ";
			}
			HYDRA_LOG(INFO, stringStream.str().c_str() )
		}

		const_cast< LogLikelihoodFCN< PdfFromNormalizedFunctor<Functor>, IteratorD, IteratorW...>* >(this)->GetPDF().SetParameters(parameters);


		auto NLL = detail::LogLikelihood2<functor_type>(this->GetPDF().GetFunctor());

		final = hydra::thrust::inner_product(select_system(system), this->begin(), this->end(),this->wbegin(),
				init,hydra::thrust::plus<GReal_t>(),NLL );

		return (GReal_t)this->GetDataSize() -final ;
	}

};


template< typename Functor,  typename Iterator, typename ...Iterators>
inline typename std::enable_if< detail::is_iterator<Iterator>::value && detail::are_iterators<Iterators...>::value,
     LogLikelihoodFCN< PdfFromNormalizedFunctor<Functor>, Iterator , Iterators... > >::type
make_loglikehood_fcn(PdfFromNormalizedFunctor<Functor> const& pdf, Iterator first, Iterator last,  Iterators... weights )
{

	return LogLikelihoodFCN< PdfFromNormalizedFunctor<Functor>, Iterator , Iterators... >(pdf, first, last, weights...);
}



template< typename Functor, typename Iterable, typename ...Iterables>
inline typename std::enable_if< (!detail::is_iterator<Iterable>::value) &&
                                ((sizeof...(Iterables)==0) || !detail::are_iterators<Iterables...>::value) &&
								(!hydra::detail::is_hydra_dense_histogram< typename std::remove_reference<Iterable>::type>::value) &&
								(!hydra::detail::is_hydra_sparse_histogram<typename std::remove_reference<Iterable>::type>::value) &&
								detail::is_iterable<Iterable>::value && detail::are_iterables<Iterables...>::value ,
	              LogLikelihoodFCN< PdfFromNormalizedFunctor<Functor>,
	                   decltype(std::declval<Iterable>().begin()),
                       decltype(std::declval<Iterables>().begin())... >>::type
make_loglikehood_fcn(PdfFromNormalizedFunctor<Functor> const& pdf, Iterable&& points, Iterables&&... weights )
{
	return make_loglikehood_fcn(pdf,
			std::forward<Iterable>(points).begin(),
			std::forward<Iterable>(points).end(),
			std::forward<Iterables>(weights).begin()...);
}



template< typename Functor, typename Histogram>
inline typename std::enable_if<detail::is_hydra_dense_histogram<Histogram>::value ||
                               detail::is_hydra_sparse_histogram<Histogram>::value,
LogLikelihoodFCN< PdfFromNormalizedFunctor<Functor>,
				  decltype(std::declval<const Histogram&>().GetBinsCenters().begin()),
                  decltype( std::declval<const Histogram&>().GetBinsContents().begin())>>::type
make_loglikehood_fcn(PdfFromNormalizedFunctor<Functor> const& pdf, Histogram const& points )
{
	return LogLikelihoodFCN< PdfFromNormalizedFunctor<Functor>,
			  decltype(std::declval<const Histogram&>().GetBinsCenters().begin()),
              decltype( std::declval<const Histogram&>().GetBinsContents().begin())>(pdf, points.GetBinsCenters().begin(),
			points.GetBinsCenters().end(), points.GetBinsContents().begin());
}





} // namespace hydra